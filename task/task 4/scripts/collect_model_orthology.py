"""Cache public OrthoDB 12.2 model-only mappings; never mutate shared inputs.

Batch SPARQL follows the source's published xref and membership examples.
Requests are serial, at least one second apart, with bounded retries.
"""
from pathlib import Path
from collections import Counter, defaultdict
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import time
import requests

TASK = Path(__file__).resolve().parents[1]
ROOT = TASK.parents[1]
SOURCE = ROOT / 'GenAge Project/models_genes/genage_models.csv'
CACHE = TASK / 'results/model_orthology/source_cache'
OUT = TASK / 'results/model_orthology'
SPECIES = {
    'Saccharomyces cerevisiae': '559292_0',
    'Caenorhabditis elegans': '6239_0',
    'Drosophila melanogaster': '7227_0',
    'Mus musculus': '10090_1',
    'Schizosaccharomyces pombe': '4896_0',
    'Podospora anserina': '515849_0',
    'Caenorhabditis briggsae': '6238_0',
    'Danio rerio': '7955_0',
    'Mesocricetus auratus': '10036_0',
}
PREFIX = '''PREFIX odb: <http://purl.orthodb.org/>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
'''


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows, fields):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def cached(name, url, params=None):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / (name + '.json')
    meta_path = CACHE / (name + '.meta.json')
    if path.exists():
        meta = json.loads(meta_path.read_text(encoding='utf-8'))
        assert meta['params'] == params and meta['url'] == url
        assert hashlib.sha256(path.read_bytes()).hexdigest() == meta['sha256']
        return json.loads(path.read_text(encoding='utf-8'))
    for attempt in range(3):
        time.sleep(1 if attempt == 0 else 5 * attempt)
        try:
            response = requests.get(url, params=params, timeout=90)
            response.raise_for_status()
            # Virtuoso can return HTTP 200 for partial timed-out results.
            if any('sql' in key.lower() and val not in ('00000', '')
                   for key, val in response.headers.items()):
                raise RuntimeError('SPARQL server reported incomplete results')
            data = response.json()
            if isinstance(data, dict) and data.get('status') == 'error':
                raise RuntimeError(str(data))
            if 'sparql' in url and 'results' not in data:
                raise RuntimeError('Missing SPARQL results')
            break
        except (requests.RequestException, ValueError, RuntimeError):
            if attempt == 2:
                raise
    path.write_bytes(response.content)
    meta = dict(url=url, params=params, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                sha256=hashlib.sha256(response.content).hexdigest(), bytes=len(response.content),
                response_headers=dict(response.headers))
    meta_path.write_text(json.dumps(meta, indent=2)+'\n', encoding='utf-8')
    return data


def sparql(name, query):
    data = cached(name, 'https://sparql.orthodb.org/sparql',
                  {'query': PREFIX + query, 'format': 'application/sparql-results+json'})
    return [{k: v['value'] for k,v in row.items()} for row in data['results']['bindings']]


def inventory():
    rows = read_csv(SOURCE)
    assert len(rows) == 2205 and set(r['organism'] for r in rows) == set(SPECIES)
    summaries = []
    for species in SPECIES:
        records = [r for r in rows if r['organism'] == species]
        ids = {r['entrez gene id'] for r in records if r['entrez gene id']}
        summaries.append(dict(species=species, records=len(records), identified_genes=len(ids),
            missing_id_records=sum(not r['entrez gene id'] for r in records),
            proposed_orthodb_organism=SPECIES[species],
            annotation_counts=json.dumps(dict(Counter(r['longevity influence'] for r in records)), sort_keys=True)))
    write_csv('inventory.csv', summaries, list(summaries[0]))
    write_csv('genage_records.csv', rows, list(rows[0]))
    return rows


def collect_mapping(max_batches=None):
    rows = inventory()
    release = cached('release', 'https://data.orthodb.org/v12/orthodb_release_id')
    assert release == 'v12.2', f'Unexpected release: {release}'
    taxons = cached('species', 'https://data.orthodb.org/v12/species', {'clade': '2759'})['data']
    selected = {t['taxon_id']: t for t in taxons if t['taxon_id'] in SPECIES.values()}
    assert set(selected) == set(SPECIES.values()), 'Missing representative assembly'
    write_csv('source_species.csv', [dict(species=s, orthodb_id=o, source_name=selected[o]['sciname'],
        assembly=selected[o]['organism_id'], definition=json.dumps(selected[o]['description'],sort_keys=True))
        for s,o in SPECIES.items()], ['species','orthodb_id','source_name','assembly','definition'])
    identifiers = sorted({r['entrez gene id'] for r in rows if r['entrez gene id']}, key=int)
    organisms = ' '.join('<http://purl.orthodb.org/odborganism/'+o+'>' for o in SPECIES.values())
    mapped = []
    batches = [identifiers[i:i+50] for i in range(0,len(identifiers),50)]
    for number, batch in enumerate(batches):
        if max_batches is not None and number >= max_batches:
            break
        values = ' '.join('<http://www.ncbi.nlm.nih.gov/gene/'+i+'>' for i in batch)
        query = '''SELECT DISTINCT ?xref ?gene ?name ?organism ?og ?level ?level_name ?og_name WHERE {
VALUES ?xref { %s }
VALUES ?organism { %s }
?gene odb:xref/odb:xrefResource ?xref; up:organism ?organism.
OPTIONAL { ?gene odb:name ?name. }
OPTIONAL { ?gene odb:memberOf ?og. ?og odb:ogBuiltAt ?level.
OPTIONAL { ?level up:scientificName ?level_name. }
OPTIONAL { ?og odb:name ?og_name. } }
}''' % (values, organisms)
        mapped.extend(sparql(f'mapping_{number:03d}', query))
        print(f'Mapping batches: {number+1}/{len(batches)}; membership rows: {len(mapped)}', flush=True)
    complete = max_batches is None or max_batches >= len(batches)
    fields = ['xref','gene','name','organism','og','level','level_name','og_name']
    write_csv('mapping_memberships.csv' if complete else 'pilot_mapping_memberships.csv',
              [{k:row.get(k,'') for k in fields} for row in mapped], fields)
    status = dict(source_release=release, complete=complete, source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        input_records=len(rows), queried_identifiers=min(len(identifiers), (max_batches or len(batches))*50),
        unique_identifiers=len(identifiers), batches=len(batches), memberships=len(mapped))
    (OUT / ('collection_status.json' if complete else 'pilot_collection_status.json')).write_text(
        json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(status,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-batches', type=int)
    parser.add_argument('--inventory-only', action='store_true')
    args = parser.parse_args()
    if args.inventory_only:
        inventory()
    else:
        collect_mapping(args.max_batches)
