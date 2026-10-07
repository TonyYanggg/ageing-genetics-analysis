"""Freeze OMA model-only exact-ID mappings and explicit genome-pair orthology."""
import argparse
import hashlib
import itertools
import json
import math
import collect_model_orthology as io

OUT = io.TASK/'results/oma_pairwise'
io.OUT = OUT
io.CACHE = OUT/'source_cache'
SPECIES = {
    'Saccharomyces cerevisiae': ('YEAST', 559292),
    'Caenorhabditis elegans': ('CAEEL', 6239),
    'Drosophila melanogaster': ('DROME', 7227),
    'Mus musculus': ('MOUSE', 10090),
    'Schizosaccharomyces pombe': ('SCHPO', 284812),
    'Podospora anserina': (None, None),
    'Caenorhabditis briggsae': ('CAEBR', 6238),
    'Danio rerio': ('DANRE', 7955),
    'Mesocricetus auratus': (None, None),
}


def pages(name, url):
    """Check all pages against server total, avoiding silent truncation."""
    result = []
    page = 1
    while True:
        key = f'{name}_{page:03d}'
        rows = io.cached(key, url, {'per_page': 10000, 'page': page})
        meta = json.loads((io.CACHE/(key+'.meta.json')).read_text())
        headers = {k.lower(): v for k,v in meta['response_headers'].items()}
        total = int(headers['x-total-count'])
        result.extend(rows)
        print(name, page, len(result), '/', total, flush=True)
        if 'rel="next"' not in headers.get('link', ''):
            assert len(result) == total, (name, len(result), total)
            return result
        assert rows
        page += 1


def main(only_mapping=False):
    version = io.cached('version_before', 'https://omabrowser.org/api/version/')
    assert version['oma_version'] == 'All.May2026', version
    genomes = pages('genomes', 'https://omabrowser.org/api/genome/')
    index = {g['code']: g for g in genomes}
    coverage = []
    for species, (code, taxon) in SPECIES.items():
        if code:
            assert index[code]['taxon_id'] == taxon
            detail = io.cached('genome_'+code, index[code]['genome_url'])
        else:
            assert not any(species.lower() in g['species'].lower() for g in genomes)
            missing_taxa=({5145,515849,2587412} if species=='Podospora anserina' else {10036})
            assert not any(g['taxon_id'] in missing_taxa for g in genomes)
            aliases=['anserina','triangularia','podospora'] if species=='Podospora anserina' else ['mesocricetus','auratus']
            assert not any(any(a in g['species'].lower() for a in aliases) for g in genomes)
            detail = {}
        coverage.append(dict(species=species, code=code or '', taxon=taxon or '',
                             status='represented' if code else 'not_in_OMA_release',
                             source_detail=json.dumps(detail, sort_keys=True)))
    io.write_csv('source_species.csv', coverage, list(coverage[0]))
    raw = io.read_csv(io.SOURCE)
    ids = sorted({r['entrez gene id'] for r in raw if r['entrez gene id']}, key=int)
    mappings = []
    for batch in range(math.ceil(len(ids)/50)):
        values = ' '.join('<https://www.ncbi.nlm.nih.gov/gene/'+i+'>' for i in ids[batch*50:(batch+1)*50])
        organisms = ' '.join('<https://omabrowser.org/oma/genome/'+str(t)+'>' for c,t in SPECIES.values() if c)
        query = '''SELECT DISTINCT ?xref ?protein ?organism WHERE {
VALUES ?xref { %s } VALUES ?organism { %s }
?gene <http://purl.org/lscr#xrefNCBIGene> ?xref.
?protein <http://semanticscience.org/resource/SIO_010079> ?gene;
 <http://purl.org/net/orth#organism> ?organism.
}''' % (values, organisms)
        data = io.cached(f'mapping_{batch:03d}', 'https://sparql.omabrowser.org/sparql',
                         {'query': query, 'format': 'application/sparql-results+json'})
        mappings.extend({k: v['value'] for k,v in row.items()} for row in data['results']['bindings'])
        print('mapping',batch,len(mappings),flush=True)
    io.write_csv('identifier_mappings.csv', mappings, ['xref','protein','organism'])
    if only_mapping:
        return
    pairs = []
    for first, second in itertools.combinations([c for c,t in SPECIES.values() if c], 2):
        rows = pages('pairs_'+first+'_'+second, f'https://omabrowser.org/api/pairs/{first}/{second}/')
        for row in rows:
            a,b = row['entry_1'],row['entry_2']
            assert {a['species']['code'],b['species']['code']} == {first,second}
            pairs.append(dict(protein1=a['omaid'],protein2=b['omaid'],
                species1=a['species']['code'],species2=b['species']['code'],
                canonical1=a['canonicalid'],canonical2=b['canonicalid'],
                main_isoform1=a['is_main_isoform'],main_isoform2=b['is_main_isoform'],
                hog1=a['oma_hog_id'],hog2=b['oma_hog_id'],
                rel_type=row['rel_type'],distance=row['distance'],score=row['score']))
    io.write_csv('source_pairs.csv', pairs, list(pairs[0]))
    after = io.cached('version_after', 'https://omabrowser.org/api/version/')
    assert after == version
    status = dict(complete=True, source_version=version, source_sha256=hashlib.sha256(io.SOURCE.read_bytes()).hexdigest(),
                  mapping_rows=len(mappings), genome_pairs=21, source_pairs=len(pairs),
                  represented_species=7, unresolved_source_species=['Podospora anserina','Mesocricetus auratus'])
    (OUT/'collection_status.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--only-mapping',action='store_true')
    main(parser.parse_args().only_mapping)
