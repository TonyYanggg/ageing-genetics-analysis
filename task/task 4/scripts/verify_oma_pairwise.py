"""Independently reconstruct gene pairs from cached API JSON and raw GenAge."""
import csv
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

TASK=Path(__file__).resolve().parents[1]
OUT=TASK/'results/oma_pairwise'
CACHE=OUT/'source_cache'


def csvrows(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))


def main(with_bridges=False):
    target=OUT/'refined' if with_bridges else OUT
    source=TASK.parents[1]/'GenAge Project/models_genes/genage_models.csv'
    status=json.loads((OUT/'collection_status.json').read_text())
    assert status['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
    checked=0
    for path in CACHE.glob('*.meta.json'):
        meta=json.loads(path.read_text())
        payload=path.with_name(path.name.removesuffix('.meta.json') if meta.get('format')=='xml' else path.name.replace('.meta.json','.json'))
        assert hashlib.sha256(payload.read_bytes()).hexdigest()==meta['sha256']
        checked+=1
    taxa={str(r['taxon']):r['species'] for r in csvrows(OUT/'source_species.csv') if r['taxon']}
    records=defaultdict(list)
    for r in csvrows(source):
        if r['entrez gene id']:records[(r['organism'],r['entrez gene id'])].append(r)
    mapping=defaultdict(set)
    for path in CACHE.glob('mapping_*.json'):
        if path.name.endswith('.meta.json'):continue
        for r in json.loads(path.read_text())['results']['bindings']:
            key=(taxa[r['organism']['value'].split('/')[-1]],r['xref']['value'].split('/')[-1])
            if key in records:mapping[r['protein']['value'].split('/')[-1]].add(key)
    if with_bridges:
        direct={key for keys in mapping.values() for key in keys}
        candidates=set()
        for path in CACHE.glob('ncbi_unmapped_dataset_*.json'):
            if path.name.endswith('.meta.json'):continue
            for report in json.loads(path.read_text()).get('reports',[]):
                gene=report.get('gene',{})
                if gene.get('type')!='PROTEIN_CODING':continue
                for key in records:
                    if key[1]!=gene['gene_id'] or key in direct:continue
                    if key[0] not in taxa.values():continue
                    if not gene.get('taxname','').startswith(key[0]):continue
                    refs=list(gene.get('swiss_prot_accessions',[]))+list(gene.get('ensembl_gene_ids',[]))
                    auth=gene.get('nomenclature_authority',{})
                    if key[0]=='Saccharomyces cerevisiae' and gene.get('locus_tag'):
                        refs.append(gene['locus_tag'])
                    if auth.get('authority') in {'FLYBASE','WormBase','SGD','MGI','PomBase'} and auth.get('identifier'):
                        refs.append(auth['identifier'])
                    for ref in refs:candidates.add((key,ref))
        derived={( (r['species'],r['entrez_id']),r['identifier']) for r in csvrows(OUT/'identifier_bridge_candidates.csv')}
        assert derived==candidates
        for key,ref in candidates:
            cache=CACHE/('oma_xref_bridge_'+ref.replace(':','_')+'.json')
            for hit in json.loads(cache.read_text()):
                if hit['xref']!=ref:continue
                if taxa.get(str(hit['genome']['taxon_id']))!=key[0]:continue
                mapping[hit['omaid']].add(key)
    badgenes={g for gs in mapping.values() if len(gs)>1 for g in gs}
    protein_gene={p:next(iter(gs)) for p,gs in mapping.items() if len(gs)==1 and not(gs & badgenes)}
    expected=set()
    expected_relations=defaultdict(set)
    expected_proteins=defaultdict(set)
    source_count=0
    genome_pair_counts=defaultdict(int)
    for path in CACHE.glob('pairs_*.json'):
        if path.name.endswith('.meta.json'):continue
        rows=json.loads(path.read_text())
        source_count+=len(rows)
        genome_pair_counts['_'.join(path.stem.split('_')[1:3])]+=len(rows)
        for r in rows:
            a,b=r['entry_1']['omaid'],r['entry_2']['omaid']
            if a in protein_gene and b in protein_gene:
                key=tuple(sorted((protein_gene[a],protein_gene[b])))
                expected.add(key)
                reverse=protein_gene[a]!=key[0]
                expected_relations[key].add(':'.join(reversed(r['rel_type'].split(':'))) if reverse else r['rel_type'])
                expected_proteins[key].add(b+'|'+a if reverse else a+'|'+b)
    assert len(genome_pair_counts)==21
    for name,n in genome_pair_counts.items():
        meta=json.loads((CACHE/('pairs_'+name+'_001.meta.json')).read_text())
        headers={k.lower():v for k,v in meta['response_headers'].items()}
        assert n==int(headers['x-total-count'])
    assert source_count==status['source_pairs']
    actual=csvrows(target/'all_genage_pairs.csv')
    keys=lambda rows:{tuple(sorted(((r['species1'],r['entrez1']),(r['species2'],r['entrez2'])))) for r in rows}
    assert len(actual)==len(keys(actual)) and keys(actual)==expected
    effects={'Increase','Decrease','Increase and Decrease'}
    eligible={p for p in expected if all(any(r['lifespan effect'] in effects for r in records[g]) for g in p)}
    directional={p for p in eligible if all(any(r['lifespan effect'] in effects and r['longevity influence'] in {'Pro-Longevity','Anti-Longevity'} for r in records[g]) for g in p)}
    assert keys(csvrows(target/'longevity_pairs.csv'))==eligible
    assert keys(csvrows(target/'directional_sensitivity_pairs.csv'))==directional
    for r in actual:
        key=((r['species1'],r['entrez1']),(r['species2'],r['entrez2']))
        assert set(r['pair_rel_types'].split(';'))==expected_relations[key]
        assert set(r['protein_pairs'].split(';'))==expected_proteins[key]
        assert set(r['genage_ids1'].split(';'))=={x['GenAge ID'] for x in records[(r['species1'],r['entrez1'])]}
        assert set(r['genage_ids2'].split(';'))=={x['GenAge ID'] for x in records[(r['species2'],r['entrez2'])]}
    assert len(csvrows(target/'record_mapping_audit.csv'))==2205
    assert len(csvrows(target/'coverage.csv'))==9
    plotted=csvrows(target/'figure_pair_counts.csv')
    assert len(plotted)==81
    for r in plotted:
        if r['status']=='assessed':
            n=sum({g[0] for g in pair}=={r['species1'],r['species2']} for pair in eligible)
            assert int(r['pairs'])==n
        else:assert r['pairs']==''
    unique_mapped=set(protein_gene.values())
    for r in csvrows(target/'coverage.csv'):
        assert int(r['mapped_genes'])==sum(g[0]==r['species'] for g in unique_mapped)
    figure=TASK/('figures/oma_refined/pairwise_evidence.svg' if with_bridges else 'figures/oma_pairwise/pairwise_evidence.svg')
    assert '<text' in figure.read_text(encoding='utf-8')
    report=dict(passed=True,cached_hashes=checked,raw_source_unchanged=True,
        full_genome_pairs=21,source_protein_pairs=source_count,reconstructed_gene_pairs=len(expected),
        reconstructed_longevity_pairs=len(eligible),reconstructed_directional_pairs=len(directional),
        interpretation='Independent computational reconstruction, not experimental or evolutionary ground truth.')
    (target/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--with-bridges',action='store_true')
    main(parser.parse_args().with_bridges)
