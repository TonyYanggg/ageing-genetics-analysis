"""Attach GenAge evidence to explicit OMA pairs, never infer transitive groups."""
from collections import Counter, defaultdict
import argparse
import hashlib
import itertools
import json
from collect_oma_pairwise import OUT, SPECIES, io

EFFECTS = {'Increase','Decrease','Increase and Decrease'}
DIRECTIONS = {'Pro-Longevity','Anti-Longevity'}


def joined(values):
    return ';'.join(sorted(set(str(v) for v in values if v)))


def main(with_bridges=False):
    destination=OUT/'refined' if with_bridges else OUT
    io.OUT=destination
    status=json.loads((OUT/'collection_status.json').read_text())
    assert status['complete']
    assert hashlib.sha256(io.SOURCE.read_bytes()).hexdigest() == status['source_sha256']
    raw=io.read_csv(io.SOURCE)
    records=defaultdict(list)
    for r in raw:
        if r['entrez gene id']:
            records[(r['organism'],r['entrez gene id'])].append(r)
    taxon_species={str(t):s for s,(c,t) in SPECIES.items() if c}
    mapping=defaultdict(set)
    inverse=defaultdict(set)
    mismatch=[]
    for r in io.read_csv(OUT/'identifier_mappings.csv'):
        species=taxon_species[r['organism'].rsplit('/',1)[-1]]
        key=(species,r['xref'].rsplit('/',1)[-1])
        if key not in records:
            mismatch.append(r)
            continue
        protein=r['protein'].rsplit('/',1)[-1]
        assert protein.startswith(SPECIES[species][0])
        mapping[key].add(protein)
        inverse[protein].add(key)
    ambiguous={p for p,keys in inverse.items() if len(keys)>1}
    direct_keys=set(mapping)
    bridge_keys=set()
    if with_bridges:
        for r in io.read_csv(OUT/'identifier_bridge_results.csv'):
            key=(r['species'],r['entrez_id'])
            if not r['oma_proteins']:continue
            assert key in records and key not in direct_keys
            for protein in r['oma_proteins'].split(';'):
                assert protein.startswith(SPECIES[key[0]][0])
                mapping[key].add(protein)
                inverse[protein].add(key)
            bridge_keys.add(key)
        ambiguous={p for p,keys in inverse.items() if len(keys)>1}
    accepted={k:ps for k,ps in mapping.items() if not (ps & ambiguous)}
    audit=[]
    for r in raw:
        key=(r['organism'],r['entrez gene id'])
        state=('missing_id' if not key[1] else 'source_species_unavailable' if not SPECIES[key[0]][0]
               else 'ambiguous_identifier' if mapping[key] & ambiguous
               else 'mapped' if accepted.get(key) else 'no_exact_xref')
        audit.append(dict(genage_id=r['GenAge ID'],species=key[0],entrez_id=key[1],symbol=r['symbol'],
            status=state,mapping_route='stable_identifier_bridge' if key in bridge_keys else 'direct_NCBI_xref' if key in direct_keys else '',
            oma_proteins=joined(mapping[key]),lifespan_effect=r['lifespan effect'],influence=r['longevity influence']))
    io.write_csv('record_mapping_audit.csv',audit,list(audit[0]))
    io.write_csv('taxonomy_mismatches.csv',mismatch,['xref','protein','organism'])
    proteins={p:k for k,ps in accepted.items() for p in ps}
    support=defaultdict(list)
    for r in io.read_csv(OUT/'source_pairs.csv'):
        if r['protein1'] not in proteins or r['protein2'] not in proteins:
            continue
        a,b=proteins[r['protein1']],proteins[r['protein2']]
        assert a[0]!=b[0]
        key=tuple(sorted((a,b)))
        oriented=dict(r)
        if a!=key[0]:
            oriented['protein1'],oriented['protein2']=r['protein2'],r['protein1']
            oriented['rel_type']=':'.join(reversed(r['rel_type'].split(':')))
        support[key].append(oriented)
    pair_rows=[]
    for (a,b),rows in sorted(support.items()):
        aa,bb=records[a],records[b]
        ea=[r for r in aa if r['lifespan effect'] in EFFECTS]
        eb=[r for r in bb if r['lifespan effect'] in EFFECTS]
        da={r['longevity influence'] for r in ea if r['longevity influence'] in DIRECTIONS}
        db={r['longevity influence'] for r in eb if r['longevity influence'] in DIRECTIONS}
        unresolved=any(r['longevity influence'] not in DIRECTIONS for r in ea+eb)
        agreement=('not_both_lifespan' if not(ea and eb) else 'discordant_or_mixed' if len(da|db)>1
                   else 'unresolved' if unresolved or not(da and db) else next(iter(da)))
        pair_rows.append(dict(species1=a[0],entrez1=a[1],symbols1=joined(r['symbol'] for r in aa),
            genage_ids1=joined(r['GenAge ID'] for r in aa),species2=b[0],entrez2=b[1],
            symbols2=joined(r['symbol'] for r in bb),genage_ids2=joined(r['GenAge ID'] for r in bb),
            both_lifespan=bool(ea and eb),both_directional=bool(da and db),direction_agreement=agreement,
            influences1=joined(r['longevity influence'] for r in aa),influences2=joined(r['longevity influence'] for r in bb),
            pair_rel_types=joined(r['rel_type'] for r in rows),protein_pairs=joined(r['protein1']+'|'+r['protein2'] for r in rows),
            source='OMA All.May2026',unit='explicit_gene_pair'))
    io.write_csv('all_genage_pairs.csv',pair_rows,list(pair_rows[0]))
    eligible=[r for r in pair_rows if r['both_lifespan']]
    io.write_csv('longevity_pairs.csv',eligible,list(pair_rows[0]))
    io.write_csv('directional_sensitivity_pairs.csv',[r for r in eligible if r['both_directional']],list(pair_rows[0]))
    coverage=[]
    for species in SPECIES:
        rs=[r for r in audit if r['species']==species]
        ids={r['entrez_id'] for r in rs if r['entrez_id']}
        mapped={r['entrez_id'] for r in rs if r['status']=='mapped'}
        coverage.append(dict(species=species,records=len(rs),identified_genes=len(ids),mapped_genes=len(mapped),
            missing_id_records=sum(not r['entrez_id'] for r in rs),source_status='represented' if SPECIES[species][0] else 'unavailable',
            mapped_fraction=len(mapped)/len(ids),mapping_states=json.dumps(dict(Counter(r['status'] for r in rs)),sort_keys=True)))
    io.write_csv('coverage.csv',coverage,list(coverage[0]))
    # Compare every old recurrent group's cross-species gene pair, without making new groups.
    old=io.read_csv(io.TASK/'results/model_orthology/record_mapping_audit.csv')
    groupkeys=defaultdict(set)
    for r in old:
        if r['mapping_status']=='mapped' and r['qualifies_lifespan']=='True':
            groupkeys[r['eukaryota_groups']].add((r['species'],r['entrez_id']))
    comparison=[]
    for group,keys in sorted(groupkeys.items()):
        if len({k[0] for k in keys})<2:
            continue
        counts=Counter()
        for a,b in itertools.combinations(sorted(keys),2):
            if a[0]==b[0]:continue
            category=('unassessable_mapping' if a not in accepted or b not in accepted
                      else 'explicit_OMA_pair' if (a,b) in support else 'no_OMA_pair_call')
            counts[category]+=1
        comparison.append(dict(orthodb_group=group,explicit_OMA_pairs=counts['explicit_OMA_pair'],
            no_OMA_pair_call=counts['no_OMA_pair_call'],unassessable_mapping=counts['unassessable_mapping']))
    io.write_csv('orthodb_group_crosscheck.csv',comparison,list(comparison[0]))
    assert not any({r['symbols1'],r['symbols2']}=={'Coq2','HSP104'} for r in pair_rows)
    summary=dict(source='OMA All.May2026',unit='unique species+Entrez gene pair',
        mapping_mode='direct plus stable-ID bridges' if with_bridges else 'direct NCBI xrefs only',
        bridge_candidate_genes=len(bridge_keys),accepted_bridge_genes=sum(k in accepted for k in bridge_keys),
        raw_records=len(raw),inventory_species=len(SPECIES),represented_species=7,
        mapped_genes=len(accepted),ambiguous_proteins=len(ambiguous),taxonomy_mismatches=len(mismatch),
        all_genage_pairs=len(pair_rows),longevity_pairs=len(eligible),
        directional_sensitivity_pairs=sum(r['both_directional'] for r in eligible),
        agreement=dict(Counter(r['direction_agreement'] for r in eligible)),
        limitation='Pairs may share genes; do not interpret pair counts as independent conserved genes or groups. Missing calls are not proof of non-orthology.')
    (destination/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--with-bridges',action='store_true')
    main(parser.parse_args().with_bridges)
