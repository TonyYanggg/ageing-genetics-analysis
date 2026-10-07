"""Analyse model-only OrthoDB memberships at explicit taxonomic levels."""
from collections import Counter, defaultdict
import hashlib
import json
from collect_model_orthology import OUT, SOURCE, SPECIES, read_csv, write_csv

EFFECTS = {'Increase', 'Decrease', 'Increase and Decrease'}
DIRECTIONS = {'Pro-Longevity', 'Anti-Longevity'}


def uniq(values):
    return ';'.join(sorted(set(str(v) for v in values if v)))


def main():
    status = json.loads((OUT/'collection_status.json').read_text())
    assert status['complete']
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == status['source_sha256']
    records = read_csv(SOURCE)
    memberships = read_csv(OUT/'mapping_memberships.csv')
    by_key = defaultdict(list)
    for row in records:
        if row['entrez gene id']:
            by_key[(SPECIES[row['organism']], row['entrez gene id'])].append(row)
    group_info = {}
    mapping = defaultdict(lambda: defaultdict(set))
    mapped_genes = defaultdict(set)
    mismatches = []
    for row in memberships:
        key = (row['organism'].rsplit('/',1)[-1], row['xref'].rsplit('/',1)[-1])
        if key not in by_key:
            mismatches.append(row)
            continue
        mapped_genes[key].add(row['gene'])
        if not row['og']:
            continue
        og = row['og'].rsplit('/',1)[-1]
        level = row['level'].rsplit('/',1)[-1]
        group_info[og] = dict(group_id=og,level=level,level_name=row['level_name'],group_name=row['og_name'])
        mapping[key][level].add(og)
    write_csv('taxonomy_mismatches.csv',mismatches,list(memberships[0]))
    ambiguous = []
    accepted = defaultdict(set)
    for key, levels in mapping.items():
        for level, groups in levels.items():
            if len(groups) != 1:
                ambiguous.append(dict(organism=key[0],entrez_id=key[1],level=level,groups=uniq(groups)))
            else:
                accepted[next(iter(groups))].add(key)
    write_csv('ambiguous_group_mappings.csv', ambiguous, ['organism','entrez_id','level','groups'])
    identity = []
    for row in records:
        key = (SPECIES[row['organism']],row['entrez gene id'])
        groups = mapping[key].get('2759',set())
        state = ('missing_entrez_id' if not key[1] else 'no_exact_species_xref' if not mapped_genes[key]
                 else 'no_eukaryota_group' if not groups else 'ambiguous_eukaryota_groups' if len(groups)>1
                 else 'mapped')
        identity.append(dict(genage_id=row['GenAge ID'],species=row['organism'],symbol=row['symbol'],
            entrez_id=key[1],lifespan_effect=row['lifespan effect'],longevity_influence=row['longevity influence'],
            qualifies_lifespan=row['lifespan effect'] in EFFECTS,mapping_status=state,
            orthodb_gene_ids=uniq(mapped_genes[key]),eukaryota_groups=uniq(groups)))
    write_csv('record_mapping_audit.csv',identity,list(identity[0]))
    coverage=[]
    for species,org in SPECIES.items():
        rows=[r for r in identity if r['species']==species]
        ids={r['entrez_id'] for r in rows if r['entrez_id']}
        mapped={r['entrez_id'] for r in rows if r['mapping_status']=='mapped'}
        coverage.append(dict(species=species,records=len(rows),unique_identified_genes=len(ids),
            missing_id_records=sum(not r['entrez_id'] for r in rows),mapped_genes=len(mapped),
            exact_xref_genes=sum(bool(mapped_genes[(org,i)]) for i in ids),
            unresolved_genes=len(ids-mapped),mapped_fraction=len(mapped)/len(ids) if ids else '',
            lifespan_evidence_records=sum(r['qualifies_lifespan'] for r in rows)))
    write_csv('coverage.csv',coverage,list(coverage[0]))
    summaries=[]
    cells=[]
    evidence=[]
    for og,keys in sorted(accepted.items()):
        info=group_info[og]
        n_effect=n_direction=n_catalogue=0
        all_signs=set()
        unresolved=False
        for species,org in SPECIES.items():
            genes=[key for key in keys if key[0]==org]
            rows=[r for key in genes for r in by_key[key]]
            effect=[r for r in rows if r['lifespan effect'] in EFFECTS]
            signs={r['longevity influence'] for r in effect if r['longevity influence'] in DIRECTIONS}
            state=('Mixed' if len(signs)==2 else 'Pro' if signs=={'Pro-Longevity'} else
                   'Anti' if signs=={'Anti-Longevity'} else 'Other' if effect else
                   'Catalogue only' if rows else 'No matched GenAge evidence')
            n_effect+=bool(effect)
            n_direction+=bool(signs)
            n_catalogue+=bool(rows)
            all_signs.update(signs)
            unresolved |= any(r['longevity influence'] not in DIRECTIONS for r in effect)
            cells.append(dict(**info,species=species,state=state,model_symbols=uniq(r['symbol'] for r in rows),
                model_entrez_ids=uniq(r['entrez gene id'] for r in rows),
                genage_ids=uniq(r['GenAge ID'] for r in rows),lifespan_genage_ids=uniq(r['GenAge ID'] for r in effect),
                influence_categories=uniq(r['longevity influence'] for r in rows),
                lifespan_effects=uniq(r['lifespan effect'] for r in rows),
                matched_gene_count=len(genes),qualifying_record_count=len(effect)))
            for r in rows:
                evidence.append(dict(group_id=og,level=info['level'],species=species,
                    genage_id=r['GenAge ID'],symbol=r['symbol'],entrez_id=r['entrez gene id'],
                    lifespan_effect=r['lifespan effect'],longevity_influence=r['longevity influence']))
        agreement=('Discordant' if len(all_signs)==2 else 'Unresolved' if unresolved or n_direction < n_effect
                   else 'Pro' if all_signs=={'Pro-Longevity'} else 'Anti' if all_signs else 'Unresolved')
        summaries.append(dict(**info,catalogue_species=n_catalogue,lifespan_species=n_effect,
            directional_species=n_direction,direction_agreement=agreement,
            matched_model_genes=len(keys),model_symbols=uniq(r['symbol'] for key in keys for r in by_key[key]),
            recurrent=n_effect>=2))
    summaries.sort(key=lambda r:(-r['lifespan_species'],r['matched_model_genes'],r['group_id']))
    write_csv('groups_all_levels.csv',summaries,list(summaries[0]))
    primary=[r for r in summaries if r['level']=='2759']
    recurrent=[r for r in primary if r['recurrent']]
    write_csv('groups_eukaryota.csv',primary,list(summaries[0]))
    write_csv('recurrent_groups_eukaryota.csv',recurrent,list(summaries[0]))
    write_csv('group_species_cells.csv',cells,list(cells[0]))
    write_csv('group_record_evidence.csv',evidence,list(evidence[0]))
    summary=dict(records=len(records),identified_species_genes=len(by_key),
        mapping_memberships=len(memberships),taxonomy_mismatches=len(mismatches),
        ambiguous_gene_level_assignments=len(ambiguous),
        eukaryota_groups=len(primary),recurrent_eukaryota_groups=len(recurrent),
        recurrence_distribution=dict(Counter(r['lifespan_species'] for r in recurrent)),
        agreement_distribution=dict(Counter(r['direction_agreement'] for r in recurrent)),
        mapped_records=sum(r['mapping_status']=='mapped' for r in identity),
        mapping_status_counts=dict(Counter(r['mapping_status'] for r in identity)))
    (OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
