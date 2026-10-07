"""Retrieve model-only membership of recurrent Eukaryota-level groups."""
from collect_model_orthology import OUT, SPECIES, read_csv, write_csv, sparql, cached


def main():
    groups = read_csv(OUT/'recurrent_groups_eukaryota.csv')
    ids = sorted(r['group_id'] for r in groups)
    organisms = ' '.join('<http://purl.orthodb.org/odborganism/'+o+'>' for o in SPECIES.values())
    all_rows=[]
    for n,start in enumerate(range(0,len(ids),20)):
        values=' '.join('<http://purl.orthodb.org/odbgroup/'+g+'>' for g in ids[start:start+20])
        query='''SELECT DISTINCT ?og ?gene ?organism ?name WHERE {
VALUES ?og { %s }
VALUES ?organism { %s }
?gene odb:memberOf ?og; up:organism ?organism.
OPTIONAL { ?gene odb:name ?name. }
}''' % (values,organisms)
        all_rows.extend(sparql(f'recurrent_members_{n:03d}',query))
        print(f'Group membership batch {n+1}/{(len(ids)+19)//20}; rows={len(all_rows)}',flush=True)
    write_csv('recurrent_group_members.csv', [{k:r.get(k,'') for k in ['og','gene','organism','name']}
                                            for r in all_rows], ['og','gene','organism','name'])
    release = cached('release_after_members','https://data.orthodb.org/v12/orthodb_release_id')
    assert release == 'v12.2'


if __name__=='__main__':
    main()
