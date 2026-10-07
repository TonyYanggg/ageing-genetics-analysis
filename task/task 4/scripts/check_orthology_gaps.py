"""Freeze supplementary source-gap checks; never merge incompatible groups."""
import json
from collect_oma_pairwise import io, OUT


def get(name,path,**params):
    return io.cached('ensembl_'+name,'https://rest.ensembl.org/'+path,
                     {'content-type':'application/json',**params})


def main():
    release=get('release','info/data')
    species=get('fungi_species','info/species',division='EnsemblFungi')
    pod=[s for s in species['species'] if s['name']=='podospora_anserina_s_mat_gca_000226545']
    assert len(pod)==1
    hamster=get('hamster_gene','lookup/id/ENSMAUG00000003914')
    assert hamster['species']=='mesocricetus_auratus' and '101830758' in hamster['description']
    raw=io.read_csv(io.SOURCE)
    supplementary=[]
    for target in ['mus_musculus','danio_rerio']:
        data=get('hamster_'+target,'homology/id/mesocricetus_auratus/ENSMAUG00000003914',
                 type='orthologues',target_species=target,sequence='none')
        for entry in data['data']:
            for h in entry['homologies']:
                assert h['source']['id']=='ENSMAUG00000003914'
                assert h['target']['species']==target
                xrefs=get('target_'+h['target']['id'],'xrefs/id/'+h['target']['id'],external_db='EntrezGene')
                ids={x['primary_id'] for x in xrefs if x['dbname']=='EntrezGene'}
                gs='Mus musculus' if target=='mus_musculus' else 'Danio rerio'
                matches=[r for r in raw if r['organism']==gs and r['entrez gene id'] in ids]
                supplementary.append(dict(source_species='Mesocricetus auratus',source_entrez='101830758',
                    target_species=gs,target_ensembl=h['target']['id'],target_entrez=';'.join(sorted(ids)),
                    relation=h['type'],taxonomy_level=h['taxonomy_level'],
                    genage_ids=';'.join(r['GenAge ID'] for r in matches),
                    qualifies_recurrence=any(r['lifespan effect'] in {'Increase','Decrease','Increase and Decrease'} for r in matches)))
    io.write_csv('supplementary_ensembl_pairs.csv',supplementary,list(supplementary[0]))
    checks=[]
    for gene in ['6194500','802462']:
        hits=get('podospora_'+gene,'xrefs/symbol/'+pod[0]['name']+'/'+gene)
        checks.append(dict(entrez_id=gene,exact_xref_hits=hits))
    identity=io.cached('podospora_orthodb_identity','https://data.orthodb.org/v12/ogdetails',{'id':'515849_0:0021bc'})
    assert any(x['id']=='6194500' for x in identity['data']['entrezgene'])
    assert any(x['id']=='B2B296' for x in identity['data']['uniprot'])
    bridge=get('podospora_uniprot_bridge','xrefs/symbol/'+pod[0]['name']+'/B2B296')
    comparas=get('comparative_databases','info/comparas')
    result=dict(ensembl_release=release,podospora_genome=pod[0],podospora_id_checks=checks,
        podospora_uniprot_bridge_hits=bridge,available_comparative_databases=comparas,
        hamster_exact_identity_confirmed=True,hamster_comparisons=['Mus musculus','Danio rerio'],
        limitations='Ensembl supplements are separate pairs, not merged HOGs. Hamster versus nonvertebrates and Podospora orthology remain unassessed/unresolved; an empty xref result is not biological absence.')
    assert release==get('release_after','info/data')
    (OUT/'source_gap_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':main()
