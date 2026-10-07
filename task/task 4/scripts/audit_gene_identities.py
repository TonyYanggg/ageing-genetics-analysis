"""Audit NCBI identity/taxon/type for every supplied GenAge identifier."""
from collections import Counter
import json
from collect_oma_pairwise import io, OUT

TAXA={
 'Saccharomyces cerevisiae':{4932,559292},'Caenorhabditis elegans':{6239},
 'Drosophila melanogaster':{7227},'Mus musculus':{10090},
 'Schizosaccharomyces pombe':{4896,284812},'Podospora anserina':{5145,515849,2587412},
 'Caenorhabditis briggsae':{6238},'Danio rerio':{7955},'Mesocricetus auratus':{10036}}


def main():
    raw=io.read_csv(io.SOURCE)
    ids=sorted({r['entrez gene id'] for r in raw if r['entrez gene id']},key=int)
    details={}
    for start in range(0,len(ids),100):
        batch=ids[start:start+100]
        data=io.cached(f'ncbi_identity_{start//100:03d}',
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
            {'db':'gene','id':','.join(batch),'retmode':'json'})
        assert 'result' in data,data
        details.update({k:v for k,v in data['result'].items() if k!='uids'})
        print('NCBI identity',start+len(batch),'/',len(ids),flush=True)
    # ESummary does not expose gene type; obtain it from Datasets for unmapped records.
    unmapped=io.read_csv(OUT/'record_mapping_audit.csv')
    pending=sorted({r['entrez_id'] for r in unmapped if r['entrez_id'] and r['status']!='mapped'},key=int)
    structured={}
    for start in range(0,len(pending),50):
        batch=pending[start:start+50]
        data=io.cached(f'ncbi_unmapped_dataset_{start//50:03d}',
            'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/'+','.join(batch)+'/dataset_report',
            {'page_size':1000})
        assert not data.get('next_page_token'), 'Unexpected pagination'
        for report in data.get('reports',[]):
            gene=report.get('gene',{})
            if gene.get('gene_id'):structured[gene['gene_id']]=gene
        print('Structured unmapped metadata',start+len(batch),'/',len(pending),flush=True)
    rows=[]
    for r in raw:
        d=details.get(r['entrez gene id'],{})
        tax=d.get('organism',{}).get('taxid','')
        state=('missing_id' if not r['entrez gene id'] else 'unresolved_ncbi' if not tax
               else 'taxon_mismatch' if int(tax) not in TAXA[r['organism']]
               else 'discontinued_or_replaced' if d.get('currentid') not in (None,'',0,'0') or d.get('status') not in (None,'',0,'0')
               else 'matched_taxon')
        rows.append(dict(genage_id=r['GenAge ID'],species=r['organism'],symbol=r['symbol'],entrez_id=r['entrez gene id'],
            identity_status=state,ncbi_symbol=d.get('name',''),ncbi_taxon=tax,
            ncbi_organism=d.get('organism',{}).get('scientificname',''),gene_type=structured.get(r['entrez gene id'],{}).get('type','not_queried_or_unavailable'),
            current_id=d.get('currentid',''),ncbi_status=d.get('status',''),
            aliases=d.get('otheraliases',''),description=d.get('description','')))
    io.write_csv('gene_identity_audit.csv',rows,list(rows[0]))
    summary=dict(records=len(rows),identity_states=dict(Counter(r['identity_status'] for r in rows)),
                 gene_types=dict(Counter(str(r['gene_type']) for r in rows)))
    (OUT/'gene_identity_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))


if __name__=='__main__':main()
