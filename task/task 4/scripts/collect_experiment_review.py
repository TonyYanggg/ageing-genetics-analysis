"""Freeze primary-paper metadata/abstracts for explicit, documented review cases."""
from collections import defaultdict
import hashlib
import json
from collect_oma_pairwise import io, OUT

CASE_IDS={'128','137','146','1833','2075','2104','2153','2163','2216','1808','695','2025',
          '578','583','1737'}


def main():
    reference_file=io.ROOT/'task/task 3/results/verified_timeline/entry_references.csv'
    references=io.read_csv(reference_file)
    selected=[r for r in references if r['genage_id'] in CASE_IDS]
    papers={}
    for row in selected:
        pmid=row['pubmed_url'].rstrip('/').split('/')[-1]
        if pmid in papers:continue
        data=io.cached('review_pmid_'+pmid,'https://www.ebi.ac.uk/europepmc/webservices/rest/search',
            {'query':'EXT_ID:'+pmid+' AND SRC:MED','format':'json','resultType':'core'})
        hits=data['resultList']['result']
        assert len(hits)==1 and hits[0]['id']==pmid
        papers[pmid]=hits[0]
        print(pmid,hits[0]['title'],flush=True)
    rows=[]
    for row in selected:
        pmid=row['pubmed_url'].rstrip('/').split('/')[-1]
        paper=papers[pmid]
        rows.append(dict(genage_id=row['genage_id'],symbol=row['symbol'],species=row['organism'],
            pmid=pmid,doi=paper.get('doi',''),title=paper['title'],abstract=paper.get('abstractText',''),
            pubmed_url='https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/',
            catalogue_reference_url=row['source_url'],catalogue_gene_url=row['gene_link'],
            reference_table_sha256=hashlib.sha256(reference_file.read_bytes()).hexdigest(),
            selection_reason=('native S6-kinase intervention comparison' if row['genage_id'] in {'578','583','1737'}
                              else 'sparse-organism coverage or taxonomic-identity discrepancy')))
    io.write_csv('experiment_review_sources.csv',rows,list(rows[0]))


if __name__=='__main__':main()
