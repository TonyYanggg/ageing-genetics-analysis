"""Retrieve three explicitly open-access main articles as JATS XML, not PDFs."""
from datetime import datetime, timezone
import hashlib
import json
import time
import xml.etree.ElementTree as ET
import requests
from collect_oma_pairwise import io, OUT

PMIDS=['26232943','26976705','28891792']


def main():
    manifest=[]
    for pmid in PMIDS:
        paper=json.loads((io.CACHE/('review_pmid_'+pmid+'.json')).read_text())['resultList']['result'][0]
        assert paper['isOpenAccess']=='Y'
        pmc=paper['pmcid']
        path=io.CACHE/('review_fulltext_'+pmc+'.xml')
        meta_path=path.with_suffix('.xml.meta.json')
        url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML'
        if path.exists():
            meta=json.loads(meta_path.read_text())
            assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['sha256']
        else:
            time.sleep(1)
            response=requests.get(url,timeout=45)
            response.raise_for_status()
            root=ET.fromstring(response.content)
            assert root.tag=='article'
            assert any(x.text==pmid for x in root.findall('.//article-id[@pub-id-type="pmid"]'))
            path.write_bytes(response.content)
            meta=dict(url=url,params=None,format='xml',retrieved_utc=datetime.now(timezone.utc).isoformat(),
                      sha256=hashlib.sha256(response.content).hexdigest(),bytes=len(response.content),
                      response_headers=dict(response.headers))
            meta_path.write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
        manifest.append(dict(pmid=pmid,pmcid=pmc,doi=paper['doi'],title=paper['title'],
            source_url=url,format='JATS XML main article',status='open_access_downloaded',
            sha256=meta['sha256'],retrieved_utc=meta['retrieved_utc'],supporting_information='not_downloaded'))
        print(pmid,pmc,len(path.read_bytes()),flush=True)
    io.write_csv('fulltext_manifest.csv',manifest,list(manifest[0]))


if __name__=='__main__':main()
