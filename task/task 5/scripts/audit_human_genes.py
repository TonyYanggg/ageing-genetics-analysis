"""Audit all GenAge human identifiers and annotate LongevityMap membership."""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import time
import requests

TASK=Path(__file__).resolve().parents[1]
ROOT=TASK.parents[1]
GENAGE=ROOT/'GenAge Project/human_genes/genage_human.csv'
LONGEVITY=ROOT/'Other downloads from HAGR/longevity_genes/longevity.csv'
OUT=TASK/'results'
CACHE=OUT/'source_cache'


def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def cached_get(name,url,params):
    CACHE.mkdir(parents=True,exist_ok=True)
    path=CACHE/(name+'.json');meta=CACHE/(name+'.meta.json')
    if path.exists():
        m=json.loads(meta.read_text())
        assert m['url']==url and m['params']==params
        assert hashlib.sha256(path.read_bytes()).hexdigest()==m['sha256']
        return json.loads(path.read_text())
    time.sleep(.5)
    response=requests.get(url,params=params,timeout=60)
    response.raise_for_status()
    data=response.json()
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':')).encode()
    path.write_bytes(payload)
    meta.write_text(json.dumps(dict(url=url,params=params,retrieved_utc=datetime.now(timezone.utc).isoformat(),
        sha256=hashlib.sha256(payload).hexdigest(),bytes=len(payload),response_headers=dict(response.headers)),indent=2)+'\n')
    return data


def main():
    genes=read_csv(GENAGE);longevity=read_csv(LONGEVITY)
    assert len(genes)==len({r['GenAge ID'] for r in genes})==len({r['symbol'] for r in genes})==307
    assert len({r['entrez gene id'] for r in genes})==307
    ids=sorted((r['entrez gene id'] for r in genes),key=int)
    current={}
    for start in range(0,len(ids),100):
        batch=ids[start:start+100]
        url='https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/'+','.join(batch)+'/dataset_report'
        data=cached_get('ncbi_human_genes_%03d'%(start//100),url,{'page_size':1000})
        assert not data.get('next_page_token')
        for report in data.get('reports',[]):
            gene=report.get('gene',{})
            if gene.get('gene_id'):current[gene['gene_id']]=gene
        print('NCBI human genes',min(start+100,len(ids)),'/',len(ids),flush=True)
    audit=[]
    for row in genes:
        gene=current.get(row['entrez gene id'],{})
        tax=str(gene.get('tax_id',''))
        state=('missing_from_current_ncbi' if not gene else 'taxon_mismatch' if tax!='9606'
               else 'missing_ensembl_id' if not gene.get('ensembl_gene_ids') else 'mapped')
        audit.append(dict(genage_id=row['GenAge ID'],genage_symbol=row['symbol'],entrez_id=row['entrez gene id'],
            genage_name=row['name'],genage_uniprot=row['uniprot'],genage_basis=row['why'],identity_status=state,
            current_symbol=gene.get('symbol',''),current_description=gene.get('description',''),current_type=gene.get('type',''),
            current_taxon=tax,ensembl_gene_ids=';'.join(gene.get('ensembl_gene_ids',[])),
            hgnc_id=gene.get('nomenclature_authority',{}).get('identifier',''),
            swiss_prot_accessions=';'.join(gene.get('swiss_prot_accessions',[])),
            synonyms=';'.join(gene.get('synonyms',[]))))
    write_csv(OUT/'human_gene_audit.csv',audit,list(audit[0]))

    # Map study-level LongevityMap gene mentions, preferring catalogue/current symbols;
    # unique current aliases are a fallback and remain explicitly labelled.
    primary=defaultdict(set);aliases=defaultdict(set)
    for row in audit:
        key=row['genage_id']
        primary[row['genage_symbol'].upper()].add(key)
        if row['current_symbol']:primary[row['current_symbol'].upper()].add(key)
        for alias in row['synonyms'].split(';'):
            if alias:aliases[alias.upper()].add(key)
    matches=[];unmatched=[];ambiguous=[]
    for entry in longevity:
        for mention in [x.strip() for x in entry['Gene(s)'].split(',') if x.strip()]:
            symbol=mention.upper();keys=primary.get(symbol,set());route='catalogue_or_current_symbol'
            if not keys:
                keys=aliases.get(symbol,set());route='unique_current_alias'
            if len(keys)==1:
                gene_id=next(iter(keys));a=next(r for r in audit if r['genage_id']==gene_id)
                matches.append(dict(genage_id=gene_id,genage_symbol=a['genage_symbol'],entrez_id=a['entrez_id'],
                    longevitymap_id=entry['id'],association=entry['Association'],population=entry['Population'],
                    variants=entry['Variant(s)'],reported_genes=entry['Gene(s)'],matched_mention=mention,
                    match_route=route,pubmed_id=entry['PubMed']))
            elif not keys:unmatched.append(dict(longevitymap_id=entry['id'],gene_mention=mention))
            else:
                ambiguous.append(dict(longevitymap_id=entry['id'],gene_mention=mention,
                    candidate_genage_ids=';'.join(sorted(keys,key=int))))
    fields=['genage_id','genage_symbol','entrez_id','longevitymap_id','association','population','variants',
            'reported_genes','matched_mention','match_route','pubmed_id']
    write_csv(OUT/'longevitymap_matches.csv',matches,fields)
    write_csv(OUT/'longevitymap_unmatched_mentions.csv',unmatched,['longevitymap_id','gene_mention'])
    write_csv(OUT/'longevitymap_ambiguous_mentions.csv',ambiguous,['longevitymap_id','gene_mention','candidate_genage_ids'])
    summary=dict(genage_rows=len(genes),unique_genage_symbols=len({r['symbol'] for r in genes}),
        unique_entrez_ids=len(ids),identity_states=dict(Counter(r['identity_status'] for r in audit)),
        gene_types=dict(Counter(r['current_type'] for r in audit)),longevitymap_records=len(longevity),
        longevitymap_gene_mentions=sum(bool(x.strip()) for r in longevity for x in r['Gene(s)'].split(',')),
        matched_longevitymap_mentions=len(matches),matched_genage_genes=len({r['genage_id'] for r in matches}),
        unmatched_longevitymap_mentions=len(unmatched),ambiguous_longevitymap_mentions=len(ambiguous),
        input_sha256={'genage_human':hashlib.sha256(GENAGE.read_bytes()).hexdigest(),
                      'longevitymap':hashlib.sha256(LONGEVITY.read_bytes()).hexdigest()},
        limitation='LongevityMap symbol matching is an annotation audit, not proof that every historical symbol is fully resolved.')
    (OUT/'human_gene_audit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
