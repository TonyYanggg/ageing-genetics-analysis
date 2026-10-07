"""Cache the public GenAge-model bibliography and its observation backlinks.

Run from any directory. Serial requests, one-second spacing, resumable cache.
Pedro recommended scraping in the supplied project email; no login is used.
"""
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urljoin, urlparse, parse_qs
import hashlib
import json
import re
import time
import requests
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'bibliography'
BASE = 'https://libage.ageing-map.org'
SOURCE = BASE + '/entries/sources/genage_model/'


def collect():
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / 'raw'
    cache.mkdir(exist_ok=True)
    session = requests.Session()
    session.headers['User-Agent'] = 'GenAgeTimelineResearch/0.2 (serial cached academic research)'
    manifest = []

    def get(url, filename):
        path = cache / filename
        meta = path.with_suffix('.meta.json')
        if not path.exists():
            for attempt in range(3):
                try:
                    response = session.get(url, timeout=45)
                    break
                except (requests.ConnectionError, requests.Timeout):
                    if attempt == 2:
                        raise
                    time.sleep(10 * (attempt + 1))
            response.raise_for_status()
            body = response.content
            BeautifulSoup(body, 'html.parser')
            path.write_bytes(body)
            meta.write_text(json.dumps({'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat()}, indent=2))
            time.sleep(1)
        body = path.read_bytes()
        record = json.loads(meta.read_text())
        record.update(file=str(path.relative_to(ROOT)), sha256=hashlib.sha256(body).hexdigest())
        manifest.append(record)
        return BeautifulSoup(body, 'html.parser'), record

    first, _ = get(SOURCE + '?page=1&s=&sort=title', 'title_index_1.html')
    expected = int(re.search(r'\((\d+)\)\s*$', first.h1.get_text(' ', strip=True)).group(1))
    last = max(int(parse_qs(urlparse(a['href']).query)['page'][0]) for a in first.select('a[href]') if 'page=' in a['href'])
    ids = set()
    for page in range(1, last + 1):
        soup = first if page == 1 else get(SOURCE + f'?page={page}&s=&sort=title', f'title_index_{page}.html')[0]
        ids.update(int(m.group(1)) for a in soup.select('a[href]') if (m := re.search(r'/entries/entry/(\d+)/', a['href'])))
        print(f'Bibliography page {page}/{last}: {len(ids)} unique references', flush=True)
    if len(ids) != expected:
        raise ValueError(f'Bibliography incomplete: {len(ids)} versus advertised {expected}')
    papers, links = [], []
    for n, ref_id in enumerate(sorted(ids), 1):
        url = BASE + f'/entries/entry/{ref_id}/'
        soup, meta = get(url, f'reference_{ref_id}.html')
        citation = soup.select_one('p.highlight-box')
        if citation is None:
            raise ValueError(f'Missing citation: {url}')
        text = citation.get_text(' ', strip=True)
        year = re.search(r'\(((?:18|19|20)\d{2})\)', text)
        pubmed = citation.select_one('a.libage-pubmed')
        papers.append(dict(libage_id=ref_id, publication_year=int(year.group(1)) if year else None,
                           title=citation.b.get_text(' ', strip=True) if citation.b else '', citation=text,
                           pubmed_url=pubmed['href'] if pubmed else '', source_url=url,
                           source_sha256=meta['sha256'], retrieved_utc=meta['retrieved_utc']))
        for a in soup.select('ol li a[href]'):
            if '/genes/details.php' not in a['href']:
                continue
            item_text = a.parent.get_text(' ', strip=True)
            if 'GenAge: Ageing-Associated Genes in Model Organisms' not in item_text:
                continue
            query = parse_qs(urlparse(a['href']).query)
            links.append(dict(libage_id=ref_id, observation_source_id=query.get('lid', [''])[0],
                              gene_link=a['href'], gene_label=a.get_text(' ', strip=True), source_url=url))
        if n % 20 == 0 or n == len(ids):
            print(f'References {n}/{len(ids)}; observation backlinks {len(links)}', flush=True)
        pd.DataFrame(papers).to_csv(OUT / 'references.csv', index=False)
        pd.DataFrame(links).to_csv(OUT / 'observation_backlinks.csv', index=False)
        (OUT / 'source_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    (OUT / 'collection_status.json').write_text(json.dumps(dict(complete=True, references=len(papers),
        advertised_references=expected, backlinks=len(links), completed_utc=datetime.now(timezone.utc).isoformat()), indent=2))


if __name__ == '__main__':
    collect()
