"""Check observation IDs absent, unassigned, or hidden in the lab's bulk export."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from scrape_genage_references import parse_detail_page, validate_against_download

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/snapshot_gap_checks'


def resolve():
    raw = OUT / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT.parents[1] / 'GenAge Project/models_genes/genage_models.csv')
    obs = pd.read_csv(ROOT / 'results/hagr_r_snapshot/longevity.csv')
    links = pd.read_csv(ROOT / 'results/bibliography/observation_backlinks.csv')
    candidates = links.merge(obs, left_on='observation_source_id', right_on='id', how='left')
    candidates = candidates[candidates.gene_id.isna() | candidates.hidden.eq(1)]
    session = requests.Session()
    session.headers['User-Agent'] = 'GenAgeTimelineResearch/0.3 (serial cached academic research)'
    results = []
    source_ids = set(candidates.observation_source_id.astype(int))
    corrections = ROOT / 'results/verified_timeline/corrected_label_assignments.csv'
    if corrections.exists():
        source_ids.update(pd.read_csv(corrections).observation_source_id.astype(int))
    for source_id in sorted(source_ids):
        url = f'https://genomics.senescence.info/genes/details.php?lid={source_id}'
        path = raw / f'lid_{source_id}.html'
        if not path.exists():
            response = session.get(url, timeout=45)
            path.write_bytes(response.content)
            path.with_suffix('.meta.json').write_text(json.dumps(dict(url=url, status_code=response.status_code,
                retrieved_utc=datetime.now(timezone.utc).isoformat()), indent=2))
            time.sleep(1)
        meta = json.loads(path.with_suffix('.meta.json').read_text())
        body = path.read_bytes()
        record = dict(observation_source_id=source_id, genage_id=None, observation_present=False,
            status_code=meta['status_code'], source_url=url, source_file=str(path.relative_to(ROOT)),
            source_sha256=hashlib.sha256(body).hexdigest(), error='')
        soup = BeautifulSoup(body, 'html.parser')
        heading = soup.find(lambda tag: tag.name == 'h2' and 'Gene information' in tag.get_text())
        match = re.search(r'ID:\s*(\d+)', heading.get_text() if heading else '')
        if meta['status_code'] != 200 or not match:
            record['error'] = 'No usable gene entry at backlink destination'
        else:
            gene, rows, _ = parse_detail_page(str(soup), requested_id=int(match.group(1)), source_url=url,
                source_sha256=record['source_sha256'], extracted_at_utc=datetime.now(timezone.utc).isoformat())
            validate_against_download(gene, genes)
            record['genage_id'] = int(gene['genage_id'])
            record['observation_present'] = str(source_id) in {row['observation_source_id'] for row in rows}
            record['symbol'] = gene['symbol']
            if not record['observation_present']:
                record['error'] = 'Destination entry does not display this observation ID'
        results.append(record)
        print(f"Observation {source_id}: entry {record['genage_id']}, present={record['observation_present']}, HTTP {record['status_code']}", flush=True)
        pd.DataFrame(results).to_csv(OUT / 'gap_checks.csv', index=False)


if __name__ == '__main__':
    resolve()
