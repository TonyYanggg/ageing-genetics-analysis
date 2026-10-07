"""Establish observation-to-entry identity directly from every Build 21 page.

Serial, cached public requests; reuse the validated pilot snapshots. LibAge
display labels never determine identity. Reruns use the saved responses.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re
import time
import pandas as pd
import requests
from bs4 import BeautifulSoup
from scrape_genage_references import parse_detail_page, validate_against_download

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parents[1] / 'GenAge Project/models_genes/genage_models.csv'
OUT = ROOT / 'results/identity_verification'


def verify():
    raw = OUT / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(SOURCE)
    session = requests.Session()
    session.headers['User-Agent'] = 'GenAgeTimelineResearch/0.3 (serial cached academic research)'
    observations, manifest, errors = [], [], []
    for n, genage_id in enumerate(genes['GenAge ID'].astype(int), 1):
        url = f'https://genomics.senescence.info/genes/details.php?id={genage_id}'
        path = raw / f'genage_{genage_id}.html'
        pilot = ROOT / f'results/pilot/raw_html/genage_{genage_id}.html'
        if pilot.exists():
            path = pilot
            metadata = dict(url=url, retrieved_utc=None, provenance='pilot snapshot; see pilot report for extraction date')
        elif path.exists():
            metadata = json.loads(path.with_suffix('.meta.json').read_text())
        else:
            for attempt in range(3):
                try:
                    response = session.get(url, timeout=45)
                    response.raise_for_status()
                    break
                except (requests.ConnectionError, requests.Timeout):
                    if attempt == 2:
                        raise
                    time.sleep(10 * (attempt + 1))
            metadata = dict(url=url, retrieved_utc=datetime.now(timezone.utc).isoformat())
            path.write_bytes(response.content)
            path.with_suffix('.meta.json').write_text(json.dumps(metadata, indent=2))
            time.sleep(1)
        body = path.read_bytes()
        digest = hashlib.sha256(body).hexdigest()
        metadata.update(file=str(path.relative_to(ROOT)), sha256=digest)
        manifest.append(metadata)
        try:
            soup = BeautifulSoup(body, 'html.parser')
            heading = soup.find(lambda tag: tag.name == 'h2' and 'Gene information' in tag.get_text())
            match = re.search(r'ID:\s*(\d+)', heading.get_text() if heading else '')
            if not match or int(match.group(1)) != genage_id:
                raise ValueError('Missing or incorrect explicit GenAge ID in page heading')
            gene, rows, _ = parse_detail_page(str(soup), requested_id=genage_id, source_url=url,
                source_sha256=digest, extracted_at_utc=datetime.now(timezone.utc).isoformat())
            validate_against_download(gene, genes)
            if not rows or any(not row['observation_source_id'] for row in rows):
                raise ValueError('Entry has no observations or an observation lacks a source ID')
            observations.extend(rows)
        except ValueError as exc:
            errors.append(dict(genage_id=genage_id, error=str(exc), source_file=str(path.relative_to(ROOT))))
        if n % 25 == 0 or n == len(genes):
            pd.DataFrame(observations).to_csv(OUT / 'verified_observations.csv', index=False)
            pd.DataFrame(errors, columns=['genage_id', 'error', 'source_file']).to_csv(OUT / 'entry_errors.csv', index=False)
            (OUT / 'source_manifest.json').write_text(json.dumps(manifest, indent=2))
            (OUT / 'status.json').write_text(json.dumps(dict(complete=n == len(genes),
                entries_checked=n, total_entries=len(genes), observations=len(observations), errors=len(errors)), indent=2))
            print(f'Checked {n}/{len(genes)} entries; {len(observations)} observations; {len(errors)} errors', flush=True)


if __name__ == '__main__':
    verify()
