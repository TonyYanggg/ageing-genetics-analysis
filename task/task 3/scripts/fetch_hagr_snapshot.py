"""Retrieve the lab's pinned public HAGR-R export, retaining URL and hash metadata."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/hagr_r_snapshot'
COMMIT = '4e83b50142ea4b8c6625a7a86b6005ff3f606ed6'


def fetch():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_file = OUT / 'manifest.json'
    saved = {item['file']: item for item in json.loads(manifest_file.read_text())} if manifest_file.exists() else {}
    for name in ['models.csv', 'longevity.csv', 'gene_biblio.csv', 'biblio.csv']:
        path = OUT / name
        if path.exists():
            assert name in saved and hashlib.sha256(path.read_bytes()).hexdigest() == saved[name]['sha256']
            continue
        url = f'https://raw.githubusercontent.com/maglab/HAGR-R/{COMMIT}/data-raw/GenAgeModel/{name}'
        response = requests.get(url, timeout=45)
        response.raise_for_status()
        path.write_bytes(response.content)
        saved[name] = dict(file=name, url=url, commit=COMMIT, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                           sha256=hashlib.sha256(response.content).hexdigest())
        manifest_file.write_text(json.dumps(list(saved.values()), indent=2))
    print('Pinned HAGR-R snapshot present; all saved hashes verified.')


if __name__ == '__main__':
    fetch()
