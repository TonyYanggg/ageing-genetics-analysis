"""Narrow integrity/reconciliation checks, not biological orthology validation."""
import hashlib
import json
from collections import Counter
from collect_model_orthology import CACHE, OUT, SOURCE, SPECIES, read_csv


def main():
    checked = 0
    for metadata in CACHE.glob('*.meta.json'):
        meta = json.loads(metadata.read_text(encoding='utf-8'))
        payload = metadata.with_name(metadata.name.replace('.meta.json', '.json'))
        assert hashlib.sha256(payload.read_bytes()).hexdigest() == meta['sha256'], payload
        checked += 1
    status = json.loads((OUT/'collection_status.json').read_text())
    assert status['complete'] and status['batches'] == 44
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == status['source_sha256']
    raw = read_csv(SOURCE)
    audit = read_csv(OUT/'record_mapping_audit.csv')
    assert len(raw) == len(audit) == 2205
    assert Counter(r['organism'] for r in raw) == Counter(r['species'] for r in audit)
    coverage = read_csv(OUT/'coverage.csv')
    assert {r['species'] for r in coverage} == set(SPECIES)
    for row in coverage:
        keys = {r['entrez_id'] for r in audit if r['species'] == row['species'] and r['mapping_status'] == 'mapped'}
        assert len(keys) == int(row['mapped_genes'])
    summary = json.loads((OUT/'analysis_summary.json').read_text())
    assert summary['mapped_records'] == sum(r['mapping_status'] == 'mapped' for r in audit)
    assert len(read_csv(OUT/'recurrent_groups_eukaryota.csv')) == 181
    validation = json.loads((OUT/'source_validation.json').read_text())
    assert validation['accepted_biological_group_count'] is None
    report = dict(status='passed', cached_response_hashes_checked=checked, raw_records=2205,
                  species=9, raw_source_hash_unchanged=True,
                  scope='Integrity and table reconciliation only; not independent full reconstruction or biological validation.')
    (OUT/'verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
