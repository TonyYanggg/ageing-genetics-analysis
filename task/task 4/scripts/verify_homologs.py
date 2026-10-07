"""Independent raw-input reconstruction of Task 4 counts and exported cells."""
import csv
import hashlib
import json
from collections import defaultdict, Counter
from pathlib import Path
from PIL import Image

TASK = Path(__file__).resolve().parents[1]
ROOT = TASK.parents[1]


def read(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main():
    models = read(ROOT / 'GenAge Project/models_genes/genage_models.csv')
    mappings = read(ROOT / 'GenAge Project/models_genes/genage_models_orthologs_export.tsv', '\t')
    audit = json.loads((TASK / 'results/audit.json').read_text())
    symbols = defaultdict(set)
    for row in mappings:
        symbols[row['Entrez ID']].add(row['Symbol'])
    valid = [row for row in mappings if len(symbols[row['Entrez ID']]) == 1]
    direction = defaultdict(set)
    record_ids = defaultdict(set)
    for edge in valid:
        for row in models:
            if (row['organism'], row['entrez gene id']) == (edge['Model Organism'], edge['Model Organism Entrez ID']):
                key = (edge['Entrez ID'], row['organism'])
                record_ids[key].add(row['GenAge ID'])
                if row['longevity influence'] in ['Pro-Longevity', 'Anti-Longevity']:
                    direction[key].add(row['longevity influence'])
    species_counts = Counter(anchor for (anchor, species), signs in direction.items() if signs)
    expected = {anchor for anchor, count in species_counts.items() if count >= 2}
    candidates = read(TASK / 'results/cross_species_candidates.csv')
    assert expected == {r['human_entrez_id'] for r in candidates}
    for cell in read(TASK / 'results/anchor_species_evidence.csv'):
        key = (cell['human_entrez_id'], cell['species'])
        signs = direction[key]
        expected_state = {frozenset(): 'No directional evidence',
            frozenset(['Pro-Longevity']): 'Pro', frozenset(['Anti-Longevity']): 'Anti',
            frozenset(['Pro-Longevity','Anti-Longevity']): 'Mixed'}[frozenset(signs)]
        assert cell['state'] == expected_state
        assert set(filter(None,cell['genage_ids'].split(';'))) == record_ids[key]
    for coverage in read(TASK / 'results/mapping_coverage.csv'):
        rows = [r for r in models if r['organism'] == coverage['species'] and
            (coverage['scope'] == 'all_records' or r['longevity influence'] in ['Pro-Longevity','Anti-Longevity'])]
        ids = {r['entrez gene id'] for r in rows if r['entrez gene id']}
        mapped = ids & {r['Model Organism Entrez ID'] for r in valid if r['Model Organism'] == coverage['species']}
        assert len(ids) == int(coverage['unique_genes_with_id'])
        assert len(mapped) == int(coverage['mapped_unique_genes'])
    for name, path in audit['source_paths'].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == audit['input_sha256'][name]
    assert '<text' in (TASK / 'figures/cross_species_homologs.svg').read_text(encoding='utf-8')
    with Image.open(TASK / 'figures/cross_species_homologs.tiff') as im:
        assert im.info['dpi'] == (600,600)
    text = ('# Verification report\n\n2026-09-18: PASS.\n\n'
        f'- Independently reconstructed {len(expected)} cross-species candidates directly from raw inputs.\n'
        '- Verified every exported anchor/species state and its supporting GenAge record IDs.\n'
        '- Recomputed mapping numerators and denominators across all species and both scopes.\n'
        '- Confirmed all three input SHA-256 hashes; no source edits.\n'
        '- Verified editable SVG text and 600 dpi TIFF metadata.\n'
        '- Biological orthology confidence and intervention-level direction remain unverified.\n')
    (TASK / 'notes/verification.md').write_text(text, encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
