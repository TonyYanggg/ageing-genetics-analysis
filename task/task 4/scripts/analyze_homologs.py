"""Reproduce the frozen-snapshot human-anchored GenAge analysis and figure."""
from pathlib import Path
import csv
import hashlib
import json
from collections import Counter, defaultdict
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

TASK = Path(__file__).resolve().parents[1]
ROOT = TASK.parents[1]
SOURCES = {
    'models': ROOT / 'GenAge Project/models_genes/genage_models.csv',
    'orthologs': ROOT / 'GenAge Project/models_genes/genage_models_orthologs_export.tsv',
    'human': ROOT / 'GenAge Project/human_genes/genage_human.csv',
}
SPECIES = ['Saccharomyces cerevisiae', 'Caenorhabditis elegans',
           'Drosophila melanogaster', 'Mus musculus']
SHORT = ['Yeast', 'Worm', 'Fly', 'Mouse']
PRIMARY = {'Pro-Longevity', 'Anti-Longevity'}
RESULTS = TASK / 'results'
FIGURES = TASK / 'figures'


def read(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def write(name, rows, fields):
    with (RESULTS / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def joined(values):
    return ';'.join(sorted(set(values)))


def hashes():
    return {key: hashlib.sha256(path.read_bytes()).hexdigest()
            for key, path in SOURCES.items()}


def main():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    before = hashes()
    models = read(SOURCES['models'])
    raw_edges = read(SOURCES['orthologs'], '\t')
    humans = {r['entrez gene id']: r for r in read(SOURCES['human'])}
    index = defaultdict(list)
    for row in models:
        if row['entrez gene id']:
            index[(row['organism'], row['entrez gene id'])].append(row)
    assert len({r['GenAge ID'] for r in models}) == len(models)
    symbols_by_id = defaultdict(set)
    for row in raw_edges:
        symbols_by_id[row['Entrez ID']].add(row['Symbol'])
    ambiguous_ids = {key for key, symbols in symbols_by_id.items() if len(symbols) > 1}
    quarantined = [r for r in raw_edges if r['Entrez ID'] in ambiguous_ids]
    write('ambiguous_human_mappings.csv', quarantined, [k for k in raw_edges[0] if k is not None])
    edges = {}
    for row in raw_edges:
        assert row['Species'] == 'Homo sapiens'
        if row['Entrez ID'] in ambiguous_ids:
            continue
        key = (row['Entrez ID'], row['Model Organism'], row['Model Organism Entrez ID'])
        if key in edges:
            assert row == edges[key], 'Conflicting duplicate mapping'
        edges[key] = row
    anchors_for_gene = defaultdict(set)
    for human_id, species, model_id in edges:
        anchors_for_gene[(species, model_id)].add(human_id)
    merged, unmatched, mismatches = [], [], []
    groups = defaultdict(list)
    for (human_id, species, model_id), edge in edges.items():
        matches = index.get((species, model_id), [])
        if not matches:
            unmatched.append(edge)
        for row in matches:
            item = dict(human_entrez_id=human_id, human_symbol=edge['Symbol'],
                        species=species, model_entrez_id=model_id,
                        model_symbol=row['symbol'], export_model_symbol=edge['Model Organism Symbol'],
                        genage_id=row['GenAge ID'], longevity_influence=row['longevity influence'],
                        lifespan_effect=row['lifespan effect'],
                        primary_evidence=row['longevity influence'] in PRIMARY,
                        human_anchors_for_model_gene=len(anchors_for_gene[(species, model_id)]))
            merged.append(item)
            groups[human_id].append(item)
            if item['model_symbol'] != item['export_model_symbol']:
                mismatches.append(item)
    write('joined_evidence.csv', merged, list(merged[0]))
    write('symbol_disagreements.csv', mismatches, list(merged[0]))
    write('unmatched_mappings.csv', unmatched, [k for k in raw_edges[0] if k is not None])
    coverage = []
    for species in sorted({r['organism'] for r in models}):
        rows = [r for r in models if r['organism'] == species]
        for scope in ['all_records', 'directional_longevity']:
            selected = [r for r in rows if scope == 'all_records' or r['longevity influence'] in PRIMARY]
            ids = {r['entrez gene id'] for r in selected if r['entrez gene id']}
            mapped = {i for i in ids if (species, i) in anchors_for_gene}
            coverage.append(dict(species=species, scope=scope, records=len(selected),
                records_missing_entrez=sum(not r['entrez gene id'] for r in selected),
                unique_genes_with_id=len(ids), mapped_unique_genes=len(mapped),
                unmapped_unique_genes=len(ids - mapped),
                mapped_fraction=len(mapped)/len(ids) if ids else '',
                species_in_export=species in {e['Model Organism'] for e in raw_edges}))
    write('mapping_coverage.csv', coverage, list(coverage[0]))
    unmapped = [r for r in models if (r['organism'], r['entrez gene id']) not in anchors_for_gene]
    write('unmapped_genage_records.csv', unmapped, list(models[0]))
    summaries, cells = [], []
    for human_id, evidence in groups.items():
        summary = dict(human_entrez_id=human_id, human_symbol=joined(r['human_symbol'] for r in evidence),
            in_human_genage=human_id in humans,
            human_genage_reason=humans.get(human_id, {}).get('why', ''))
        directions = set()
        n_species = 0
        for species, short in zip(SPECIES, SHORT):
            all_rows = [r for r in evidence if r['species'] == species]
            primary = [r for r in all_rows if r['primary_evidence']]
            signs = {r['longevity_influence'] for r in primary}
            state = ('Mixed' if len(signs) == 2 else
                     'Pro' if signs == {'Pro-Longevity'} else
                     'Anti' if signs == {'Anti-Longevity'} else 'No directional evidence')
            n_species += bool(signs)
            directions.update(signs)
            summary[short.lower() + '_state'] = state
            summary[short.lower() + '_genes'] = joined(r['model_symbol'] for r in primary)
            cells.append(dict(human_entrez_id=human_id, human_symbol=summary['human_symbol'],
                species=species, state=state, directional_genes=joined(r['model_symbol'] for r in primary),
                directional_entrez_ids=joined(r['model_entrez_id'] for r in primary),
                genage_ids=joined(r['genage_id'] for r in all_rows),
                all_longevity_annotations=joined(r['longevity_influence'] for r in all_rows)))
        summary.update(directional_species=n_species,
            cross_species_candidate=n_species >= 2,
            direction_agreement=('Pro' if directions == {'Pro-Longevity'} else 'Anti')
                if n_species >= 2 and len(directions) == 1 else
                'Mixed directions' if n_species >= 2 else 'Not cross-species',
            unique_model_genes=len({(r['species'], r['model_entrez_id']) for r in evidence}),
            shares_model_gene_with_other_anchor=any(r['human_anchors_for_model_gene'] > 1 for r in evidence))
        summaries.append(summary)
    summaries.sort(key=lambda r: (-r['directional_species'], r['human_symbol'], r['human_entrez_id']))
    candidates = [r for r in summaries if r['cross_species_candidate']]
    write('human_anchor_summary.csv', summaries, list(summaries[0]))
    write('cross_species_candidates.csv', candidates, list(summaries[0]))
    write('anchor_species_evidence.csv', cells, list(cells[0]))
    distribution = Counter(r['directional_species'] for r in summaries)
    write('species_count_distribution.csv', [{'directional_species': n, 'human_anchors': distribution[n]}
          for n in range(5)], ['directional_species', 'human_anchors'])
    displayed = [r for r in candidates if r['directional_species'] >= 3][:25]
    if not displayed:
        displayed = candidates[:25]
    write('figure_displayed_candidates.csv', displayed, list(summaries[0]))
    plot(coverage, distribution, displayed)
    assert before == hashes(), 'Source files changed'
    assert sum(distribution.values()) == len(groups)
    assert not unmatched, 'Resolve unmatched mappings before interpreting results'
    audit = dict(input_sha256=before, source_paths={k: str(v.relative_to(ROOT)) for k,v in SOURCES.items()},
        python=platform.python_version(), matplotlib=matplotlib.__version__,
        model_records=len(models), raw_mapping_rows=len(raw_edges), unique_mapping_edges=len(edges),
        quarantined_mapping_rows=len(quarantined), ambiguous_human_ids=sorted(ambiguous_ids),
        duplicate_mapping_rows=len(raw_edges)-len(quarantined)-len(edges), joined_evidence_rows=len(merged),
        unmatched_mapping_rows=len(unmatched), symbol_disagreements=len(mismatches),
        duplicate_model_keys={str(k): len(v) for k,v in index.items() if len(v)>1},
        human_anchors=len(groups), candidates=len(candidates),
        candidate_direction_agreement=dict(Counter(r['direction_agreement'] for r in candidates)),
        species_count_distribution=dict(distribution), displayed_candidates=len(displayed),
        inputs_unchanged=True)
    (RESULTS / 'audit.json').write_text(json.dumps(audit, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(audit, indent=2))


def plot(coverage, distribution, displayed):
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 7, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .6})
    fig = plt.figure(figsize=(183/25.4, min(230, 103+len(displayed)*4.1)/25.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, max(1.8, len(displayed)*.17)],
                           left=.15, right=.97, bottom=.10, top=.87, hspace=.48, wspace=.55)
    a, b, c = fig.add_subplot(grid[0,0]), fig.add_subplot(grid[0,1]), fig.add_subplot(grid[1,:])
    selected = [next(r for r in coverage if r['species']==s and r['scope']=='directional_longevity') for s in SPECIES]
    values = [100*r['mapped_fraction'] for r in selected]
    a.barh(SHORT, values, color='#608C9C', height=.6)
    a.invert_yaxis()
    a.set_xlim(0, 113)
    a.set_xticks([0, 25, 50, 75, 100])
    a.set_xlabel('Directional genes mapped (%)')
    for i,r in enumerate(selected):
        a.text(values[i]+2, i, f"{r['mapped_unique_genes']}/{r['unique_genes_with_id']}", va='center', fontsize=6)
    a.set_title('a   Mapping coverage', loc='left', fontweight='bold', pad=10)
    b.bar(range(5), [distribution[n] for n in range(5)], color='#608C9C', width=.65)
    for n in range(5):
        b.text(n, distribution[n], str(distribution[n]), ha='center', va='bottom', fontsize=7)
    b.set_ylim(0, max(distribution.values())*1.2)
    b.set_xticks(range(5))
    b.set_xlabel('Species with directional evidence')
    b.set_ylabel('Human gene anchors')
    b.set_title('b   Cross-species recurrence', loc='left', fontweight='bold', pad=10)
    states = {'No directional evidence': 0, 'Pro': 1, 'Anti': 2, 'Mixed': 3}
    matrix = [[states[r[s.lower()+'_state']] for s in SHORT] for r in displayed]
    c.imshow(matrix, cmap=ListedColormap(['#F1F1F1','#A4CAD5','#E7BB8C','#B9ADC8']),
             vmin=0, vmax=3, aspect='auto', interpolation='nearest')
    c.set_xticks(range(4), SHORT)
    c.xaxis.tick_top()
    c.set_yticks(range(len(displayed)), [r['human_symbol'] for r in displayed])
    c.tick_params(length=0, pad=5)
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            c.text(j,i,['–','P','A','M'][val], ha='center',va='center',fontsize=7)
    c.set_xticks(np.arange(-.5,4,1),minor=True)
    c.set_yticks(np.arange(-.5,len(displayed),1),minor=True)
    c.grid(which='minor', color='white', linewidth=1.2)
    c.tick_params(which='minor',length=0)
    for spine in c.spines.values():
        spine.set_visible(False)
    c.set_title('c   Human anchors with broad model-species evidence', loc='left', fontweight='bold', pad=26)
    fig.text(.15,.042,'P  Pro-longevity     A  Anti-longevity     M  Mixed     –  No directional evidence', fontsize=7)
    fig.text(.15,.018,'Human anchors overlap; shared mapping does not establish conserved function.', fontsize=6.5)
    fig.suptitle('Longevity evidence across human-anchored homolog mappings', x=.15, ha='left', fontsize=10)
    for ext in ['svg','pdf','png','tiff']:
        options = {'dpi': 600 if ext=='tiff' else 200}
        if ext=='tiff':
            options['pil_kwargs'] = {'compression':'tiff_lzw'}
        fig.savefig(FIGURES / ('cross_species_homologs.'+ext), **options)
    plt.close(fig)


if __name__ == '__main__':
    main()
