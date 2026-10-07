"""Audit and derive entry/reference timelines from the completed bibliography.

No network requests. Ambiguous names remain unresolved rather than guessed.
"""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parents[1] / 'GenAge Project/models_genes/genage_models.csv'
DATA = ROOT / 'results/bibliography'
OUT = ROOT / 'results/timeline'


def derive():
    status = json.loads((DATA / 'collection_status.json').read_text())
    if not status['complete']:
        raise ValueError('Collection is incomplete')
    genes = pd.read_csv(SOURCE).fillna('').rename(columns={'GenAge ID': 'genage_id'})
    refs = pd.read_csv(DATA / 'references.csv')
    links = pd.read_csv(DATA / 'observation_backlinks.csv', dtype={'observation_source_id': str}).fillna('')
    genes['gene_label'] = genes['symbol'] + ', ' + genes['name'] + ' (' + genes['organism'] + ')'
    normalize = lambda s: ' '.join(s.split())
    genes['gene_label'] = genes.gene_label.map(normalize)
    links['gene_label'] = links.gene_label.map(normalize)
    unique = genes[~genes.gene_label.duplicated(keep=False)]
    links = links.merge(unique[['gene_label', 'genage_id', 'symbol', 'organism']], on='gene_label', how='left', validate='many_to_one')
    links['mapping_method'] = 'exact_full_label'
    # LibAge retains older gene names and symbol capitalization. Accept a
    # symbol/species match only when it identifies exactly one Build 21 entry.
    for index in links.index[links.genage_id.isna()]:
        label = links.at[index, 'gene_label']
        candidates = genes.iloc[0:0]
        for organism in genes.organism.unique():
            ending = ' (' + organism + ')'
            if label.endswith(ending):
                symbol = label[:-len(ending)].split(', ', 1)[0]
                candidates = genes[genes.organism.eq(organism) & genes.symbol.str.casefold().eq(symbol.casefold())]
                break
        if len(candidates) == 1:
            row = candidates.iloc[0]
            links.loc[index, ['genage_id', 'symbol', 'organism']] = [row.genage_id, row.symbol, row.organism]
            links.at[index, 'mapping_method'] = 'unique_symbol_species_casefold'
        else:
            links.at[index, 'mapping_method'] = 'unresolved'
    # Observation IDs from the pilot provide an independent exact-identity check.
    pilot_obs = pd.read_csv(ROOT / 'results/pilot/genage_observations_pilot.csv', dtype={'observation_source_id': str})
    pilot_refs = pd.read_csv(ROOT / 'results/pilot/genage_references_pilot.csv', dtype={'observation_source_id': str})
    pilot_map = pilot_obs.set_index('observation_source_id').genage_id.to_dict()
    mismatches = []
    for row in links.itertuples():
        if row.observation_source_id in pilot_map and pd.notna(row.genage_id):
            if row.genage_id != pilot_map[row.observation_source_id]:
                mismatches.append(dict(observation_source_id=row.observation_source_id,
                    libage_id=row.libage_id, label_genage_id=int(row.genage_id),
                    pilot_genage_id=int(pilot_map[row.observation_source_id]), gene_label=row.gene_label))
    expected_pairs = set(zip(pilot_refs.observation_source_id, pilot_refs.libage_id))
    actual_pairs = set(zip(links.observation_source_id, links.libage_id))
    missing_pilot = expected_pairs - actual_pairs
    OUT.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(mismatches, columns=['observation_source_id', 'libage_id', 'label_genage_id', 'pilot_genage_id', 'gene_label']).to_csv(OUT / 'pilot_identity_conflicts.csv', index=False)
    unresolved = links[links.genage_id.isna()]
    unresolved.to_csv(OUT / 'unresolved_links.csv', index=False)
    mapped = links[links.genage_id.notna()].copy()
    mapped['genage_id'] = mapped.genage_id.astype(int)
    mapped = mapped.merge(refs, on='libage_id', how='left', validate='many_to_one', suffixes=('_link', ''))
    mapped.to_csv(OUT / 'mapped_observation_references.csv', index=False)
    entry_refs = mapped.drop_duplicates(['genage_id', 'libage_id'])
    earliest = entry_refs.groupby('genage_id').publication_year.min().rename('earliest_associated_year')
    coverage = genes.merge(earliest, on='genage_id', how='left')
    coverage.to_csv(OUT / 'entry_coverage.csv', index=False)
    counts = coverage.groupby('organism').agg(total_entries=('genage_id', 'size'), dated_entries=('earliest_associated_year', 'count'))
    counts.to_csv(OUT / 'coverage_by_organism.csv')
    dated = coverage.dropna(subset=['earliest_associated_year']).copy()
    years = range(int(refs.publication_year.min()), int(refs.publication_year.max()) + 1)
    grid = pd.MultiIndex.from_product([sorted(genes.organism.unique()), years], names=['organism', 'year'])
    annual = dated.groupby(['organism', 'earliest_associated_year']).size()
    annual.index.names = ['organism', 'year']
    annual = annual.reindex(grid, fill_value=0).rename('entries_first_referenced').reset_index()
    annual['cumulative_entries'] = annual.groupby('organism').entries_first_referenced.cumsum()
    annual.to_csv(OUT / 'annual_entry_timeline.csv', index=False)
    papers = refs.drop_duplicates('libage_id').groupby('publication_year').size().reindex(years, fill_value=0)
    papers.rename('unique_libage_references').rename_axis('year').to_csv(OUT / 'annual_reference_timeline.csv')
    report = dict(total_entries=len(genes), dated_entries=len(dated), references=len(refs),
        unresolved_links=len(unresolved), missing_pilot_pairs=sorted(missing_pilot),
        pilot_identity_conflicts=mismatches,
        missing_reference_years=int(refs.publication_year.isna().sum()),
        build21_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        definition='Counts GenAge entries, not deduplicated biological genes; earliest associated reference is not necessarily discovery.',
        mapping_status='UNVALIDATED label-based candidates; do not use for scientific conclusions',
        final_ready=bool(len(unresolved)==0 and not missing_pilot and not mismatches and len(dated)==len(genes)))
    (OUT / 'audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    derive()
