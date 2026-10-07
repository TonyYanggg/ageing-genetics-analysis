"""Derive timelines using GenAge observation IDs, never LibAge display names."""
from pathlib import Path
import hashlib
import json
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parents[1] / 'GenAge Project/models_genes/genage_models.csv'
BIB = ROOT / 'results/bibliography'
IDENTITY = ROOT / 'results/identity_verification'
SNAPSHOT = ROOT / 'results/hagr_r_snapshot'
OUT = ROOT / 'results/verified_timeline'


def derive():
    for item in json.loads((SNAPSHOT / 'manifest.json').read_text()):
        assert hashlib.sha256((SNAPSHOT / item['file']).read_bytes()).hexdigest() == item['sha256']
    genes = pd.read_csv(SOURCE).rename(columns={'GenAge ID': 'genage_id'})
    models = pd.read_csv(SNAPSHOT / 'models.csv')
    identity_check = genes.merge(models, left_on='genage_id', right_on='id_gene', how='outer', suffixes=('_build', '_snapshot'), indicator=True)
    assert identity_check['_merge'].eq('both').all()
    assert identity_check.symbol_build.eq(identity_check.symbol_snapshot).all()
    assert identity_check.organism_build.eq(identity_check.organism_snapshot).all()
    longevity = pd.read_csv(SNAPSHOT / 'longevity.csv')
    obs = longevity[longevity.gene_id.notna()].rename(columns={'id': 'observation_source_id', 'gene_id': 'genage_id'}).copy()
    obs['genage_id'] = obs.genage_id.astype(int)
    obs['observation_source_id'] = obs.observation_source_id.astype(int).astype(str)
    obs = obs.merge(genes[['genage_id', 'symbol', 'organism']], on='genage_id', validate='many_to_one')
    obs['source_sha256'] = hashlib.sha256((SNAPSHOT / 'longevity.csv').read_bytes()).hexdigest()
    live = pd.read_csv(IDENTITY / 'verified_observations.csv', dtype={'observation_source_id': str})
    pilot_obs = pd.read_csv(ROOT / 'results/pilot/genage_observations_pilot.csv', dtype={'observation_source_id': str})
    live = pd.concat([live, pilot_obs]).drop_duplicates('observation_source_id')
    checked = live.merge(obs[['observation_source_id', 'genage_id']], on='observation_source_id', how='left', suffixes=('_live', '_snapshot'), validate='one_to_one')
    assert checked.genage_id_live.eq(checked.genage_id_snapshot).all(), 'Official snapshot disagrees with live observation identity'
    refs = pd.read_csv(BIB / 'references.csv')
    links = pd.read_csv(BIB / 'observation_backlinks.csv', dtype={'observation_source_id': str})
    OUT.mkdir(parents=True, exist_ok=True)
    identities = obs[['observation_source_id', 'genage_id', 'symbol', 'organism', 'source_sha256']].drop_duplicates()
    if identities.observation_source_id.duplicated().any():
        raise ValueError('An observation ID maps to multiple identity records; investigate before joining')
    mapped = links.merge(identities, on='observation_source_id', how='left', validate='many_to_one')
    unresolved = mapped[mapped.genage_id.isna()]
    unresolved.to_csv(OUT / 'unresolved_backlinks.csv', index=False)
    gaps = pd.read_csv(ROOT / 'results/snapshot_gap_checks/gap_checks.csv', dtype={'observation_source_id': str})
    for row in gaps[gaps.genage_id.notna()].itertuples():
        assert row.observation_present
        expected_id = identities.set_index('observation_source_id').loc[row.observation_source_id, 'genage_id']
        assert row.genage_id == expected_id
    exceptions = unresolved.merge(gaps, on='observation_source_id', how='left', suffixes=('_bibliography', '_destination'), validate='one_to_one')
    assert exceptions.status_code.eq(200).all() and exceptions.error.eq('No usable gene entry at backlink destination').all()
    # Explicit entry-not-found content distinguishes dangling links from network failures.
    for row in exceptions.itertuples():
        assert 'Entry not found' in (ROOT / row.source_file).read_text(encoding='utf-8')
    exceptions.to_csv(OUT / 'excluded_dangling_backlinks.csv', index=False)
    no_reference = obs[~obs.observation_source_id.isin(links.observation_source_id)]
    no_reference.to_csv(OUT / 'observations_without_bibliography_reference.csv', index=False)
    assert set(no_reference.observation_source_id) == {'2158'}, 'Unexpected observation without reference requires review'
    assert json.loads((ROOT / 'results/snapshot_gap_checks/raw/observation_2158.json').read_text()) == []
    if len(unresolved) or len(no_reference):
        print(f'Coverage exceptions: {len(unresolved)} unresolved backlinks; {len(no_reference)} observations without a bibliography reference')
    mapped = mapped[mapped.genage_id.notna()].copy()
    mapped['genage_id'] = mapped.genage_id.astype(int)
    mapped = mapped.merge(refs, on='libage_id', how='left', validate='many_to_one', suffixes=('_link', ''))
    pilot = pd.read_csv(ROOT / 'results/pilot/genage_references_pilot.csv', dtype={'observation_source_id': str})
    expected = set(zip(pilot.genage_id, pilot.observation_source_id, pilot.libage_id))
    actual = set(zip(mapped.genage_id, mapped.observation_source_id, mapped.libage_id))
    missing_pilot = expected - actual
    if missing_pilot:
        raise ValueError(f'Pilot triple mismatch: {missing_pilot}')
    if mapped.publication_year.isna().any():
        raise ValueError('A mapped reference lacks a publication year')
    mapped.to_csv(OUT / 'verified_observation_references.csv', index=False)
    previous = pd.read_csv(ROOT / 'results/timeline/mapped_observation_references.csv', dtype={'observation_source_id': str})
    comparison = previous[['observation_source_id', 'libage_id', 'genage_id', 'gene_label']].merge(
        mapped[['observation_source_id', 'libage_id', 'genage_id']], on=['observation_source_id', 'libage_id'],
        suffixes=('_label', '_verified'), validate='one_to_one')
    corrections = comparison[comparison.genage_id_label.ne(comparison.genage_id_verified)]
    corrections.to_csv(OUT / 'corrected_label_assignments.csv', index=False)
    assert set(corrections.observation_source_id).issubset(set(gaps.loc[gaps.observation_present.eq(True), 'observation_source_id']))
    entry_refs = mapped.drop_duplicates(['genage_id', 'libage_id'])
    entry_refs.to_csv(OUT / 'entry_references.csv', index=False)
    earliest = entry_refs.groupby('genage_id').publication_year.min().rename('earliest_associated_year')
    coverage = genes.merge(earliest, on='genage_id', how='left')
    coverage.to_csv(OUT / 'entry_coverage.csv', index=False)
    missing_entries = set(coverage.loc[coverage.earliest_associated_year.isna(), 'genage_id'])
    assert missing_entries == {2081, 2173}, 'Unexpected undated entry requires review'
    for gene_id in missing_entries:
        page = BeautifulSoup((ROOT / f'results/snapshot_gap_checks/raw/gene_{gene_id}.html').read_bytes(), 'html.parser')
        assert f'Gene information (ID: {gene_id})' in page.get_text(' ', strip=True)
        assert not page.select('ul.references-list')
    coverage.groupby('organism').agg(total_entries=('genage_id', 'size'), dated_entries=('earliest_associated_year', 'count')).to_csv(OUT / 'coverage_by_organism.csv')
    dated = coverage.dropna(subset=['earliest_associated_year'])
    years = range(int(refs.publication_year.min()), int(refs.publication_year.max()) + 1)
    grid = pd.MultiIndex.from_product([sorted(genes.organism.unique()), years], names=['organism', 'year'])
    annual = dated.groupby(['organism', 'earliest_associated_year']).size()
    annual.index.names = ['organism', 'year']
    annual = annual.reindex(grid, fill_value=0).rename('entries_first_referenced').reset_index()
    annual['cumulative_entries'] = annual.groupby('organism').entries_first_referenced.cumsum()
    assert int(annual.groupby('organism').cumulative_entries.last().sum()) == len(dated)
    annual.to_csv(OUT / 'annual_entry_timeline.csv', index=False)
    used_refs = refs[refs.libage_id.isin(mapped.libage_id)]
    used_refs.groupby('publication_year').size().reindex(years, fill_value=0).rename('unique_libage_references').rename_axis('year').to_csv(OUT / 'annual_reference_timeline.csv')
    # Preserve every tied earliest publication for each entry's provenance.
    earliest_rows = entry_refs.merge(earliest, on='genage_id')
    earliest_rows[earliest_rows.publication_year.eq(earliest_rows.earliest_associated_year)].to_csv(OUT / 'earliest_reference_provenance.csv', index=False)
    report = dict(total_entries=len(genes), dated_entries=len(dated), references=len(used_refs),
        verified_observations=len(obs), mapped_links=len(mapped), unresolved_links=len(unresolved),
        excluded_dangling_links=len(exceptions), live_sample_entries=int(live.genage_id.nunique()),
        corrected_label_assignments=len(corrections),
        live_sample_observations=len(live), snapshot_commit='4e83b50142ea4b8c6625a7a86b6005ff3f606ed6',
        observations_without_reference=len(no_reference), missing_pilot_pairs=[], pilot_identity_conflicts=[],
        missing_reference_years=0, build21_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [SNAPSHOT / 'models.csv', SNAPSHOT / 'longevity.csv', BIB / 'references.csv', BIB / 'observation_backlinks.csv']},
        mapping_status='Official HAGR-R observation-to-entry IDs; all entry identities match Build 21; live sample and targeted gap checks agree',
        preview_ready=bool(not missing_pilot and len(exceptions)==len(unresolved)),
        final_ready=bool(not len(unresolved) and not len(no_reference) and len(dated)==len(genes)),
        definition='GenAge entries by earliest associated publication, not verified historical gene discovery; distinct LibAge references for publications')
    (OUT / 'audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    derive()
