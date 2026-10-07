# Model-only orthology source audit

Superseded as the primary method by `../oma_pairwise/README.md`. Retained as evidence for source selection; do not use the 181 broad source groups as the current Task 4 result.

Status: implementation and source collection complete for the provisional OrthoDB screen; biological interpretation **not approved**. See `../../notes/source_validation.md`. All nine GenAge organisms are included; no human dataset is used.

## Inputs and provenance

Immutable GenAge Build 21 CSV: 2,205 records. SHA-256 is recorded in `collection_status.json`. OrthoDB 12.2 selected assemblies are documented in `source_species.csv`. `source_cache/` contains raw JSON and paired metadata files with exact request parameters, retrieval time, headers and checksums. OMA spot checks are retrieval-frozen live responses. Historical human-anchor files remain outside this folder.

## Reproduce from project root

Run these Python scripts in order (requests, matplotlib and numpy required):

```powershell
python 'task/task 4/scripts/collect_model_orthology.py'
python 'task/task 4/scripts/analyze_model_orthology.py'
python 'task/task 4/scripts/collect_model_group_members.py'
python 'task/task 4/scripts/audit_orthology_source.py'
python 'task/task 4/scripts/plot_model_orthology.py'
python 'task/task 4/scripts/verify_model_audit.py'
```

Collectors reuse checksum-verified caches; network access is needed only for missing responses. Do not silently replace caches with a different release.

## Outputs

- `inventory.csv`, `genage_records.csv`: all-nine inventory and derived source copy.
- `mapping_memberships.csv`, `collection_status.json`, `source_species.csv`: collected identifiers, memberships and assembly provenance.
- `record_mapping_audit.csv`, `coverage.csv`, `ambiguous_group_mappings.csv`, `taxonomy_mismatches.csv`: mapping decisions, denominators and exclusions.
- `groups_all_levels.csv`, `groups_eukaryota.csv`, `recurrent_groups_eukaryota.csv`: provisional source-defined group summaries, NOT validated biological candidates. Nested levels must not be summed.
- `group_record_evidence.csv`, `group_species_cells.csv`, `recurrent_group_members.csv`, `recurrent_atlas_cells.csv`: record-level evidence and source membership, separate from lifespan evidence.
- `figure_displayed_groups.csv`, `figure_recurrence_distribution.csv`: plotted source data. Draft PNG/SVG/PDF/TIFF exports are in `../../figures/model_orthology/`.
- `analysis_summary.json`, `source_validation.json`, `verification.json`: numerical summary, cross-resource conflict, and narrow integrity checks.

## Limitations and next step

The 181 recurrent Eukaryota-level source groups are an audit finding, not a defensible count of conserved longevity genes. Coq2/HSP104 grouping disagrees with OMA. Resolve source/group resolution across all nine organisms before accepting results. Missing IDs and absent mappings are unknown, not biological absence. P/A sensitivity analysis and experimental-context review remain unfinished.

## History

2026-09-18: first model-only pipeline, all-nine coverage audit and source conflict documented; shared raw files and historical analysis preserved.
