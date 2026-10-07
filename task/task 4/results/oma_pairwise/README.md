# Explicit model-organism ortholog pairs

**Superseded direct-ID baseline.** The current canonical analysis is in [`refined/`](refined/README.md), with 228 pairs after exact stable-ID recovery. The 219-pair files and narrative below are preserved as a baseline audit. `source_cache/` is shared with the refined pipeline and now also contains identity, exact cross-reference and experiment-review responses. Do not confuse baseline counts with current outputs.

## Status and interpretation

Historical **direct-ID baseline** for the subsequent refined OMA audit, not the main focused review figure. The refined 228-pair table is under `refined/`. This baseline has 219 unique gene pairs with qualifying GenAge lifespan evidence on both sides; 136 meet the narrower two-sided Pro/Anti criterion. These pairs share genes: neither count is a count of independent conserved genes or orthogroups. The complete all-nine biological survey remains unfinished because of source/identifier gaps.

## Inputs

- Immutable GenAge Build 21 model CSV: 2,205 records across nine organisms. Source SHA-256 in `collection_status.json`.
- OMA All.May2026 / API 1.11: exact NCBI Gene cross-references, seven selected reference genomes, all 21 genome-pair endpoints. Complete pagination gives 124,274 source protein pairs. Before/after release checks match.
- Ensembl release 116: separate hamster-to-mouse/zebrafish calls and Podospora identifier checks. No Ensembl HOGs are merged with OMA groups. The two hamster target genes have no matching GenAge lifespan records and add no recurrent pair.

No human gene data or human anchors are used. Source genomes and taxon/strain choices appear in `source_species.csv`. All requests, parameters, response payloads, retrieval times and checksums are frozen in `source_cache/`.

## Reproduction

From the project root, with Python, requests, numpy and matplotlib:

```powershell
python 'task/task 4/scripts/collect_oma_pairwise.py'
python 'task/task 4/scripts/analyze_oma_pairwise.py'
python 'task/task 4/scripts/check_orthology_gaps.py'
python 'task/task 4/scripts/plot_oma_pairwise.py'
python 'task/task 4/scripts/verify_oma_pairwise.py'
```

Cached responses are hash-checked and reused. Internet access is required only for missing responses. The analyzer also reads the preserved earlier model-only mapping audit for the source comparison; historical human-anchor outputs are not inputs.

## Outputs

| Output | Purpose/status |
|---|---|
| `longevity_pairs.csv` | Primary OMA-supported pairs with lifespan evidence on both sides |
| `directional_sensitivity_pairs.csv` | Two-sided Pro/Anti subset, not the discovery universe |
| `all_genage_pairs.csv` | All linked GenAge pairs; currently all 219 meet lifespan eligibility |
| `record_mapping_audit.csv`, `coverage.csv` | All 2,205 records and all nine organisms, exclusions and source gaps |
| `identifier_mappings.csv`, `source_pairs.csv`, `source_species.csv` | Collected source evidence; source pairs retain API orientation |
| `orthodb_group_crosscheck.csv` | Every old recurrent source group compared with explicit OMA calls; missing calls are not disproven orthology |
| `supplementary_ensembl_pairs.csv`, `source_gap_checks.json` | Separate supplementary calls and unresolved-ID checks |
| `analysis_summary.json`, `collection_status.json`, `verification.json` | Counts, collection completeness and independent reconstruction |
| `figure_pair_counts.csv` | All 81 plotted matrix cells; NA distinguishes source gaps |
| `../../figures/oma_pairwise/pairwise_evidence.*` | PNG/SVG/PDF/TIFF coverage and pair-count summary |

Derived relation types are oriented from species1 to species2. Multiple isoforms and duplicate records cannot inflate the gene-pair count. Genes involved in ambiguous many-gene-to-one-protein identifier links are excluded, not arbitrarily selected. Full annotations and GenAge record IDs remain traceable.

## Verification and limitations

The independent verifier reads the immutable input and raw cached JSON rather than reusing the analysis functions. It reconstructs all 219 gene pairs and the 136-pair sensitivity set, verifies complete pagination, source hashes, annotation provenance, relation orientation and figure counts. Numerical reproducibility does not establish experimental conservation.

Podospora and hamster are absent from OMA; all nine remain in the inventory. Ensembl confirms the hamster gene and two model orthologs, but these are not new recurrent-longevity pairs. Podospora exact NCBI-ID queries return no hits. Additional cross-references and comparative resources require further resolution, not inferred absence. Other unmatched and ambiguous genes also remain in the audit. Experimental manipulation/context validation remains pending.

## History

2026-09-18: introduced explicit pairs to address the ancient-group problem; recovered C. briggsae daf-2; retained both earlier analyses unchanged. Method decision and figure caption: `../../notes/pairwise_resolution.md`.
