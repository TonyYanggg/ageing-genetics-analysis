# Refined explicit protein-ortholog pairs

## Status

Verified broad source-and-coverage audit, now supplying the [complete four-organism Task 4 overview](../../../figures/complete_four_species/complete_ortholog_overview.png) and its [206-row direct-pair table](../../complete_four_species/all_qualifying_pairs.csv). The earlier broad figures remain available as supporting history. This package contains 228 unique model-gene pairs with GenAge lifespan evidence on both sides across seven OMA-represented organisms. There are 145 two-sided Pro/Anti pairs: 75 Anti/Anti, 35 Pro/Pro and 35 discordant/mixed; the other 83 have unresolved directional classification. Counts overlap in genes and are not independent orthogroups. Complete package does not mean complete biological coverage or mechanistic validation.

## Inputs and method

Immutable Build 21 model CSV (2,205 records, nine organisms); OMA All.May2026/API 1.11 (seven represented genomes, 124,274 protein pairs across all 21 genome pairs). Parent `source_cache/` freezes responses, parameters, retrieval dates, release checks and SHA-256 hashes. Exact NCBI Gene mappings are supplemented by 63 taxon-matched stable-ID bridges from NCBI Datasets to literal OMA cross-references. The bridge audit is in the parent folder. Model-organism authority IDs and yeast systematic locus tags are accepted; symbol similarity is not.

Ambiguous gene-to-protein links are excluded before counting unique unordered species+Entrez pairs. Record and isoform multiplicity cannot inflate counts. Qualifying lifespan annotations are Increase, Decrease, or Increase and Decrease on both sides; P/A labels define only a secondary subset. Relation cardinalities are oriented from species1 to species2. No graph components are promoted to orthogroups. No human catalogue/anchor input is used.

## Reproduction

From project root, Python with requests, numpy, matplotlib, Pillow and PyMuPDF (`fitz`, PDF QA only). No generating commit is available: this workspace is not a Git repository. Existing installed packages were used; no dependencies were installed for this refinement.

```powershell
# Existing frozen collection is reused; these source steps need network only for missing caches.
python 'task/task 4/scripts/audit_gene_identities.py'
python 'task/task 4/scripts/recover_oma_identifiers.py'
python 'task/task 4/scripts/collect_experiment_review.py'
python 'task/task 4/scripts/collect_review_fulltext.py'

# Analysis and verification below run from local frozen inputs.
python 'task/task 4/scripts/analyze_oma_pairwise.py' --with-bridges
python 'task/task 4/scripts/build_resolution_audit.py'
python 'task/task 4/scripts/plot_oma_pairwise.py' --with-bridges
python 'task/task 4/scripts/plot_model_examples.py'
python 'task/task 4/scripts/verify_oma_pairwise.py' --with-bridges
python 'task/task 4/scripts/verify_review_package.py'
```

Initial OMA and supplementary collection commands are retained in the [parent baseline README](../README.md). The analyzer also reads the preserved OrthoDB audit for comparison, not for pair inclusion. The manually curated `../../../notes/experiment_review.tsv` is a reviewed input, not machine-extracted experimental truth.

## Output provenance

| Output | Status | Generator and purpose |
|---|---|---|
| `longevity_pairs.csv`, `all_genage_pairs.csv` | Current | `analyze_oma_pairwise.py --with-bridges`; explicit evidence pairs |
| `directional_sensitivity_pairs.csv` | Supporting | Same analyzer; two-sided P/A subset |
| `coverage.csv`, `record_mapping_audit.csv`, `analysis_summary.json` | Supporting | Same analyzer; all-nine inventory and mapping |
| `resolution_ledger.csv`, `remaining_gaps.csv`, `resolution_summary.json` | Supporting | `build_resolution_audit.py`; all 2,205 records and 67 remaining records classified |
| `orthodb_group_crosscheck.csv` | Supporting | Analyzer; comparison with superseded broad source groups |
| `figure_pair_counts.csv` | Supporting | `plot_oma_pairwise.py --with-bridges`; all 81 matrix cells |
| `figure_example_pairs.csv` | Supporting | `plot_model_examples.py`; exact plotted pairs and selection reasons |
| `verification.json`, `review_package_verification.json` | Supporting | Independent reconstruction and package checks |
| `../../../figures/oma_refined/` | Supporting broad audit | Python plot scripts; editable SVG/PDF, PNG and 600-dpi TIFF |

## Interpretation and limitations

The source-specific protein analysis is not a comprehensive all-nine orthology atlas. Podospora/hamster source gaps, noncoding loci, current pseudogene annotations, retired/missing IDs, ambiguous links and catalogue/experiment conflicts remain explicit. A mapped gene without a recovered evidence pair does not lack orthologs or lifespan function. Gene-type checks apply to baseline unmapped records, not the complete catalogue. Primary-paper review is purposive and limited in depth, not exhaustive. See [review note, captions and acceptance work](../../../notes/refined_review.md) and [history](../../../CHANGELOG.md).
