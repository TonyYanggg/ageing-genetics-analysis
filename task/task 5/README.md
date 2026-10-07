# Task 5 — Human ageing genes and disease associations

## Request

Map disease associations for the complete GenAge Build 21 human-gene catalogue. On 19 September 2026, the user explicitly chose **all 307 human genes**, regardless of whether they have orthologs or longevity evidence in model organisms. Pedro's supplied wording proposed cautious human extrapolation and possible disease associations, but did not specify this exact complete-catalogue design. The user's direct scope decision governs the analysis.

## Status

Ready for review. All 307 catalogue genes are retained; 306 current Ensembl targets were queried against direct Open Targets Platform 26.06 associations and MT-CO1 is explicitly marked unqueryable through this route. The checksummed source cache, 422,250 unique target–term relations, evidence datatype/source scores, all-gene summary, independent verification and four-panel review figure are complete.

## Proposed outcome

A resource-bounded map of disease associations for every GenAge human gene, including explicit zero/unmapped states and evidence-type provenance. LongevityMap membership and association direction are annotations, not inclusion filters. Summary figures may describe database evidence but must not imply that disease association establishes an effect on human lifespan.

## Inputs

| Input | Version or date | Role | Integrity |
|---|---|---|---|
| `../../GenAge Project/human_genes/genage_human.csv` | Build 21, 28 August 2023 | Complete primary gene universe: all 307 entries | Raw/immutable; SHA-256 `3515935b6c4d40cbccf35e129aa52206cc95c1fe2b426adec784c2e41927174d` |
| `../../Other downloads from HAGR/longevity_genes/longevity.csv` | LongevityMap Build 3, 24 June 2017 | Separate direct-human-longevity evidence annotation | Raw/immutable; 550 study-level records; SHA-256 `1202393b6bc02c63115ef5b39fb7e438ef7252a9dbaf8721e101a2893abf013c` |
| Open Targets Platform GraphQL API | Data 26.06; API 26.6.3; retrieved 19 September 2026 | Direct target–disease/phenotype relations plus overall, datatype and datasource scores | `enableIndirect: false`; responses and request metadata cached and checksummed under `results/source_cache/opentargets_26_06/` |

## Method and reproduction

Begin with all 307 GenAge genes; do not filter on model-organism orthology or LongevityMap membership. `scripts/audit_human_genes.py` validates current identifiers through NCBI Datasets without rewriting the immutable catalogue and separately matches LongevityMap study records. `scripts/collect_opentargets_associations.py` freezes Open Targets 26.06, retrieves direct associations in resumable pages, preserves overall/datatype/datasource scores and reconciles offset-page overlap caused by tied low scores. Terms are classified with ontology namespace plus Open Targets top-level category while all original term IDs, names and category fields remain in the full table. `scripts/verify_opentargets_outputs.py` reconstructs every pair from checksummed cached responses. `scripts/make_review_figure.py` generates the source-data tables and all figure exports in Python.

Reproduce and verify from the project root:

```powershell
python "task/task 5/scripts/audit_human_genes.py"
python "task/task 5/scripts/collect_opentargets_associations.py"
python "task/task 5/scripts/verify_opentargets_outputs.py"
python "task/task 5/scripts/make_review_figure.py"
```

Task 4 mappings are outside the inclusion logic. Disease association does not establish an effect on human lifespan. “All relations” means all qualifying direct associations in the declared source release, not all claims in the biomedical literature. Descriptive disease mapping and statistical enrichment are different deliverables; no enrichment is implied by this scope.

## Canonical outputs

| Output | Contents |
|---|---|
| `results/human_gene_audit.csv` | All 307 GenAge genes with frozen and current identifiers; four symbol updates are preserved rather than overwritten |
| `results/longevitymap_matches.csv` | 408 matched study-level mentions covering 126 GenAge genes; annotation only, not a filter |
| `results/opentargets_gene_disease.csv` | 422,250 unique direct gene–term relations with original term/category fields and overall score |
| `results/opentargets_datatype_scores.csv` | 454,625 relation–datatype score rows |
| `results/opentargets_datasource_scores.csv` | 456,430 relation–datasource score rows |
| `results/gene_disease_summary.csv` | One row for each of 307 genes, including MT-CO1's unqueryable state and term-class counts |
| `results/opentargets_collection_summary.json` | Frozen release metadata, output hashes, counts and pagination-reconciliation audit |
| `results/verification_report.json` | Independent PASS report reconstructed from checksummed cached responses |
| `figures/task5_human_gene_disease_review.{svg,pdf,tiff,png}` | Checked double-column review figure; SVG text remains editable |
| `results/figure_1_source_data/` | Four panel-specific source-data tables |
| `notes/figure_1_legend.md` | Review caption and interpretation boundary |

## Verification

`scripts/verify_opentargets_outputs.py` passed on 19 September 2026. It independently confirmed all 307 gene rows, reconstructed all 306 queryable targets from checksummed cache files, preserved one explicit unqueryable gene, matched all 422,250 unique output pairs to the cached source union, checked per-gene/class totals, validated all evidence-score foreign keys and confirmed output hashes. Thirty-seven genes required pagination reconciliation; every one reached the API-declared unique count. Visual QA confirmed readable PNG/TIFF output and editable SVG text.

## Descriptive result

Of 422,250 direct relations, 393,383 were classified as disease/clinical condition, 25,715 as phenotype/trait, 2,376 as biological process, 422 as medical procedure and 354 as non-human/other. Every one of the 306 queryable genes had at least one returned relation, ranging from 40 for GTF2H2 to 6,459 for EGFR; this breadth makes raw counts poor evidence of specificity. The median Open Targets overall score was 0.0326, the 90th percentile was 0.0660 and 1,383 returned relations had score zero. Literature and animal-model evidence were the two most prevalent datatypes.

## Decisions and limitations

- Keep this human disease-association question separate from model-species orthology (Task 4) and the closed enrichment workstream (Task 6).
- GenAge catalogue membership is broader than direct human-longevity association; LongevityMap evidence must remain a separate field.
- The reported “few samples” concern applies to direct human-longevity evidence, not to the fixed 307-gene catalogue denominator.
- Open Targets “direct” means no disease-ontology propagation in this query; it does not mean direct causal human evidence. The evidence includes literature and animal-model data.
- Raw relation counts reflect study density, evidence-source breadth and ontology structure as well as biology; they must not be ranked as ageing importance.
- MT-CO1 lacks a current Ensembl target in the NCBI mapping used here. Its status is unqueryable, not association-free.
- Database association scores are not probabilities, causal effect sizes or measures of longevity relevance.

## History

See [CHANGELOG.md](CHANGELOG.md). This workspace has no Git history or generating commit.

## Contents

This folder contains task-local scripts, immutable-response caches, derived tables, figure source data, review exports, captions and verification. Shared HAGR files remain unchanged.
