# Task 4 — Longevity orthologs across four major model organisms

## Request and decision

The current question is: **which orthologous genes have lifespan evidence across budding yeast, nematode, fruit fly and mouse, and what do the actual experiments show?** These four organisms contain 2,138/2,205 GenAge Build 21 model records (97.0%). On 26 September 2026 the user requested the same main-figure style for all genes. The primary visual output is therefore a paginated four-column atlas of every gene in the qualifying pair set: 320 distinct genes and all 206 direct pairs. The six-family figure remains a concise experiment-level companion. This does not include all 2,138 records irrespective of ortholog evidence, or imply that the other five organisms lack relevant genes. Task 4 has **no human-gene component**.

## Status

**Review-ready and independently verified for the all-gene visual atlas; not yet approved by Pedro.**

The [complete visual atlas](figures/all_gene_atlas/all_gene_ortholog_atlas.pdf) uses the same boxed gene-card style across 36 pages, with four organism columns on every page. All 133 candidate networks are included, including all two-species cases. Each of 320 genes appears once; multiple paralogs are stacked rather than omitted. The 20 previously checked experiments retain their summaries; 300 other cards explicitly show catalogue annotations and linked references. The [legend](notes/all_gene_atlas_legend.md) explains this distinction. A [searchable local index](figures/all_gene_atlas/index.html), [gene/card index](results/all_gene_atlas/gene_cards.csv), [exact displayed pairs](results/all_gene_atlas/all_displayed_pairs.csv), [row/page index](results/all_gene_atlas/page_index.csv) and [reference provenance](results/all_gene_atlas/catalogue_references.csv) accompany the PDF. Every page has PNG, editable SVG, PDF and 600-dpi TIFF exports at 183 × 247 mm. The 320 gene keys represent 321 GenAge records because one worm gene has two entries.

The [six-family companion figure](figures/expanded/expanded_ortholog_examples.png) ([PDF](figures/expanded/expanded_ortholog_examples.pdf), editable SVG and 600-dpi TIFF in the same folder) contains 20 genes and 21 selected direct pairs. Its [legend](notes/expanded_figure_legend.md), [selection rationale](notes/expanded_selection.md), [evidence table](results/expanded/displayed_evidence.csv), [selected pairs](results/expanded/selected_oma_pairs.csv) and [candidate screen](results/expanded/candidate_screen.csv) remain available. It is not the exhaustive visualization requested on 26 September. Pedro's positive response concerned the original three-example version, not this expansion or the complete atlas.

The [complete overview](figures/complete_four_species/complete_ortholog_overview.png) and [legend](notes/complete_figure_legend.md) now serve as supporting analysis. They summarize all 206 qualifying direct pairs among 320 distinct model genes. The [complete sortable pair table](results/complete_four_species/all_qualifying_pairs.csv), [six-comparison source table](results/complete_four_species/species_pair_summary.csv), [gene multiplicity table](results/complete_four_species/gene_pair_multiplicity.csv), [all-nine coverage table](results/complete_four_species/all_nine_coverage.csv) and [summary](results/complete_four_species/summary.json) retain the complete catalogue and source limits.

The [three-example illustration](figures/focused/focused_ortholog_examples.png) remains supporting context for TOR, S6 kinase and insulin/IGF receptors. Its [paper-level legend](notes/focused_figure_legend.md), [ten-record evidence ledger](notes/focused_experiment_evidence.tsv) and [ten selected OMA pairs](results/focused/selected_oma_pairs.csv) retain targeted primary-experiment checks. It is illustrative, not the complete gene-level result.

The broader [refined OMA analysis](results/oma_pairwise/refined/README.md) remains source provenance: 228 explicit pairs across seven OMA-represented organisms, of which 206 have both genes in the principal four. Its all-nine inventory and gap ledger remain available. The prior broad [nine-organism pairwise figure](figures/oma_refined/pairwise_evidence.png) is supporting history. Pair counts overlap in genes and do not count independent orthogroups. The earlier human-anchor and OrthoDB screens are historical audits, not current results.

## Complete result and interpretation

The six organism-pair counts are 97 yeast–worm, 34 yeast–fly, 14 yeast–mouse, 33 worm–fly, 16 worm–mouse and 12 fly–mouse. Of the 206 direct pairs, 66 are Anti/Anti, 34 Pro/Pro, 32 discordant or mixed and 74 direction-unresolved at the catalogue level. Sixty-seven of the 320 genes have multiple partners, including 47 with partners in more than one species. The complete table identifies each gene, GenAge record, lifespan-effect annotation, direction label, direct OMA relationship and source. The [independent verifier](scripts/verify_complete_four_species.py) reconstructs counts, gene identities, coverage, multiplicity and figure values.

## Output provenance and contents

| Output | Status | Generated by / inputs | Role |
|---|---|---|---|
| `figures/all_gene_atlas/` | Current | `scripts/build_all_gene_atlas.py`; complete pair table, expanded evidence/screen and frozen Task 3 reference mapping | All 320 genes in the main-figure style; bookmarked PDF, all page exports and searchable index |
| `results/all_gene_atlas/`, `notes/all_gene_atlas_legend.md` | Supporting | Same builder plus manually written legend | All gene cards, all source pairs, page index, references, source hashes and QA |
| `figures/expanded/expanded_ortholog_examples.*` | Supporting | `scripts/build_expanded_examples.py`; complete pair table, raw GenAge CSV and the two evidence ledgers | Concise six-family biological companion |
| `results/expanded/displayed_evidence.csv`, `selected_oma_pairs.csv`, `candidate_screen.csv`, `summary.json` | Supporting | Same builder; source paths and hashes in summary | Exact display provenance and complete candidate screen |
| `notes/expanded_figure_legend.md`, `expanded_selection.md`, `expanded_additional_evidence.tsv` | Supporting | Manual primary-paper checks with explicit review depths | Legend, selection rationale and new evidence |
| `results/expanded/papers/` | Supporting | `scripts/collect_expanded_papers.py`; Europe PMC publication API | Cached metadata/abstracts, including considered alternatives |
| `figures/complete_four_species/`, `results/complete_four_species/` | Supporting | `scripts/build_complete_four_species.py`; frozen refined OMA pairs, GenAge and coverage ledger | Complete 206-pair catalogue and overview |
| `figures/focused/`, `results/focused/` | Supporting | `scripts/plot_focused_examples.py`; original evidence ledger and frozen pairs | Preserved three-row template Pedro liked |

`scripts/` holds reproducible builders/collectors/verifiers, `results/` derived evidence and audits, `figures/` exports, and `notes/` methods, legends and evidence decisions. Earlier analysis stages remain labelled historical provenance below.

## Selected biological companion and interpretation

| Example | Organisms shown | Direct OMA pairs among shown genes | Important qualification |
|---|---|---:|---|
| TOR | Yeast, worm, fly, mouse | 6 | Yeast TOR1 edges are one-to-many; yeast replicative lifespan and animal survival differ |
| S6 kinase | Worm, fly, mouse | 2 | No direct worm–mouse pair; mouse S6K1 benefit is female-specific |
| Insulin/IGF receptor | Worm, fly, mouse | 2 | Mouse Igf1r is one of two OMA-linked daf-2 receptor paralogs; no direct fly–mouse pair in the selection |
| AMPK | Yeast, worm, fly | 2 | Selected yeast replicative result differs from the animal experiments; no direct yeast–worm pair in the selection |
| ATG7 | Yeast, worm, fly | 3 | Worm cell tests the longevity benefit of dietary restriction in eat-2 mutants |
| SOD2 | Yeast, worm, fly, mouse | 6 | Worm deletion extends lifespan; yeast sch9Δ benefit is lost; fly and connective-tissue-deficient mouse survival decreases |

Each occupied cell displays a gene, intervention, outcome and primary PMID; its GenAge ID and review depth are in the evidence ledger. Selection balances interpretability, primary traceability and biological breadth; it is not a top-rank or prevalence estimate. Arrows and colours encode reported intervention effects, not GenAge Pro/Anti labels. Worm atg-7's original Increase/Unclear annotation remains in the table alongside the selected eat-2 experiment. The AMPK and SOD2 rows deliberately retain differing results. TOR and S6K share a pathway and fly paper. Empty cells do not mean no ortholog or no lifespan effect.

For a manuscript, the bounded takeaway is: “Selected orthologous genes have experimental lifespan evidence across major model organisms, spanning nutrient signalling, autophagy and oxidative stress, but intervention and organism context can change the direction and interpretation of the effect.” The 206-pair catalogue documents the broader database evidence and is not a count of conserved causal mechanisms.

## Inputs and provenance

- Immutable `../../GenAge Project/models_genes/genage_models.csv`: GenAge Build 21, 28 August 2023, 2,205 model-organism entries across nine species.
- Frozen OMA All.May2026 explicit model-to-model protein pairs and identity resolution in `results/oma_pairwise/`. The [refined source README](results/oma_pairwise/refined/README.md) records release, caches, hashes, coverage and exclusions. Human genes are not inputs.
- Selected primary papers are linked in the [expanded legend](notes/expanded_figure_legend.md); the original ten records remain in the focused ledger. New manual evidence is in `notes/expanded_additional_evidence.tsv`; structured publication records are cached in `results/expanded/papers/`. Most checks are abstract/primary-record level; mouse S6K1 and yeast/worm AMPK have the documented targeted full-text checks. This is not a systematic study-quality or replication assessment.
- Catalogue-linked references for all atlas genes come from `../task 3/results/verified_timeline/entry_references.csv`, the existing frozen GenAge/LibAge mapping. The atlas retains source URLs, hashes and timestamps in `results/all_gene_atlas/catalogue_references.csv`. These are not new primary-paper checks for the 300 catalogue-only cards.

The all-nine inventory includes 911 yeast, 889 worm, 202 fly, 136 mouse, 61 fission-yeast, three Podospora, and one each *C. briggsae*, zebrafish and hamster entry. Sparse organisms and OMA source gaps are recorded in [coverage.csv](results/oma_pairwise/refined/coverage.csv) and the [resolution ledger](results/oma_pairwise/refined/resolution_ledger.csv). A missing pair is not a biological negative.

## Reproduce and verify

From the project root, with the existing Python, matplotlib and Pillow environment:

```powershell
python -X utf8 'task/task 4/scripts/build_all_gene_atlas.py'
python -X utf8 'task/task 4/scripts/verify_all_gene_atlas.py'
python -X utf8 'task/task 4/scripts/build_expanded_examples.py'
python -X utf8 'task/task 4/scripts/verify_expanded_examples.py'
python 'task/task 4/scripts/build_complete_four_species.py'
python 'task/task 4/scripts/verify_complete_four_species.py'
python 'task/task 4/scripts/plot_focused_examples.py'
python 'task/task 4/scripts/verify_focused_examples.py'
```

The all-gene atlas verifier passed on 26 September 2026: all 320 distinct gene keys, 321 GenAge IDs, 206 exact pairs and 133 networks match the frozen inputs. It checked the actual SVG/PDF text on all 36 pages, decoded every PNG/TIFF, checked physical dimensions and 600-dpi TIFF metadata, verified all source hashes, and reconciled all references and page-index entries. The builder checked 3,558 text elements for containment and overlap. A final clean rebuild from the frozen inputs was followed by the same full verification suite and visual inspection of contact sheets, dense/full-size PNG pages and PDF renders; the dual-name worm gene wraps without clipping. Search was exercised with a gene symbol, a GenAge ID, no match and an empty query restoring all 320 entries. See the [QA record](notes/all_gene_atlas_qa.md) and [machine verification](results/all_gene_atlas/verification.json). No Git generating commit exists.

The expanded builder uses the complete pair table, raw GenAge records, preserved original evidence and new evidence ledger. It generates the candidate screen, selected pairs, combined evidence, summary and four figure exports without network access. Its verifier independently reconstructs candidates by taxon/Entrez key and checks exact selected pairs, original evidence preservation, gene labels, primary PMIDs and exports. `collect_expanded_papers.py` is the optional network collector; cached records suffice for offline verification. The complete and focused verifiers also pass. PNG and Python-rendered PDF were visually inspected; the builder checks all 117 text elements for box containment and overlap. Existing [OMA reproduction and audit instructions](results/oma_pairwise/refined/README.md) remain available. The project is not a Git repository; no commit or branch exists.

2026-09-24 layout revision: the same text and data now sit inside aligned title, grid, edge-note and caveat boxes. The width remains 183 mm; height is 166 mm to preserve readable type. The plot script checks text containment and overlap before writing exports. No evidence, annotation, reference or interpretation changed.

## Limitations and handoff

The catalogue and its new visual atlas are complete for **recoverable direct OMA protein pairs among the principal four**; they are not a universal orthology atlas or an effect-size meta-analysis. All 320 participating genes are displayed, but 300 cards remain catalogue annotations without newly reviewed intervention summaries. The 67 records outside the four-organism comparison remain in the all-nine coverage table; Podospora and hamster lack OMA source representation. Missing mappings and pairs are not biological negatives. Network rows are layout units, not OMA orthogroups. The reviewed cards retain sex specificity, conflicting Igf1r cohorts, the S6K1 correction and intervention/background qualifications. Pedro has not reviewed the atlas. Manuscript presentation remains an editorial decision.

## History and retained work

[CHANGELOG.md](CHANGELOG.md) records the scope changes. [analysis scope](notes/analysis_scope.md) and [revised analysis plan](notes/revised_analysis_plan.md) describe the present direction. The [refined review](notes/refined_review.md) documents broader provenance and experimental caveats. Earlier human-anchor outputs under `results/` and `figures/cross_species_homologs.*` are preserved for reproducibility but must not be used as evidence for the present model-to-model question. Task 5 (human genes and disease associations) remains separate; Task 6 is closed unless reopened by the user.
