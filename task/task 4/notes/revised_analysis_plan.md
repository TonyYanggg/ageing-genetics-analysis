# Task 4 analysis plan — all genes in the biological figure style

Updated 2026-09-26: the user requested the main-figure style for all genes. The primary visual deliverable is the complete four-column atlas, including all 320 genes and 206 source pairs. The six-family experiment illustration and complete count overview are companions; the all-nine audit remains source/coverage provenance. No breadth or selected-family cutoff applies to the atlas.

## Question and scope

Which orthologous genes have lifespan evidence across budding yeast, worm, fly and mouse, and what do the experiments show? These organisms contain 2,138/2,205 Build 21 model records. The atlas displays every gene participating in the 206 qualifying direct pairs, not every GenAge record regardless of orthology evidence. Existing primary summaries are preserved for 20 genes; 300 other cards display explicitly labelled catalogue annotations and references. The other five species remain in the coverage audit. Human genes and disease associations are outside Task 4.

Maintain three separate evidentiary levels: an explicit orthology relationship, recurring lifespan association in GenAge, and comparable causal function. Only the first two are supported by the complete pair analysis. Do not treat overlapping pairs as source-defined orthogroups or independent conserved genes.

## Implemented sequence

1. Start with immutable GenAge Build 21 and the frozen OMA All.May2026 explicit protein-pair analysis. Keep the all-nine source coverage, unresolved mapping and record ledgers in `results/oma_pairwise/refined/`; OMA covers seven species.
2. Filter the 228 source-supported direct pairs to those with both genes among the principal four. Preserve all 206 qualifying rows, the original OMA gene relationships, GenAge IDs, lifespan-effect labels and direction classifications. Report distinct genes and partner multiplicity separately.
3. Retain mapping coverage, species-pair counts and direction categories in the supporting overview. Screen all 206 pairs into connected candidate networks for search only: 133 networks, 22 spanning at least three species. Never infer missing edges or source-defined orthogroups by connectivity.
4. Preserve the approved TOR, S6 kinase and insulin/IGF receptor evidence. Add AMPK, ATG7 and SOD2 after targeted primary-paper review. The [selection rationale](expanded_selection.md) documents the six choices and eligible alternatives; `expanded_additional_evidence.tsv` records new experiments. The combined ledger retains original GenAge annotations separately from displayed experiment outcomes. Catalogue-wide results remain database associations.
5. Use `build_all_gene_atlas.py` to paginate all 133 candidate networks in the same 183 × 247 mm four-column card style. Include all two-species networks and all paralogs; write every explicit direct pair beneath its network. Export the full bookmarked PDF, all individual PNG/SVG/PDF/600-dpi TIFF pages, a searchable gene index and exact card/pair/reference tables. Use `verify_all_gene_atlas.py` to reconcile every gene, pair and page, then inspect contact sheets and representative dense full-size pages.
6. Use the [atlas legend](all_gene_atlas_legend.md) for all-gene inclusion, card types, pair codes, reference provenance and limitations. The [expanded legend](expanded_figure_legend.md) retains experiment detail: missing direct edges, mouse sex specificity, Igf1r replication conflict, S6K1 correction, AMPK assay differences, worm ATG7 eat-2 background and mouse Sod2 tissue specificity. No selected direction is generalized into a universal pro- or anti-ageing function.

## Acceptance and handoff

The all-gene objective is met only when the actual rendered atlas contains every one of the 320 gene keys and 206 source pairs, including all two-species cases, with no gene omitted or pair added. Every card must retain evidence provenance and clearly distinguish reviewed experiments from catalogue-only evidence. Every page and export must pass text/content/geometry checks, with a visual sweep of all pages. The searchable index, source tables, legend and current project instructions must point to the atlas. Existing complete-analysis and example verifiers must also pass. Pedro has not reviewed the atlas; it does not claim uniform conserved function or primary validation of all 300 additional genes.

Historical scope correction and the former all-nine plan are documented in `../CHANGELOG.md`, `scope_revision.md` and `refined_review.md`.
