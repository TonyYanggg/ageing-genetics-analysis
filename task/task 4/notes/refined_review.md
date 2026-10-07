# Refined model-to-model review package

**Status after Pedro's scope decision (2026-09-24):** this document describes the preserved broad OMA audit and its former Figures 1–2. The [focused three-example figure and legend](focused_figure_legend.md) are now the main Task 4 review handoff. The numerical audit below remains valid supporting provenance, not the review's headline presentation.

## Claim and current boundary

GenAge Build 21 contains repeated lifespan associations among **explicit protein-ortholog pairs**. This is a catalogue-based, source-specific result, not an exhaustive list of conserved longevity mechanisms. All nine organisms enter the audit; OMA represents seven. There are 228 unique lifespan-evidence pairs, including 145 pairs with Pro/Anti annotations on both sides. These overlapping pairs must not be described as 228 independent genes or orthogroups.

## Identity recovery and exclusions

The frozen direct-NCBI-xref baseline has 219 pairs. The refined analysis uses exact UniProt, Ensembl, model-organism authority identifiers (including PomBase), and yeast systematic locus tags supplied by NCBI Datasets. OMA substring search results are accepted only on literal identifier equality and matching taxon. No gene-symbol guesses, human anchors, or transitive grouping are used. Sixty-three additional gene identities are recovered; this yields 2,136 mapped unique species+Entrez keys and 2,138 mapped records. The extra mappings add nine evidence pairs, not 63 pairs.

All 2,205 records remain in the resolution ledger. The 67 remaining records comprise 18 ambiguous gene-to-protein links, 18 noncoding loci, eight current pseudogene annotations, eight discontinued identifiers, three missing IDs, four confirmed heterologous constructs, two unresolved construct donors, two gene-paper conflicts, one taxon conflict, two OMA species-gap records and one remaining protein identifier gap. Categories are mutually exclusive record-level triage, not counts of independent biological genes. Current gene-type metadata was collected for baseline unmapped records only. RNA orthology is not assessed by this protein pipeline.

## Sparse-organism checks

- C. briggsae daf-2 contributes four explicit OMA pairs. Its cited abstract establishes orthology/dauer results, but does not independently validate the catalogue's lifespan assay details.
- Zebrafish tert has an explicit pair with mouse Tert. Both have GenAge Pro-Longevity annotations; this is not a new mechanistic validation.
- Hamster Csnk1e is absent from OMA. Separate Ensembl 116 calls identify mouse and zebrafish orthologs, neither with a matching GenAge entry. Broader hamster comparisons remain unassessed, not negative. Historical and current Ensembl IDs differ; do not mix assemblies silently.
- Podospora is absent from OMA. Exact lookup recovers Entrez 6194500 as PODANS_6_5980 (translation initiation factor), but its cited lifespan paper concerns PaAif2/PaAmid2. Gene-paper attribution is quarantined. The other conflicting record labels cox1 while its paper studies cyc1. The remaining COX1 paper concerns removal of a mitochondrial intron, not loss of the encoded protein. Ensembl fungi calls for the recovered locus return no orthologs to the two queried yeasts; this is not proof of biological absence.

## Experimental review

The manually curated `experiment_review.tsv` records the review depth, manipulation, context, source PMID and claim limit for each case. Primary metadata/abstracts and three openly licensed JATS main articles are cached. Most cases are abstract-level; targeted full-text review is explicitly labelled. Supporting information has not been reviewed. Selection is purposive: sparse organisms, identity conflicts, and a native S6-kinase example, not an unbiased systematic literature review.

The S6-kinase example illustrates why matching Anti-Longevity labels need qualification: worm inhibition, fly dominant-negative expression and mouse deletion are different interventions. Mouse survival benefit is female-specific in the retrieved results. The [2011 correction text](https://citeseerx.ist.psu.edu/document?doi=f39eea152043b5246fdfa4632451e33eaa76436f&repid=rep1&type=pdf) reports that a worm *rsks-1;aak-2* double-mutant strain lacked the *rsks-1* allele. Replacement experiments confirmed AMPK dependence of the worm lifespan effect, but growth and reproduction effects were AMPK-independent. It does not report a change to mouse survival. Sources: [worm study](https://pubmed.ncbi.nlm.nih.gov/17266680/), [second worm study](https://pubmed.ncbi.nlm.nih.gov/17266679/), [fly study](https://pubmed.ncbi.nlm.nih.gov/15186745/), [mouse primary results](https://pmc.ncbi.nlm.nih.gov/articles/PMC4954603/).

## Figure 1 contract and caption

Conclusion: explicit protein orthologs have recurring catalogue lifespan evidence, while coverage and evidence are uneven across all nine organisms. Quantitative-grid archetype; Python/matplotlib; 183 × 145 mm; editable SVG/PDF, PNG preview, 600-dpi TIFF. Descriptive catalogue counts, not biological replicates: no inferential test, confidence interval or error bars are appropriate.

**Figure 1 | Source-specific ortholog pairs with lifespan evidence across GenAge model organisms.** (a) Mapped unique species+Entrez gene keys divided by all identifiable keys per organism, after ambiguous mappings are excluded. Stable-ID recovery is included; missing-ID records remain in the audit but not this denominator. (b) Unique explicit OMA gene pairs with an Increase, Decrease, or Increase and Decrease lifespan annotation on both sides. P/A labels do not determine inclusion. The symmetric matrix displays each pair twice. Zero means no qualifying pair was recovered by this workflow, not absence of an ortholog or lifespan effect. NA indicates organisms not represented in OMA. All nine organisms remain in the inventory; supplementary Ensembl checks are not pooled into this matrix. Pairs share genes and are not independent orthogroups. Source: GenAge Build 21 and OMA All.May2026, retrieved 18 September 2026. Source tables: `refined/coverage.csv`, `refined/figure_pair_counts.csv`, `refined/longevity_pairs.csv`.

## Figure 2 contract and caption

Conclusion: actual model-gene pairs illustrate annotation agreement, discordance and contributions from sparse organisms, without turning pair connectivity into orthogroups. Qualitative table, Python/matplotlib, 183 × 153 mm; editable SVG/PDF, PNG and 600-dpi TIFF. No pooled experimental effect sizes or inferential statistics.

**Figure 2 | Selected model-gene pairs with GenAge longevity annotations.** Ten explicit edges illustrate S6-kinase annotation agreement, SOD annotation discordance, and contributions from C. briggsae and zebrafish. Selection is purposive, not a ranking or independent validation sample: all qualifying pairs among GenAge IDs 578/583/1737 and 984/621/622, plus all pairs incident to IDs 2025 or 1792, are displayed. Relation cardinality follows the source and is oriented from model gene 1 to model gene 2. Anti/Anti and Pro/Pro refer to catalogue influence labels, not identical experimental interventions. No worm–mouse S6-kinase pair was recovered; no edge is added through transitivity. Native S6-kinase experiments differ in manipulation and context; the mouse result is sex-specific. The correction to the linked paper concerns a worm double-mutant strain, as described above. SOD and sparse-organism rows are catalogue-level examples, not fully reviewed experimental comparisons. Figure 1 retains the all-nine inventory and source-gap states. Source data: `refined/figure_example_pairs.csv`; experimental context and citations: `experiment_review.tsv` and `experiment_review_sources.csv`.

## Acceptance and residual limitations

- Both figures exist with exact source tables; example selection, all 228 direction classifications, record ledger and review PMID provenance are independently checked.
- SVG text remains editable; PDFs contain embedded TrueType fonts and selectable text; TIFFs are 600 dpi at declared dimensions. PNGs are visually inspected and PDF bounds/renders are checked separately.
- Final offline reproduction passed on 18 September 2026: network socket connections were disabled, then identity audit, exact-ID recovery, paper metadata/full-text cache reuse, refined analysis, resolution ledger, both plots and both verifiers were run in sequence. The same 228 pairs, 145 directional pairs, 67 remaining records and ten example pairs were recovered. All 216 cached payload hashes and the immutable GenAge hash pass. Baseline verification still reconstructs 219/136 pairs. Python compilation passes.
- The S6-kinase correction was located and reviewed on 25 September 2026; see the experimental review above. The prior access limitation is retained in the earlier change history, not as a current open issue.
- Unsupported-species and noncoding coverage remain limitations; no unsupported absence or universal conservation claim is made. Full primary-paper re-curation of every pair would be a further research stage, not a consequence of successful numerical checks.

This is a bounded review package, not proof of biological perfection or exhaustive all-nine coverage. The associated task status records whether final QA has completed.
