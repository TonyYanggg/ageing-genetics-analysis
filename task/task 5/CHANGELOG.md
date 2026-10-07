# Change history

## 2026-09-19 — Verified Open Targets review package

- Audited all 307 GenAge human genes through current Entrez records; 306 mapped to an Ensembl target and MT-CO1 remained explicitly unqueryable through that route.
- Matched 408 LongevityMap study-level gene mentions covering 126 catalogue genes without using them as an inclusion filter.
- Froze Open Targets Platform data 26.06/API 26.6.3 and collected all direct associations for the 306 queryable genes.
- Reconciled unstable offset pagination for 37 genes by unioning checksummed alternative-page retrievals until every target reached its API-declared unique count.
- Produced 422,250 unique gene–term relations, 454,625 datatype-score rows, 456,430 datasource-score rows and a 307-row all-gene summary.
- Added an independent verifier; all gene, pair, class, foreign-key and hash checks passed.
- Added a Python-generated four-panel review figure in SVG/PDF/TIFF/PNG, panel source data, figure contract, caption and visual QA.

## 2026-09-19 — Complete human-catalogue scope

- User selected all 307 GenAge Build 21 human genes as the Task 5 universe, regardless of model-organism orthology.
- LongevityMap was retained as a separate human-longevity evidence annotation, not an inclusion filter.
- Defined disease mapping as resource-bounded and evidence-preserving; genes with no qualifying association remain in the denominator.
- Began disease-source and identifier audit. No disease findings or figure are yet accepted.

## 2026-09-18

- Qualified attribution using the referenced conversation; original email remains unverified.
- Clarified the later exploratory disease-association scope and distinguished ageing catalogue membership from direct human longevity evidence.
- Documented candidate inputs, missing disease source, planned method and absence of scientific outputs.
