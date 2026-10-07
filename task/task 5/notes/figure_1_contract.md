# Task 5 review figure contract

## Core conclusion

The complete GenAge human catalogue is broadly connected to Open Targets disease evidence, but association density and evidence composition vary sharply across genes; database coverage is therefore descriptive and must not be interpreted as strength of human-longevity evidence.

## Figure archetype

Asymmetric quantitative grid with the all-gene ranked distribution as the hero panel.

## Target and output

- Purpose: scientific review with Pedro before manuscript-level interpretation.
- Backend: Python (the task's analysis pipeline is Python-specific).
- Final size: double-column, 183 mm wide by approximately 155 mm high.
- Exports: editable SVG, editable-text PDF, 600 dpi TIFF and 300 dpi PNG preview.

## Panel map and evidence hierarchy

- **a — hero evidence:** every gene ranked by number of unique direct Open Targets relations on a logarithmic scale, distinguishing genes with and without separate LongevityMap evidence. Shows uneven database density without filtering on longevity evidence.
- **b — scope and composition:** 307 GenAge genes, 306 queryable targets, one explicitly unqueryable mitochondrial gene, and the term-class composition of all unique direct relations. Establishes the denominator and prevents silent gene loss.
- **c — score distribution:** reverse cumulative distribution of the Open Targets overall association score, with median and 90th percentile marked. Shows that the score distribution is broad and avoids an arbitrary high/low cutoff.
- **d — evidence composition:** number of unique gene–disease/phenotype pairs and genes represented by each Open Targets evidence datatype. Shows that “association” combines heterogeneous evidence.

## Statistics and source data

- No inferential test is claimed or needed; this is a complete resource-bounded catalogue, not a sample comparison.
- Counts are unique direct target–term pairs in Open Targets Platform release 26.06.
- Score summaries are descriptive quantiles; Open Targets scores are not probabilities or effect sizes.
- Panel source tables are generated alongside the figure and remain traceable to the verified association and score tables.

## Reviewer risks

- Raw relation counts are strongly affected by study density and ontology breadth.
- Open Targets' disease endpoint includes phenotype, process, procedure and non-human terms; these must be separated visibly.
- LongevityMap overlap is an annotation, not an inclusion criterion or validation label.
- MT-CO1 is unqueryable by this Ensembl-target route; zero must not be reported as no disease association.
- Direct database association does not imply causality or an effect on human lifespan.

## Image integrity

No microscopy or manipulated raster data are used. All panels are generated directly from tabular outputs by the task-local Python script.
