# Task 5 analysis contract

## Question

What disease and phenotype associations are reported for every gene in the GenAge Build 21 human catalogue, and what evidence supports each association?

## Universe

All 307 GenAge human entries. Inclusion does not depend on model-organism orthology, Task 4 results or LongevityMap membership. Entrez Gene ID is the frozen catalogue key; current Ensembl/approved symbols are mapped and audited fields.

## Evidence layers

1. GenAge catalogue membership: ageing relevance curated by HAGR.
2. LongevityMap: study-level human longevity associations, including significant, nonsignificant and population-specific results.
3. Disease resource: versioned direct gene–disease/phenotype associations with evidence-source provenance.

These layers answer different questions and must not be collapsed into one “human longevity gene” label.

## Counting

- Primary row: unique GenAge gene plus disease ontology ID, retaining original disease name and evidence details.
- Deduplicate repeated source rows only with traceable provenance.
- Count every gene once in all-gene denominators, including zero-association and mapping-failure states.
- Report diseases, phenotypes and measurements separately when classification is available.
- Scores describe a database's evidence aggregation and are not comparable to effect sizes.

## Claim boundary

The analysis can state that a database reports evidence connecting a GenAge gene to a disease or phenotype. It cannot by itself establish causality, direction of effect, therapeutic benefit, ageing mechanism or an effect on human lifespan.
