# Task 5 analysis scope and verification

## Scope

- Include all 307 genes in `GenAge Project/human_genes/genage_human.csv`, regardless of model-organism orthology or Task 4 membership.
- Do not use Task 4 mappings as an inclusion filter. Task 4 may be joined later only as optional, clearly labelled context on explicit user instruction.
- Treat LongevityMap as a separate annotation of direct human-longevity studies. Do not equate GenAge catalogue membership with direct human-longevity association.
- Map diseases and phenotypes using a declared, versioned resource. “All” means all qualifying associations in that resource and release, not all biomedical knowledge.

## Evidence rules

- Preserve source identifiers, disease ontology IDs, evidence data type/source, scores and retrieval provenance wherever available.
- Prefer direct gene–disease associations for the main table. Do not silently mix ontology-propagated indirect associations with direct evidence.
- Separate diseases from measurements, normal traits and generic phenotypes where the source ontology supports the distinction. If classification is imperfect, keep the original term and state the limitation.
- Do not interpret association scores as probabilities, effect sizes or causal proof. Disease association does not establish an effect on human lifespan.
- Keep every GenAge gene in the audit, including genes with zero associations, retired identifiers, ambiguous mappings or unsupported target types.
- Avoid arbitrary score thresholds in the complete source table. Any presentation threshold must be labelled and accompanied by the full table.

## Outputs and verification

- Keep scripts, cached responses, derived tables, figures, notes and checks in Task 5.
- Preserve shared HAGR inputs unchanged and record their hashes.
- Verify gene-universe reconciliation, exact identifier mappings, pagination/completeness, duplicate collapse, ontology classification, zero-association genes and every plotted value independently.
- Do not mark the task complete until the source scope, full gene–association table, all-gene summary, limitations, figure/caption and reproducible verification exist.
