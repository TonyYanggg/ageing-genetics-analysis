# Task 4 analysis and figure contract

Status: historical contract for the preserved **exploratory human-anchor screen only**. The main Task 4 design is now [revised_analysis_plan.md](revised_analysis_plan.md); see [scope_revision.md](scope_revision.md) for attribution and scope corrections. Running the existing scripts reproduces this screen, not the revised task.

Date: 2026-09-18. Backend: Python, selected by the assistant at the user's request.

## Question and scope

Which human gene anchors link GenAge longevity-associated model genes in at least two of budding yeast, worm, fly and mouse, and where do their annotated longevity directions agree?

Use the immutable Build 21 model CSV and the existing human-targeted website ortholog TSV. Join on species plus Entrez ID, never symbol alone. Retain every matching GenAge record and preserve its ID. Missing IDs cannot be matched. Exact mapping duplicates are removed with counts recorded. Human GenAge membership is contextual annotation, not a fifth species with experimental lifespan evidence.

The primary evidence set includes only Pro-Longevity and Anti-Longevity records. Necessary for fitness, Unclear and Unannotated remain in the audit and full joined table but do not establish directional longevity evidence. Each human anchor/species cell collects all primary directions: Pro, Anti, Mixed, or No directional evidence. At least two species with primary evidence defines a cross-species candidate. Agreement requires every represented primary-evidence cell to contain the same single direction. This is annotation agreement, not proof of a conserved causal effect or equivalent interventions.

Human anchors are overlapping sets, not independent orthogroups. Many-to-many mappings and paralogs are retained and flagged. No new orthology is inferred transitively. No mapping is not evidence of no homolog. Export database version, algorithm and confidence are unknown. File timestamps do not establish retrieval dates. Preserve the export as a frozen local snapshot with SHA-256 provenance.

Identifier audit rule: where one human Entrez ID has multiple exported symbols, quarantine all mappings for that anchor before analysis. The snapshot has conflicts at IDs 3020 (H3F3A/H3F3B) and 6606 (SMN1/SMN2). Coverage and anchor counts use only retained mappings; raw and quarantined counts are reported separately.

## Figure contract

Provisional conclusion: the frozen mapping identifies human-anchored candidates with directional longevity evidence across model species, subject to incomplete and unequal mapping coverage. Numerical claims are filled only after analysis.

Archetype: quantitative grid, with a principal evidence matrix and two compact supporting panels. Panel a shows species-specific mapping coverage among unique genes with primary evidence; panel b shows counts of human anchors by number of primary-evidence species; panel c shows all candidates with at least three species (or a transparently selected top 25 if there are more). The complete table always includes all anchors.

Double-column width 183 mm; height adjusted to displayed rows, maximum 230 mm. Sans-serif text 7 pt, editable SVG/PDF text, PNG preview and 600 dpi TIFF. Pro/Anti/Mixed are encoded by both colour and letters; missing evidence uses a dash. No statistical tests, error bars or independent replicate claims: this is a descriptive catalogue analysis. Source CSVs accompany every panel. Final submission requirements and biological validation remain pending.

## Verification

Check input hashes before/after; mapping uniqueness and join coverage; duplicate record keys; species and direction denominators; exact coverage and distribution totals; independent candidate-count reconstruction from raw inputs; SVG editable text and TIFF resolution; visual review of the rendered PNG. Investigate unmatched rows and symbol disagreements and retain their audit tables.
