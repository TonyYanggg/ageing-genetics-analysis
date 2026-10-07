# Task 4 analysis scope and verification

## Current objective

The latest user decision (26 September 2026) is the same main-figure style for ALL genes in the four-organism qualifying pair set. The primary visual deliverable is `figures/all_gene_atlas/all_gene_ortholog_atlas.pdf`, with the same four columns and boxed cards across 36 pages. Include all 320 species/Entrez genes (321 GenAge records), all 206 exact source pairs, all 133 candidate networks and all two-species cases. Do not substitute the six selected families or only the 22 broad networks for this exhaustive visualization. `index.html` locates genes across pages. The six-family experiment figure and count overview remain supporting companions. Preserve all original figures. Pedro has approved only the original three-example format.

There is no human component in Task 4. Human gene–disease associations belong to Task 5. Functional enrichment (Task 6) is closed unless the user reopens it.

## Scientific safeguards

- GenAge Build 21 includes nine organisms and 2,205 model entries. The first four have 2,138 (97.0%), which justifies their review focus but not an absence claim for the other five. Retain the all-nine coverage and gap ledgers as supporting provenance.
- Use OMA All.May2026 explicit model-to-model protein-ortholog pairs. Preserve actual gene identifiers, GenAge IDs, relation cardinalities and release provenance. Never infer additional pairs or orthogroups by unconstrained transitive merging, shared symbols, pathway membership or a human anchor.
- Distinguish orthology, recurring GenAge lifespan association, and comparable causal effect. A GenAge Pro/Anti label is a compressed annotation, not a standardized experimental result.
- For every atlas gene retain GenAge IDs, original lifespan/influence fields, source pairs and catalogue-linked PMIDs. For the 20 reviewed experiment cards also retain the exact intervention, outcome, primary PMID, review depth and qualifications. The other 300 cards must be explicitly labelled GenAge annotation; never fabricate an intervention, add experimental arrows, or imply their papers have been reviewed. All papers in the reference index retain their original source provenance. Preserve the worm ATG7 annotation/selected-experiment distinction, yeast AMPK replicative context, yeast SOD2 sch9Δ requirement and connective-tissue-specific mouse Sod2 deficiency.
- The 133 connected candidate networks (22 spanning at least three species) are layout/search units. Never call them inferred OMA orthogroups or fill absent direct pairs. Every direct pair must be explicit in the atlas, with source-oriented cardinalities. Preserve all paralogs as separate cards within their organism column. Gene numbers Y/W/F/M reset within each network. Keep the six-family selection rationale as history, not the atlas inclusion criterion.
- Keep the selected mouse S6K1 sex specificity, mouse receptor paralogy and TOR yeast one-to-many relation visible. The 2011 S6K1 correction concerns a worm double-mutant strain; it confirms the lifespan epistasis result in a rebuilt strain and does not report a change to mouse survival. TOR and S6K are one pathway, not independent mechanism counts.
- A missing or unselected cell does not show lack of ortholog, GenAge evidence or lifespan function. OMA covers seven of the nine organisms; source gaps and noncoding/identifier limitations remain in the refined audit.
- Preserve historical human-anchor, OrthoDB and broad OMA outputs as labelled provenance. Do not overwrite raw files in shared `GenAge Project/` or `Other downloads from HAGR/`.

## Deliverables and verification

Run `scripts/build_all_gene_atlas.py` and `scripts/verify_all_gene_atlas.py` for the atlas and all card/pair/page/reference tables. Verification must cover all 320 gene keys, all 206 pair rows, all 133 networks, every page's actual rendered SVG/PDF card text and every PNG/SVG/PDF/TIFF export, not a selected sample alone. Inspect contact sheets for all pages plus full-size dense and representative PNG/PDF pages. Run the complete-analysis, expanded-example and original-focused verifiers to preserve supporting work. Builders and verifiers work offline from frozen inputs. The broader all-nine provenance remains unchanged.

Pedro approved the original format, not the complete atlas or expanded selection. Neither illustration nor catalogue proves evolutionary conservation of longevity mechanisms. Update these instructions, the task README, index, legend and change history if the figure or status materially changes.
