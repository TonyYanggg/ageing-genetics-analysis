# All-gene atlas completion audit — 26 September 2026

Objective: use the same main-figure style for every gene in the complete four-organism Task 4 comparison. The requirement is exhaustive visualization of the qualifying pair set, including the two-species cases; no six-family or three-species-only selection replaces it.

## Review figure contract

- **Core conclusion:** every GenAge gene participating in the 206 recovered direct OMA relationships among the four principal organisms can be located and traced without implying that connected rows are conserved mechanisms.
- **Archetype:** paginated quantitative grid / evidence atlas.
- **Target and size:** reviewer-facing supplementary atlas, 183 × 247 mm per page.
- **Backend:** Python only for construction, export and visual QA.
- **Evidence hierarchy:** gene cards show the evidence record; direct-pair strips show the exact source relationships; row labels and page/index navigation organize the material without adding biological edges.
- **Statistics:** none; this is a descriptive catalogue and contains no pooled effects or inferential tests.
- **Primary reviewer risks:** mistaking connected rows for orthogroups, treating GenAge Pro/Anti labels as standardized interventions, reading blank cells as biological absence, or assuming all papers were newly checked. These risks are stated in the header, footer and legend.

| Requirement | Current-state evidence | Result |
|---|---|---|
| Every qualifying gene | `gene_cards.csv` contains exactly the 320 species/Entrez keys independently reconstructed from all 206 frozen pairs; raw GenAge identities also reconciled | PASS |
| Every direct pair | `all_displayed_pairs.csv` preserves all 206 complete source rows exactly once; each edge label appears on its actual page | PASS |
| All records and paralogs retained | 321 GenAge IDs in 320 gene cards; all 133 networks independently reconstructed; multiple members in a species are separately stacked | PASS |
| Same main-figure style | Four fixed organism columns, outlined gene cards, row labels, relationship strips and bounded explanatory text throughout 36 pages | PASS; visual inspection |
| No invented experiment claims | 20 reviewed cards exactly match the existing evidence ledger; 300 catalogue cards have no invented intervention, primary PMID claim or experiment arrow | PASS |
| References are traceable | Every catalogue PMID reconciles to the frozen GenAge/LibAge reference mapping; original source URLs, hashes and timestamps are retained | PASS |
| Every page contains its gene-card content | All card text matched against actual SVG text and extracted combined/per-page PDF text | PASS, all 36 pages |
| Readable, unclipped layout | 3,558 text elements fit assigned boxes without overlap; PDF word bounds remain on-page; contact sheets cover all 36 pages | PASS |
| Export bundle | All 144 per-page PNG/SVG/PDF/TIFF files plus the 36-page bookmarked PDF checked; raster files fully decoded, physical dimensions 183 × 247 mm, TIFF 600 dpi | PASS |
| Navigation | 133 PDF bookmarks, 320 HTML gene entries, exact CSV gene/page index; search exercised for TOR1, GenAge 1939, a non-match and empty query | PASS |
| Existing work preserved | Complete-analysis, expanded-example and original-focused verifiers all pass; all input hashes in atlas summary reconcile | PASS |
| Documentation reflects current request | Root README, task index, task README, AGENTS, analysis plan, legend and changelog identify the complete atlas as the current visual output | Updated |

## Visual inspection

All four contact sheets were inspected, covering pages 1–36. Full-size PNG inspection included the original examples, a multiple-paralog page and the dual-name record on page 31. Python-rendered PDF inspection included pages 3, 15 and 36, including the densest five-member worm column. The page 31 dual-name label initially exceeded its box; explicit separate name lines fixed it while preserving both identifiers. No gene or relationship was removed to make the layout fit. Some pages have extra white space because networks are kept together rather than split between pages.

## Scope of verification

This audit establishes full visualization and source fidelity. It does not establish comparable causal effects, conserved mechanisms, experimental quality or new primary-paper validation for the 300 catalogue-only genes. The scientific qualifications for the reviewed experiments remain in the row notes and linked legends. Source gaps outside the recovered pair set remain in the all-nine audit.

## Reproducibility

`build_all_gene_atlas.py` performs a full offline build by default. During layout recovery, unchanged per-page exports were reused only after their complete SVG text sequence matched the newly measured layout; the combined PDF and all tables were regenerated. The final verifier then checked every existing raster/vector export and rendered card against current source data. The default full-build path does not depend on reuse.

Machine-readable evidence: `results/all_gene_atlas/verification.json`, `summary.json`, `text_bounds.json`, `gene_cards.csv`, `all_displayed_pairs.csv` and `page_index.csv`. Existing supporting verifiers were run successfully in the same work session. No Git repository or generating commit exists.

## Final review-readiness pass

The atlas was rebuilt cleanly from the frozen inputs after the initial completion audit. The post-build verifier again passed every source-to-card, source-to-edge, SVG/PDF text, raster decode, dimensions, DPI, bookmark and index check. All four contact sheets and representative full-size pages were reinspected, including the densest paralog layout and the dual-name worm record. No clipping, overlap, inconsistent type hierarchy or ambiguous colour-only encoding was found. Sparse pages retain white space intentionally because a connected candidate network is never split between pages. The atlas is ready to send for scientific review; Pedro has not yet approved its biological scope or use in the manuscript.
