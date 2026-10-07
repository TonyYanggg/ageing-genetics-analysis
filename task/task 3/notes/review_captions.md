# Task 3 review captions and checks

Status: review drafts, not approved final figures. Prepared 2026-09-18.

## Overall timeline

**GenAge entries by earliest associated reference year.** **a,** Cumulative number of GenAge Build 21 model-organism entries assigned to the earliest dated reference associated with each entry in the captured GenAge/LibAge sources. Each GenAge ID contributes once, including when several references share its earliest year. **b,** Annual numbers of distinct LibAge reference IDs linked to mapped entries; a reference linked to multiple entries is counted once. The analysis dates 2,203 of 2,205 entries (99.91%) and includes 708 distinct references spanning 1988–2018. These are descriptive database counts; no statistical tests or uncertainty intervals are applied.

Dates represent the earliest available associated publication, not independently verified historical discovery dates or dates of database entry. The captured range does not imply that discoveries ended in 2018. Twenty-two backlinks with entry-not-found destinations were excluded. ZC395.10 and Mir29a remain in the coverage denominator but have no dates in the database-derived curves. External candidate papers are documented separately and have not been added to these plots. One ASF1 observation lacks a reference; that entry's date uses its other available reference evidence.

Source data: `../results/verified_timeline/annual_entry_timeline.csv` and `annual_reference_timeline.csv`. Entry-level evidence: `entry_references.csv` and `earliest_reference_provenance.csv` in the same directory.

## Organism timelines

**Earliest associated reference timelines across nine GenAge organisms.** Each panel shows cumulative GenAge entries using the same date definition as the overall timeline. Labels give dated entries divided by all Build 21 entries for that organism. All nine organisms are retained, including those represented by a single entry. Panels share the publication-year range but use separate y-axis scales; their apparent heights and slopes therefore cannot be compared directly as absolute counts or rates. The same missing-reference and historical-interpretation limitations apply.

Source data: `../results/verified_timeline/annual_entry_timeline.csv`; denominators from `entry_coverage.csv`.

## Review exports and reproduction

The main figure is `../figures/review/timeline_review`; the supporting organism figure is `../figures/review/timeline_by_organism_review`. Both have PNG (300 dpi), SVG with editable text, and PDF exports. Original exploratory exports are retained. Final journal dimensions and submission formats remain to be chosen after figure review.

From the project root: `python "task/task 3/scripts/plot_timeline.py" --review`.

`../figures/review/timeline_provenance.json` records generating-script, input and output hashes. The review mode changes presentation only and uses the existing database-derived source tables.

## Checks completed for this revision

- Independently recomputed earliest years from entry/reference rows and compared every dated entry with the coverage table.
- Recomputed all 279 organism/year annual and cumulative counts directly from entry-level years; all matched the plotting inputs.
- Recomputed all 31 annual distinct-reference counts from entry/reference rows; all matched.
- Verified all input/output hashes listed in the review manifest; checked SVG text elements and PDF file signatures.
- Visually inspected both PNG exports: no clipped labels or overlaps observed; coverage labels and separate y-axis-scale warning are readable.

These checks verify the saved data-to-figure chain, not the correctness of every underlying paper or historical discovery claim. PDF signatures were checked, but PDF rendering and text selection were not separately tested. Historical milestone annotations are not included; any such additions require their own primary-source validation.

## Next scientific decision

Review the overall plot as the main figure and the nine-organism plot as supporting material. Before adding biological explanations for the large steps, trace the largest annual increases to their contributing papers. Their size alone does not establish a change in discovery rate or a particular experimental cause.
