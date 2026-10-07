# Task 3 — verified mapping and exploratory timeline

## Result

The first exploratory timeline dates **2,203 of 2,205 Build 21 entries (99.91%)**, using **2,732 observation/reference links** and **708 distinct LibAge references**, spanning 1988–2018. All nine organisms are included.

This is an entry timeline based on the earliest **dated associated reference**, not a reconstruction of definitive historical gene-discovery dates. The end of the captured reference range must not be interpreted as the end of discoveries in ageing genetics.

## Authoritative mapping

The lab's public [HAGR-R GenAgeModel export](https://github.com/maglab/HAGR-R/tree/4e83b50142ea4b8c6625a7a86b6005ff3f606ed6/data-raw/GenAgeModel) provides explicit observation-to-gene IDs in `longevity.csv`. The snapshot is pinned to commit `4e83b50142ea4b8c6625a7a86b6005ff3f606ed6`; downloaded bytes, URLs, dates and hashes are retained in `../results/hagr_r_snapshot/`.

- All 2,205 exported entry IDs, symbols and organisms match the local Build 21 catalogue.
- A live-page sample of 282 observations from 205 entries agrees with the exported observation-to-gene IDs. This combines the sequential live sample with the earlier pilot, deduplicated by observation ID.
- All 25 pilot gene/observation/reference triples are recovered exactly.
- The earlier label-based approach assigned 10 links to incorrect entry IDs. Each corrected destination was checked directly against the current GenAge page. The corrections are listed in `../results/verified_timeline/corrected_label_assignments.csv`.
- Five observations marked hidden in the export are displayed on the live pages; direct checks confirmed their IDs and they are included.

The full live-page collector stopped after a server error at entry 202. Its saved pages remain useful verification evidence. The official bulk export replaced the need to finish that full crawl; this report does not claim every gene page was individually checked.

## Coverage and exclusions

| Organism | Dated entries | Build 21 entries |
|---|---:|---:|
| Caenorhabditis briggsae | 1 | 1 |
| Caenorhabditis elegans | 888 | 889 |
| Danio rerio | 1 | 1 |
| Drosophila melanogaster | 202 | 202 |
| Mesocricetus auratus | 1 | 1 |
| Mus musculus | 135 | 136 |
| Podospora anserina | 3 | 3 |
| Saccharomyces cerevisiae | 911 | 911 |
| Schizosaccharomyces pombe | 61 | 61 |

- **22 dangling bibliography backlinks:** 21 observation IDs are absent from the snapshot and one has no assigned gene. All 22 actual destinations returned an explicit “Entry not found” page. They are retained in `excluded_dangling_backlinks.csv`, with source snapshots, and excluded from mapped entry counts. Two bibliography references have no remaining mapped links, reducing the publication count from 710 to 708.
- **Two undated entries in the database-derived analysis:** ZC395.10 (ID 2081) and Mir29a (ID 2173). Neither has an observation/reference block in the current GenAge page or an observation in the export. Primary papers were recovered externally for both entries, but they are kept out of the canonical curves so that database-derived and manually curated evidence are not mixed. See [`missing_reference_investigation.md`](missing_reference_investigation.md).
- **One observation without a reference:** ASF1 observation 2158 is displayed, but its LibAge reference endpoint returns an empty list and the export has no bibliography ID. ASF1 has another dated reference, so the entry is included using available evidence. Its true earliest associated date cannot be established from an empty reference.

Do not infer that every removed link is biologically incorrect: the exclusion is about inability to associate it with a current entry using these sources.

## Figures and interpretation

`../figures/timeline_exploratory.png` and `.svg`:

- **a:** cumulative entries by earliest dated associated publication year; each GenAge entry is counted once.
- **b:** annual distinct LibAge references attached to mapped entries; a paper linked to many genes is counted once.

`../figures/timeline_by_organism_exploratory.png` and `.svg`:

- The same entry definition, separately for all nine organisms. The y-axis scale differs between panels and is labelled accordingly.

The largest annual increases in first-referenced entries are 2010 (462), 2007 (390), and 2008 (324). These are properties of the captured database and associated literature, not unbiased estimates of scientific discovery activity.

GenAge IDs define the counting unit. Multiple entries for the same biological gene are not silently merged. No inferential tests, uncertainty intervals or smoothing are used; these are descriptive database counts.

## Validation performed

The figure inputs and calculations were checked at four levels:

1. **Entry identity:** all 2,205 HAGR-R IDs, symbols and organisms match Build 21.
2. **Mapping:** 282 live observations from 205 entries agree with the pinned export; all 25 pilot gene/observation/reference triples were recovered; ten misleading label assignments were checked directly and corrected.
3. **Calculation:** cumulative series are monotonic, their final total is 2,203, annual publication counts sum to 708, mapped observation/reference pairs are unique, and every included reference has a year. Tied earliest references remain in the provenance table.
4. **Reproducibility and display:** source and output SHA-256 hashes are recorded, figures were regenerated from the saved CSVs, and PNG/SVG exports were visually inspected for labels, clipping, panel order and axis scales.

This is strong validation of extraction, mapping and plotting consistency. It is not an independent manual audit of every one of the 2,203 entry dates, and it cannot establish the true historical discovery year beyond the operational definition of earliest associated GenAge reference.

## Provenance and reproduction

The current source-data tables are in `../results/verified_timeline/`, especially `entry_references.csv`, `earliest_reference_provenance.csv`, `entry_coverage.csv`, and the two annual timeline CSVs. `audit.json` includes input hashes and validation results. Earlier `../results/timeline/` files are superseded diagnostic candidates and are not figure inputs.

From the project root, with the saved source snapshots present:

```text
python "task/task 3/scripts/fetch_hagr_snapshot.py"
python "task/task 3/scripts/derive_verified_timeline.py"
python "task/task 3/scripts/plot_timeline.py"
```

The first command verifies or retrieves the pinned export. The other two are local. `collect_bibliography.py` and `resolve_snapshot_gaps.py` reproduce the cached web collections when needed. Shared downloaded source files were not modified.

The figures are exploratory drafts for review, not Pedro-approved final figures. PNG and editable-text SVG were visually inspected for clipping, labels, panel ordering and axis scales. Publication exports and historical milestone annotation remain a later refinement.
