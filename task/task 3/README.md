# Task 3 — Timeline of ageing-gene discoveries

## Request

Create a review figure illustrating the timeline of discoveries represented in GenAge. Pedro explained that entries can have multiple references and that reference information is available on the GenAge website rather than in the downloadable dataset; he suggested scraping the site.

## Status

Exploratory figures ready for review. Authoritative observation IDs date **2,203 of 2,205 entries**, using **708 references** and **2,732 observation/reference links**. Ten label-based assignments were corrected and checked directly. Primary papers were externally recovered for the two undated entries, but remain in a separate sensitivity layer. See [`notes/verified_timeline_report.md`](notes/verified_timeline_report.md) for validation and [`notes/missing_reference_investigation.md`](notes/missing_reference_investigation.md) for the two-paper decision.

The current presentation drafts are `figures/review/timeline_review` and `figures/review/timeline_by_organism_review`, each in PNG, SVG and PDF. See [`notes/review_captions.md`](notes/review_captions.md) for captions, reproduction and revision checks. Original exploratory exports remain available. No Pedro-approved final figure exists yet.

## Inputs

| Input | Version | Role | Integrity |
|---|---|---|---|
| GenAge/LibAge pages and bibliography | Cached 2026-09-17–18 | Reference years, observation links, and targeted identity checks | Saved raw responses, dates, URLs and hashes |
| Official HAGR-R `models.csv` and `longevity.csv` | Commit `4e83b50142ea4b8c6625a7a86b6005ff3f606ed6` | Authoritative observation-to-entry IDs | Saved in `results/hagr_r_snapshot/` with manifest and hashes |
| Shared GenAge Build 21 downloads | Build 21 | Stable gene identifiers and baseline records | Raw and immutable |

## Figure contract

The working claim, evidence chain, export plan, and review risks are recorded in [`notes/figure_contract.md`](notes/figure_contract.md).

## Work plan

1. **Complete:** Inspect the GenAge entry and LibAge reference structures.
2. **Complete:** Build and validate a small extraction prototype for entries, observations, references, and publication years.
3. **Complete:** Collect the LibAge bibliography and its GenAge observation backlinks using serial, cached requests, following Pedro's suggestion in the project email.
4. **Complete:** Replace unreliable display labels with official observation-to-entry IDs, cross-check the live sample and targeted exceptions.
5. **Complete:** Build exploratory cumulative-entry, annual-publication and organism-specific views, with source-data tables.
6. **Complete:** Prepare separate review exports with coverage labels, captions and independently recalculated plotting counts.
7. **Next:** Review the main/supporting figure choice and trace major annual increases to contributing papers before adding historical interpretations.

The pilot extractor is `scripts/scrape_genage_references.py`. It defaults to 12 representative entries spanning all nine organisms, caches raw HTML, writes separate observation and reference tables, and refuses more than 25 IDs.

```text
python "task/task 3/scripts/scrape_genage_references.py"
```

## Pilot outputs

| Output | Status | Purpose |
|---|---|---|
| `results/pilot/raw_html/` | Supporting | Dated source snapshots used for parser validation |
| `results/pilot/raw_reference_html/` | Supporting | Cached LibAge responses for each sampled observation |
| `results/pilot/genage_observations_pilot.csv` | Supporting | One row per extracted GenAge observation |
| `results/pilot/genage_references_pilot.csv` | Supporting | One row per observation-reference association |
| `notes/pilot_extraction_report.md` | Supporting | Coverage, validation, and interpretation report |

The validated pilot contains 12 GenAge entries, 25 observations, 25 observation-reference associations, and 22 unique references. It spans all nine organisms in Build 21 and has no missing publication years. These are validation outputs, not the full timeline dataset.

## Source access note

The pilot was collected on 2026-09-17 from the [GenAge gene pages](https://genomics.senescence.info/genes/) and their linked [LibAge GenAge bibliography](https://libage.ageing-map.org/entries/sources/genage_model/). Pedro explicitly suggested scraping the references in the supplied project email. Full collection follows that scope with serial requests, caching, and one-second spacing, mindful of the server-load guidance in the [HAGR terms](https://genomics.senescence.info/legal.html).

## Full collection and audit

Run `python "task/task 3/scripts/collect_bibliography.py"` to collect the bibliography and its observation backlinks. Cached responses are reused on restart; request failures stop the run. `results/bibliography/collection_status.json` is written only after completion. Intermediate CSV files are partial until that marker exists.

The source's default author ordering omitted three references at page boundaries in the first collection attempt (707 distinct IDs versus 710 advertised). Title ordering is used instead, with a mandatory total-count check.

The original `derive_timeline.py` preserves the failed label-based audit in `results/timeline/`; those candidate counts are superseded. The current pipeline is:

```text
python "task/task 3/scripts/fetch_hagr_snapshot.py"
python "task/task 3/scripts/derive_verified_timeline.py"
python "task/task 3/scripts/plot_timeline.py"
```

The first command verifies or retrieves the pinned official snapshot; the others use local cached sources. `resolve_snapshot_gaps.py` reproduces targeted destination checks. The all-page collector `verify_genage_identities.py` provided a partial live sample before the bulk export was found; a full crawl is no longer required for this analysis.

| Output | Status | Provenance |
|---|---|---|
| `figures/review/timeline_review.*` and `timeline_by_organism_review.*` | Current review presentation | `plot_timeline.py --review`; same verified annual tables, plus entry coverage for denominators; hashes in review manifest |
| `results/bibliography/raw/` | Supporting source snapshots | Public bibliography pages and paper details, with per-file URL/date metadata |
| `results/bibliography/source_manifest.json` | Supporting | Cached source paths and SHA-256 hashes |
| `results/bibliography/references.csv` | Complete source-derived table | 710 LibAge references; generated by `collect_bibliography.py` |
| `results/bibliography/observation_backlinks.csv` | Complete source-derived table | 2,754 GenAge backlinks; their labels are not authoritative gene identities |
| `results/timeline/` | Superseded diagnostic candidates | Earlier `derive_timeline.py` label matching; not figure inputs |
| `results/verified_timeline/` | Current exploratory source data | `derive_verified_timeline.py`; official IDs joined to Build 21 and bibliography |
| `figures/timeline_exploratory.png` / `.svg` | Current draft | `plot_timeline.py`; verified annual entry and publication tables |
| `figures/timeline_by_organism_exploratory.png` / `.svg` | Current supporting draft | `plot_timeline.py`; verified annual entry table, all nine organisms |

The original label conflict is resolved using the pinned official export. All 2,205 entry identities match Build 21; 282 observations from 205 live-sampled entries agree with the export. All 25 pilot triples match, and all ten corrected assignments were checked directly. Counts are GenAge entries, not deduplicated biological genes.

Twenty-two dangling links were excluded after their destinations returned explicit entry-not-found pages. ZC395.10 and Mir29a have no database reference blocks and remain undated in the canonical analysis. External searching recovered high-confidence 2015 and 2016 primary papers; adding them changes only 2 of 2,205 entries (0.09%) and does not alter the interpretation. ASF1 has one observation with an empty reference response, so its date reflects available references only. The captured publication range is 1988–2018; this does not imply discoveries ended in 2018. See the reports for all limitations.

## Provenance requirements

Before producing the first full timeline, record the extraction date, requested URLs, scraper version, failures, duplicate-reference rules, validation sample, and the decision on earliest-reference versus all-reference display. The first current figure must have a row in an output-provenance table following [`../TASK_README_TEMPLATE.md`](../TASK_README_TEMPLATE.md).

No Git commit is available because this project is not currently a Git repository.

## History

See [`CHANGELOG.md`](CHANGELOG.md).
