# Task 3 change log

- **2026-09-18:** Added separate PNG/SVG/PDF review exports using `plot_timeline.py --review`, organism coverage labels and full captions. Independently recomputed all entry years, 279 organism/year counts and 31 annual reference counts from entry-level tables; all matched. Verified manifest hashes and inspected PNGs. Revised the figure claim to earliest associated reference dates. Raw data and original exploratory exports remain unchanged.

- **2026-09-18:** Investigated the two entries lacking GenAge/LibAge references. Recovered a direct 2015 `daf-41/ZC395.10` lifespan paper and the 2016 brain miR-29 paper whose table exactly matches the `Mir29a` GenAge summary. Kept both outside the canonical database-derived curves, recorded them as external candidates, and documented that adding them changes only 2/2,205 entries (0.09%). Added an explicit validation summary for the figures.

- **2026-09-18:** Found and pinned the lab's official HAGR-R export. All 2,205 entry identities match Build 21, and 282 sampled live observations from 205 entries agree with its observation mapping. Directly checked ten corrected label assignments and five old hidden flags. Excluded 22 confirmed dangling links; documented two undated entries and one observation with no reference. Generated first PNG/SVG timeline drafts from 2,203 dated entries and 708 references, with source tables and a verified timeline report. Earlier label-based outputs remain as superseded diagnostics.

- **2026-09-18:** Completed collection of 710 references and 2,754 GenAge backlinks. All 25 pilot observation/reference pairs were recovered, but backlink 2725 was labelled din-1 despite opening aak-2 (GenAge ID 1). Label-based mapping also left 553 unresolved backlinks. Saved the audit and candidate diagnostic tables; blocked plotting pending authoritative identity mapping.

- **2026-09-18:** Started serial, cached collection through the LibAge bibliography and GenAge observation backlinks, following Pedro's scraping suggestion. Added a local Build 21 mapping and timeline audit. Author-sorted pagination returned 707/710 distinct references; switched to title ordering and retained the failed attempt's cache for audit.

- **2026-09-17:** Defined the timeline figure contract, evidence chain, export plan, and review risks.
- **2026-09-17:** Confirmed that GenAge detail pages expose observation-level primary references with LibAge and PubMed identifiers.
- **2026-09-17:** Added a cached, rate-limited pilot extractor capped at 25 GenAge entries. Full-catalogue extraction has not yet been attempted.
- **2026-09-17:** Validated the pilot on 12 entries spanning all nine Build 21 organisms: 25 observations, 25 observation-reference links, 22 unique references, and no missing publication years. The `age-1` references reproduced the expected 1988, 1990, 1995, 2005, and 2008 dates.
