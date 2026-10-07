# Full bibliography collection audit — 2026-09-18

**Historical stage report:** the identity problem described below has since been resolved for the exploratory figures using the official HAGR-R export and targeted checks. See [`verified_timeline_report.md`](verified_timeline_report.md) for current results. This document preserves why label-based mapping was rejected.

## Result

The public GenAge model-organism bibliography was collected completely against its advertised reference count: **710 of 710 references**, with **2,754 GenAge observation backlinks**. All references have parsed publication years (1988–2018). No duplicate nonempty PubMed URLs were found. These years describe this bibliography snapshot, not the full historical span of ageing genetics.

The source snapshots and per-page URL/retrieval metadata are in `../results/bibliography/raw/`. `../results/bibliography/source_manifest.json` records their SHA-256 hashes. The collector uses serial requests and one-second spacing, with limited backoff for connection failures. It resumes from cached pages.

## Checks and findings

- Title-sorted pagination returned all 710 distinct reference IDs. Author-sorted pagination returned only 707 because three IDs repeated across page boundaries; those pages were retained as an audit trail but excluded from the successful source manifest.
- All 25 observation/reference pairs in the validated pilot are present in the complete bibliography.
- Label matching tentatively dates 1,707 of 2,205 Build 21 entries. **This is not validated gene coverage.** There are 553 unresolved backlinks (500 distinct labels), including historical names, formatting differences, and duplicate entry names.
- Crucially, one label match contradicts the pilot: LibAge reference 4515 lists observation link `lid=2725` as `din-1` (candidate GenAge ID 177). Following the actual link opens **aak-2, GenAge ID 1**, matching the pilot. The destination HTML was saved in `../results/route_probe/lid_2725.html`.

Evidence URLs:

- https://libage.ageing-map.org/entries/entry/4515/
- https://genomics.senescence.info/genes/details.php?lid=2725

The incorrect label shows that even an apparently exact symbol/name match is insufficient evidence of gene identity. Do not repair only this example and assume all other labels are correct.

## Interpretation and next step

The bibliography collection is complete, but gene assignment is not. The tables in `../results/timeline/` are **candidate diagnostic outputs only** and must not support scientific conclusions. `audit.json` records the failed identity check; the plotting script refuses to generate figures while it remains unresolved.

Next, establish observation-to-entry identity from the actual GenAge destinations or from each Build 21 gene detail page and its listed observation IDs. Cross-check that mapping against the bibliography, report any unmatched references or entries, and only then regenerate earliest-reference and cumulative counts. Retain the original LibAge labels as source evidence, not as authoritative identifiers.

The annual bibliography-reference table can describe the captured bibliography, but cannot substitute for a validated gene-discovery timeline. No canonical or exploratory figure was produced in this stage.

## Reproduction

From the project root:

```text
python "task/task 3/scripts/collect_bibliography.py"
python "task/task 3/scripts/derive_timeline.py"
```

The first command uses the complete cache; the second is entirely local. The source Build 21 file remains unchanged. No Git commit is available because the project is not a Git repository.
