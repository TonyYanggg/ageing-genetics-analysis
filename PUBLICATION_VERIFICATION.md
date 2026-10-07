# Publication verification — 7 October 2026

The public archive was prepared from the existing research workspace without modifying its original datasets, scripts or figures. Temporary Python bytecode and internal working instructions were omitted. Scientific scope requirements were retained as task-local analysis-scope notes. Public documentation uses repository-relative navigation and preserves evidence limitations.

## Checks performed on the publication copy

| Check | Outcome |
|---|---|
| Copy integrity | All 4,040 non-documentation research files initially matched their original SHA-256 hashes |
| Baseline inventory | Recalculated 2,205 raw records across nine organisms; matched saved organism counts |
| Reference timeline | Recalculated earliest reference dates for 2,203 entries and 708 distinct references; matched annual/cumulative series and pinned input hashes |
| Lifespan-extension verifier | PASS: 830 plotting records and four maximum labels reconstructed from raw data |
| Complete orthology verifier | PASS: 206 exact pairs, 320 genes, six species comparisons, all-nine coverage and exports |
| Focused examples verifier | PASS: ten records and ten source pairs, labels and exports |
| Expanded examples verifier | PASS: six families, 20 genes, 21 exact pairs, candidate screen and evidence ledger |
| Full atlas verifier | PASS: all 320 genes, 206 pairs, 133 networks, 36 pages and four export formats; 3,558 text elements checked |
| Human associations verifier | PASS: 307 genes retained; 306 queryable targets; 422,250 unique pairs matched to checksummed sources; evidence scores and output hashes reconciled |
| Public text review | No internal conversation references in project documentation or research scripts; source taxonomy identifiers were preserved unchanged |
| Credential-pattern scan | No matches for the tested common token/private-key patterns; this is a bounded scan, not a security certification |
| Large-file handling | Nine original files at least 10 MiB are stored through Git LFS; source bytes are preserved |

Verification scripts refreshed their task-local reports on the public copy. A regenerated report's machine-specific source path was made repository-relative. These documentation/report updates do not change biological results. The original local archive remains unchanged.

## Reproduction limits

Checks used the installed Python 3.13 environment and versions recorded in `requirements.txt`. A clean dependency installation and macOS/Linux execution were not tested. External collection scripts were not rerun; checks used the frozen responses already in the archive. Earlier live-source checks and visual figure reviews remain dated historical evidence in each task's documentation.

The initial public Git commit records publication of this archive, not the dates on which the original analyses were performed. No earlier Git history has been reconstructed or backdated.
