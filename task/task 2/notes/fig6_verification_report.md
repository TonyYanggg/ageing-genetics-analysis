# Figure 6 independent verification report

Audit date: 2026-08-20  
Verification script: `scripts/verify_fig6_data.py`

## Raw source verified

- Exact filename: `genage_models.csv`
- Full local path: `GenAge Project\models_genes\genage_models.csv`
- Dimensions: **2,205 rows x 8 columns**
- SHA-256: `8717b30212a2f6039a7b70c28dff8bc5696c166e04a150f95302f90b8585c2a0`
- Exact columns, in file order: `GenAge ID`, `symbol`, `name`, `organism`, `entrez gene id`, `avg lifespan change (max obsv)`, `lifespan effect`, `longevity influence`
- Relevant raw fields: organism = `organism`; symbol/name = `symbol` and `name`; Entrez ID = `entrez gene id`; longevity influence = `longevity influence`; lifespan effect = `lifespan effect`; magnitude = `avg lifespan change (max obsv)`.

The scripts used for Figures 1, 5 and 6 all point to the same local raw path, and the independently recomputed organism/effect counts exactly match the saved Figures 1 and 5 tables. The saved tables do not embed a source checksum, so the raw SHA-256 above is the definitive byte-level identity for this audit.

### Raw lifespan-effect values before cleaning

| raw value (repr) | count | characters | leading/trailing whitespace | empty/missing in source |
| --- | --- | --- | --- | --- |
| 'Increase and Decrease' | 144 | 21 | False | False |
| 'Increase' | 1155 | 8 | False | False |
| 'Decrease' | 892 | 8 | False | False |
| '' | 14 | 0 | False | True |

No capitalization or whitespace normalization is needed: all nonempty categories already exactly match `Increase`, `Decrease`, or `Increase and Decrease`. The 14 empty strings are explicitly recoded to `Missing/Other` only for tabulation; no record is dropped.

### Organism and raw effect counts

| organism | total | Increase | Decrease | Increase and Decrease | Missing/Other |
| --- | --- | --- | --- | --- | --- |
| Caenorhabditis elegans | 889 | 608 | 246 | 31 | 4 |
| Mus musculus | 136 | 56 | 68 | 10 | 2 |
| Saccharomyces cerevisiae | 911 | 298 | 516 | 90 | 7 |
| Drosophila melanogaster | 202 | 132 | 56 | 13 | 1 |
| Mesocricetus auratus | 1 | 1 | 0 | 0 | 0 |
| Podospora anserina | 3 | 3 | 0 | 0 | 0 |
| Schizosaccharomyces pombe | 61 | 56 | 5 | 0 | 0 |
| Danio rerio | 1 | 0 | 1 | 0 | 0 |
| Caenorhabditis briggsae | 1 | 1 | 0 | 0 | 0 |

Every organism total and effect count is an exact match to the prior Figures 1 and 5 tables. In particular, *S. pombe* has **61** raw GenAge records.

## Filtering verified

The Figure 6 dataset was reconstructed directly from the raw CSV with three simultaneous rules: organism is one of the four plotted organisms; the raw lifespan-effect value is exactly `Increase`; and `avg lifespan change (max obsv)` parses to a finite numeric value. A value is “usable” only if its stripped raw field is nonblank and converts to a finite number. No deduplication, aggregation, imputation, capping or winsorisation was applied.

Raw magnitude-field inspection:

- blank fields: **933**
- whitespace-only but nonempty fields: **0**
- literal NaN-like strings (`NaN`, `NA`, `N/A`, `null`, `none`, case-insensitive): **0**
- nonblank nonnumeric parse failures: **0**
- infinite values: **0**
- usable finite numeric values: **1,272**
- positive values: **1,272**
- zero values: **0**
- negative values: **0**

There are no unusual nonblank strings to list. All missing/unusable cases are blank raw fields.

### Magnitude-data availability by organism and effect category

| organism | organism_total_records | lifespan_effect_category | category_records | usable_numeric_magnitude_records | missing_magnitude_records | usable_magnitude_percent |
| --- | --- | --- | --- | --- | --- | --- |
| Caenorhabditis elegans | 889 | Increase | 608 | 453 | 155 | 74.506579 |
| Caenorhabditis elegans | 889 | Decrease | 246 | 163 | 83 | 66.260163 |
| Caenorhabditis elegans | 889 | Increase and Decrease | 31 | 0 | 31 | 0 |
| Caenorhabditis elegans | 889 | Missing/Other | 4 | 0 | 4 | 0 |
| Mus musculus | 136 | Increase | 56 | 46 | 10 | 82.142857 |
| Mus musculus | 136 | Decrease | 68 | 25 | 43 | 36.764706 |
| Mus musculus | 136 | Increase and Decrease | 10 | 0 | 10 | 0 |
| Mus musculus | 136 | Missing/Other | 2 | 0 | 2 | 0 |
| Saccharomyces cerevisiae | 911 | Increase | 298 | 217 | 81 | 72.818792 |
| Saccharomyces cerevisiae | 911 | Decrease | 516 | 218 | 298 | 42.248062 |
| Saccharomyces cerevisiae | 911 | Increase and Decrease | 90 | 0 | 90 | 0 |
| Saccharomyces cerevisiae | 911 | Missing/Other | 7 | 0 | 7 | 0 |
| Drosophila melanogaster | 202 | Increase | 132 | 114 | 18 | 86.363636 |
| Drosophila melanogaster | 202 | Decrease | 56 | 33 | 23 | 58.928571 |
| Drosophila melanogaster | 202 | Increase and Decrease | 13 | 0 | 13 | 0 |
| Drosophila melanogaster | 202 | Missing/Other | 1 | 0 | 1 | 0 |
| Mesocricetus auratus | 1 | Increase | 1 | 0 | 1 | 0 |
| Mesocricetus auratus | 1 | Decrease | 0 | 0 | 0 | NA |
| Mesocricetus auratus | 1 | Increase and Decrease | 0 | 0 | 0 | NA |
| Mesocricetus auratus | 1 | Missing/Other | 0 | 0 | 0 | NA |
| Podospora anserina | 3 | Increase | 3 | 1 | 2 | 33.333333 |
| Podospora anserina | 3 | Decrease | 0 | 0 | 0 | NA |
| Podospora anserina | 3 | Increase and Decrease | 0 | 0 | 0 | NA |
| Podospora anserina | 3 | Missing/Other | 0 | 0 | 0 | NA |
| Schizosaccharomyces pombe | 61 | Increase | 56 | 0 | 56 | 0 |
| Schizosaccharomyces pombe | 61 | Decrease | 5 | 0 | 5 | 0 |
| Schizosaccharomyces pombe | 61 | Increase and Decrease | 0 | 0 | 0 | NA |
| Schizosaccharomyces pombe | 61 | Missing/Other | 0 | 0 | 0 | NA |
| Danio rerio | 1 | Increase | 0 | 0 | 0 | NA |
| Danio rerio | 1 | Decrease | 1 | 1 | 0 | 100 |
| Danio rerio | 1 | Increase and Decrease | 0 | 0 | 0 | NA |
| Danio rerio | 1 | Missing/Other | 0 | 0 | 0 | NA |
| Caenorhabditis briggsae | 1 | Increase | 1 | 1 | 0 | 100 |
| Caenorhabditis briggsae | 1 | Decrease | 0 | 0 | 0 | NA |
| Caenorhabditis briggsae | 1 | Increase and Decrease | 0 | 0 | 0 | NA |
| Caenorhabditis briggsae | 1 | Missing/Other | 0 | 0 | 0 | NA |

The exact 830 record-level rows used by the main Figure 6 were written to `task\task 2\results\tables\fig6_verified_plotting_data.csv`. Each output row retains the zero-based raw row index and the one-based CSV line number (including the header as line 1).

## Sample sizes

| Organism | Current Figure 6 n | Independently recomputed n | Match |
| --- | --- | --- | --- |
| Saccharomyces cerevisiae | 217 | 217 | True |
| Caenorhabditis elegans | 453 | 453 | True |
| Drosophila melanogaster | 114 | 114 | True |
| Mus musculus | 46 | 46 | True |

The independent total is **830 = 217 + 453 + 114 + 46**. All four comparisons match exactly.

## Summary statistics

The standard deviation is the sample standard deviation (`ddof=1`); quartiles use pandas' default linear interpolation. These values were calculated from the raw-derived Figure 6 records, not from the existing summary file.

| organism | n | mean | median | sd | q1 | q3 | iqr | minimum | maximum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Saccharomyces cerevisiae | 217 | 35.788018 | 29 | 45.344476 | 20 | 40 | 20 | 10 | 600 |
| Caenorhabditis elegans | 453 | 33.913466 | 25 | 51.605627 | 16 | 40 | 24 | 4 | 1,000 |
| Drosophila melanogaster | 114 | 35.514825 | 30 | 24.761633 | 20 | 44.75 | 24.75 | 5 | 150 |
| Mus musculus | 46 | 22.841304 | 20 | 9.6615981 | 16 | 30 | 14 | 7 | 46 |

### Column-by-column comparison with `fig6_lifespan_magnitude_summary.csv`

| file | field | key | independent | existing | classification |
| --- | --- | --- | --- | --- | --- |
| fig6_lifespan_magnitude_summary.csv | n | Saccharomyces cerevisiae / Increase | 217 | 217 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Saccharomyces cerevisiae / Increase | 35.788018 | 35.788 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Saccharomyces cerevisiae / Increase | 29 | 29 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Saccharomyces cerevisiae / Increase | 45.344476 | 45.3445 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Saccharomyces cerevisiae / Increase | 20 | 20 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Saccharomyces cerevisiae / Increase | 40 | 40 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Saccharomyces cerevisiae / Increase | 20 | 20 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Saccharomyces cerevisiae / Increase | 10 | 10 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Saccharomyces cerevisiae / Increase | 600 | 600 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Saccharomyces cerevisiae / Decrease | 218 | 218 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Saccharomyces cerevisiae / Decrease | 40.018349 | 40.0183 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Saccharomyces cerevisiae / Decrease | 40 | 40 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Saccharomyces cerevisiae / Decrease | 16.194089 | 16.1941 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Saccharomyces cerevisiae / Decrease | 25 | 25 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Saccharomyces cerevisiae / Decrease | 50 | 50 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Saccharomyces cerevisiae / Decrease | 25 | 25 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Saccharomyces cerevisiae / Decrease | 10 | 10 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Saccharomyces cerevisiae / Decrease | 85 | 85 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Saccharomyces cerevisiae / Increase and Decrease | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Saccharomyces cerevisiae / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Caenorhabditis elegans / Increase | 453 | 453 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Caenorhabditis elegans / Increase | 33.913466 | 33.9135 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Caenorhabditis elegans / Increase | 25 | 25 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Caenorhabditis elegans / Increase | 51.605627 | 51.6056 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Caenorhabditis elegans / Increase | 16 | 16 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Caenorhabditis elegans / Increase | 40 | 40 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Caenorhabditis elegans / Increase | 24 | 24 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Caenorhabditis elegans / Increase | 4 | 4 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Caenorhabditis elegans / Increase | 1,000 | 1000 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Caenorhabditis elegans / Decrease | 163 | 163 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Caenorhabditis elegans / Decrease | 25.177914 | 25.1779 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Caenorhabditis elegans / Decrease | 22.200001 | 22.2 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Caenorhabditis elegans / Decrease | 14.244978 | 14.245 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Caenorhabditis elegans / Decrease | 13.5 | 13.5 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Caenorhabditis elegans / Decrease | 33.5 | 33.5 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Caenorhabditis elegans / Decrease | 20 | 20 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Caenorhabditis elegans / Decrease | 7 | 7 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Caenorhabditis elegans / Decrease | 84 | 84 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Caenorhabditis elegans / Increase and Decrease | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Caenorhabditis elegans / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Drosophila melanogaster / Increase | 114 | 114 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Drosophila melanogaster / Increase | 35.514825 | 35.5148 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Drosophila melanogaster / Increase | 30 | 30 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Drosophila melanogaster / Increase | 24.761633 | 24.7616 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Drosophila melanogaster / Increase | 20 | 20 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Drosophila melanogaster / Increase | 44.75 | 44.75 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Drosophila melanogaster / Increase | 24.75 | 24.75 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Drosophila melanogaster / Increase | 5 | 5 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Drosophila melanogaster / Increase | 150 | 150 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Drosophila melanogaster / Decrease | 33 | 33 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Drosophila melanogaster / Decrease | 32.690909 | 32.6909 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Drosophila melanogaster / Decrease | 30 | 30 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Drosophila melanogaster / Decrease | 17.566114 | 17.5661 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Drosophila melanogaster / Decrease | 20 | 20 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Drosophila melanogaster / Decrease | 42.599998 | 42.6 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | iqr | Drosophila melanogaster / Decrease | 22.599998 | 22.6 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | minimum | Drosophila melanogaster / Decrease | 4 | 4 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Drosophila melanogaster / Decrease | 75 | 75 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Drosophila melanogaster / Increase and Decrease | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Drosophila melanogaster / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Mus musculus / Increase | 46 | 46 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Mus musculus / Increase | 22.841304 | 22.8413 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Mus musculus / Increase | 20 | 20 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Mus musculus / Increase | 9.6615981 | 9.6616 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Mus musculus / Increase | 16 | 16 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Mus musculus / Increase | 30 | 30 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Mus musculus / Increase | 14 | 14 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Mus musculus / Increase | 7 | 7 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Mus musculus / Increase | 46 | 46 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Mus musculus / Decrease | 25 | 25 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Mus musculus / Decrease | 33.84 | 33.84 | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Mus musculus / Decrease | 33 | 33 | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Mus musculus / Decrease | 19.20434 | 19.2043 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q1 | Mus musculus / Decrease | 18 | 18 | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Mus musculus / Decrease | 43 | 43 | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Mus musculus / Decrease | 25 | 25 | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Mus musculus / Decrease | 6 | 6 | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Mus musculus / Decrease | 83 | 83 | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Mus musculus / Increase and Decrease | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Mus musculus / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Schizosaccharomyces pombe / Increase | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Schizosaccharomyces pombe / Increase | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Schizosaccharomyces pombe / Decrease | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Schizosaccharomyces pombe / Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | n | Schizosaccharomyces pombe / Increase and Decrease | 0 | 0 | exact match |
| fig6_lifespan_magnitude_summary.csv | mean | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | median | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | sd | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q1 | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | q3 | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | iqr | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | minimum | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |
| fig6_lifespan_magnitude_summary.csv | maximum | Schizosaccharomyces pombe / Increase and Decrease | NA | NA | exact match |

`exact match` means numerical identity (allowing equivalent integer/decimal representations). `rounding difference` means the existing value is the independent value formatted to six significant digits. No comparison is classified as a discrepancy.

## Maximum records

The maximum was found independently by sorting each organism's raw-derived Figure 6 records by numeric magnitude descending, with raw row index as a deterministic tie-breaker. Every maximum is unique; the four proposed labels are therefore correct. None of the four maximum-gene symbols is repeated within its organism in either the full raw dataset or the Figure 6 subset.

| organism | symbol | entrez gene id | parsed_numeric_lifespan_change | raw_row_index_zero_based | csv_line_number | raw_lifespan_change | lifespan effect | maximum_tie_count | symbol_records_in_raw_organism | symbol_records_in_figure6 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Saccharomyces cerevisiae | SIR2 | 851520 | 600.0 | 604 | 606 | 600 | Increase | 1 | 1 | 1 |
| Caenorhabditis elegans | age-1 | 174762 | 1000.0 | 6 | 8 | 1000 | Increase | 1 | 1 | 1 |
| Drosophila melanogaster | Vinc | 31201 | 150.0 | 2057 | 2059 | 150 | Increase | 1 | 1 | 1 |
| Mus musculus | Ghrh | 14601 | 46.0 | 1822 | 1824 | 46 | Increase | 1 | 1 | 1 |

Every tied maximum would be retained in `task\task 2\results\tables\fig6_verified_maximum_records.csv`; the observed tie count is one for every organism.

### Top 10 Figure 6 records per organism


### Saccharomyces cerevisiae

| rank | symbol | entrez gene id | parsed_numeric_lifespan_change | raw_row_index_zero_based | csv_line_number | raw_lifespan_change | lifespan effect | GenAge ID |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | SIR2 | 851520 | 600 | 604 | 606 | 600 | Increase | 606 |
| 2 | SCH9 | 856612 | 300 | 588 | 590 | 300 | Increase | 590 |
| 3 | GRE3 | 856504 | 100 | 1885 | 1887 | 100 | Increase | 1911 |
| 4 | MET17 | 851010 | 100 | 1939 | 1941 | 100 | Increase | 1965 |
| 5 | CYR1 | 853452 | 90 | 148 | 150 | 90 | Increase | 149 |
| 6 | RPL22A | 850750 | 90 | 1698 | 1700 | 90 | Increase | 1724 |
| 7 | UGA3 | 851384 | 75 | 1709 | 1711 | 75 | Increase | 1735 |
| 8 | UBR2 | 850711 | 67 | 1995 | 1997 | 67 | Increase | 2021 |
| 9 | FZF1 | 852638 | 66 | 1708 | 1710 | 66 | Increase | 1734 |
| 10 | SCP1 | 854549 | 65 | 590 | 592 | 65 | Increase | 592 |

### Caenorhabditis elegans

| rank | symbol | entrez gene id | parsed_numeric_lifespan_change | raw_row_index_zero_based | csv_line_number | raw_lifespan_change | lifespan effect | GenAge ID |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | age-1 | 174762 | 1,000 | 6 | 8 | 1000 | Increase | 7 |
| 2 | daf-2 | 175410 | 200 | 162 | 164 | 200 | Increase | 163 |
| 3 | let-363 | 172167 | 150 | 31 | 33 | 150 | Increase | 32 |
| 4 | osm-5 | 180585 | 150 | 484 | 486 | 150 | Increase | 486 |
| 5 | unc-13 | 172497 | 150 | 700 | 702 | 150 | Increase | 702 |
| 6 | unc-31 | 178233 | 150 | 702 | 704 | 150 | Increase | 704 |
| 7 | rps-23 | 178188 | 126 | 567 | 569 | 126 | Increase | 569 |
| 8 | daf-4 | 175781 | 120 | 165 | 167 | 120 | Increase | 166 |
| 9 | ucr-1 | 175772 | 110 | 695 | 697 | 110 | Increase | 697 |
| 10 | che-3 | 172593 | 100 | 120 | 122 | 100 | Increase | 121 |

### Drosophila melanogaster

| rank | symbol | entrez gene id | parsed_numeric_lifespan_change | raw_row_index_zero_based | csv_line_number | raw_lifespan_change | lifespan effect | GenAge ID |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Vinc | 31201 | 150 | 2057 | 2059 | 150 | Increase | 2083 |
| 2 | Nmdmc | 47895 | 120 | 2072 | 2074 | 120 | Increase | 2098 |
| 3 | Indy | 40049 | 92 | 364 | 366 | 92 | Increase | 366 |
| 4 | trh | 38065 | 90 | 2062 | 2064 | 90 | Increase | 2088 |
| 5 | 5-HT2A | 40575 | 90 | 2063 | 2065 | 90 | Increase | 2089 |
| 6 | InR | 42549 | 85 | 366 | 368 | 85 | Increase | 368 |
| 7 | snz | 31704 | 85 | 855 | 857 | 85 | Increase | 879 |
| 8 | fabp | 3772232 | 81 | 1936 | 1938 | 81 | Increase | 1962 |
| 9 | Gadd45 | 35646 | 77 | 826 | 828 | 77 | Increase | 849 |
| 10 | E(z) | 39203 | 76 | 821 | 823 | 76 | Increase | 844 |

### Mus musculus

| rank | symbol | entrez gene id | parsed_numeric_lifespan_change | raw_row_index_zero_based | csv_line_number | raw_lifespan_change | lifespan effect | GenAge ID |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Ghrh | 14601 | 46 | 1822 | 1824 | 46 | Increase | 1848 |
| 2 | Gdf15 | 23886 | 43 | 1845 | 1847 | 43 | Increase | 1871 |
| 3 | Ghr | 14600 | 40 | 301 | 303 | 40 | Increase | 302 |
| 4 | Igf1r | 16001 | 40 | 363 | 365 | 40 | Increase | 365 |
| 5 | Pou1f1 | 18736 | 40 | 504 | 506 | 40 | Increase | 506 |
| 6 | Eps8 | 13860 | 37 | 794 | 796 | 37 | Increase | 817 |
| 7 | Fgf21 | 56636 | 36 | 1774 | 1776 | 36 | Increase | 1800 |
| 8 | Txn1 | 22166 | 35 | 689 | 691 | 35 | Increase | 691 |
| 9 | Insr | 16337 | 33.299999 | 371 | 373 | 33.29999924 | Increase | 373 |
| 10 | Adcy5 | 224129 | 30 | 5 | 7 | 30 | Increase | 6 |

## S. pombe

| Quantity | Independently recomputed |
| --- | ---: |
| Total records | 61 |
| Increase | 56 |
| Decrease | 5 |
| Increase and Decrease | 0 |
| Missing/Other lifespan effect | 0 |
| Records with a usable numerical magnitude | 0 |
| Records without a usable magnitude | 61 |

Thus, all **61 / 61 (100%)** *S. pombe* records lack a usable magnitude. Its exclusion is justified solely by absent magnitude data, not by an absence of GenAge records.

## Plot integrity

The current `plot_final_extension_only` function was tested only after the independent dataset had been created. Its save routine was replaced in memory with a no-op, so no figure was written. The scatter x arrays and boxplot input arrays were captured and compared, in order, with the independently filtered record-level magnitudes.

| organism | verified rows | scatter x observations | boxplot input observations | x values identical | jitter changes x |
| --- | --- | --- | --- | --- | --- |
| Saccharomyces cerevisiae | 217 | 217 | 217 | True | False |
| Caenorhabditis elegans | 453 | 453 | 453 | True | False |
| Drosophila melanogaster | 114 | 114 | 114 | True | False |
| Mus musculus | 46 | 46 | 46 | True | False |

- captured scatter observations: **830**
- captured boxplot-input observations: **830**
- verified CSV rows: **830**
- unique raw row indices: **830**
- deduplication/aggregation operations in the current main plotting function: **none detected**
- jitter: changes only the y display coordinate; every captured x magnitude remained exactly equal to the raw-derived value

The scatter and boxplot both receive the same 830 record-level values. No row is omitted, plotted twice, deduplicated or aggregated before plotting.

### Logarithmic-axis checks

| organism | minimum | zero_values | negative_values |
| --- | --- | --- | --- |
| Caenorhabditis elegans | 4.0 | 0 | 0 |
| Mus musculus | 7.0 | 0 | 0 |
| Drosophila melanogaster | 5.0 | 0 | 0 |
| Saccharomyces cerevisiae | 10.0 | 0 | 0 |

Every plotted magnitude is strictly greater than zero, so the log x-axis is numerically legitimate.

### Figure-file immutability during the audit

| Figure 6 file | SHA-256 before | SHA-256 after | unchanged |
| --- | --- | --- | --- |
| task\task 2\figures\exploratory\fig6_lifespan_extension_magnitude.png | cb81860b934b3217f283c002acabe201f465a615e9b514a23427c8726e32ce91 | cb81860b934b3217f283c002acabe201f465a615e9b514a23427c8726e32ce91 | True |
| task\task 2\figures\exploratory\fig6_lifespan_extension_magnitude.svg | 843a60e4fc0b3e42dacbb4cc022f09a9e464e3dcf3eef0526878196d54d352a0 | 843a60e4fc0b3e42dacbb4cc022f09a9e464e3dcf3eef0526878196d54d352a0 | True |
| task\task 2\figures\exploratory\fig6_signed_lifespan_effect_magnitude.png | 08a32b64198779eb0403ff9de96f9d3a80e1205256aa67214c61986ce95a6fa7 | 08a32b64198779eb0403ff9de96f9d3a80e1205256aa67214c61986ce95a6fa7 | True |
| task\task 2\figures\exploratory\fig6_signed_lifespan_effect_magnitude.svg | 241f3d190ac8500c54501f550eaacf09c57fa7ad70f3d535dd8de0330c7d64b0 | 241f3d190ac8500c54501f550eaacf09c57fa7ad70f3d535dd8de0330c7d64b0 | True |
| task\task 2\figures\exploratory\fig6_preview_A_extension_only.png | c606fb3b827b985a63b2866eae61966c1ea9484b6195591b8cbc0c5bb8f2a05d | c606fb3b827b985a63b2866eae61966c1ea9484b6195591b8cbc0c5bb8f2a05d | True |
| task\task 2\figures\exploratory\fig6_preview_A_extension_only.svg | 494e48d79d53b8199aeae782af6ff082d819b85ddf275f481d0ba267577bcc7d | 494e48d79d53b8199aeae782af6ff082d819b85ddf275f481d0ba267577bcc7d | True |
| task\task 2\figures\exploratory\fig6_preview_B_signed_effects.png | c70acc344c6c523e85731c552bfea7691ecd1a99d9005162b14a3166c707f71f | c70acc344c6c523e85731c552bfea7691ecd1a99d9005162b14a3166c707f71f | True |
| task\task 2\figures\exploratory\fig6_preview_B_signed_effects.svg | 06fd163151174452afc1a74bd8a6b4d458662952f8ff606cc724a2f6617a7f86 | 06fd163151174452afc1a74bd8a6b4d458662952f8ff606cc724a2f6617a7f86 | True |

## Units and interpretation

### Confirmed from source/documentation

- The downloaded CSV calls the raw field `avg lifespan change (max obsv)` but does not include a unit in the column header.
- The local `release.html` records GenAge Build 21 and 2,205 model-organism genes; it does **not include** a definition of the magnitude field or its unit.
- The official [GenAge model-organism download page](https://genomics.senescence.info/genes/models.html) identifies Build 21 and the 2,205-record download.
- The official [GenAge model-organism search](https://genomics.senescence.info/genes/search.php?lifespan_effect=increase&long_influence=&sort=5) labels the display column “Max average lifespan change” and appends `%` to these records, including age-1 1000%, SIR2 600%, SCH9 300%, daf-2 200% and Vinc 150%.
- Official detail pages also express the field under percentage-change headings, including [age-1](https://genomics.senescence.info/genes/details.php?gene=age-1&organism=elegans), [SIR2](https://genomics.senescence.info/genes/details.php?gene=SIR2&organism=cerevisiae), and [Ghrh](https://genomics.senescence.info/genes/details.php?gene=Ghrh&organism=musculus).

Therefore, `%` is justified as the source database's display unit for this field. The axis should be interpreted as **GenAge-reported percentage-coded lifespan-change summary values**, not as percentages independently recalculated from the underlying experiments.

### Assumed/inferred or unavailable

- The downloaded CSV itself does not define the exact formula behind `avg lifespan change (max obsv)`.
- The inspected local and official pages do not provide enough methodological detail to reconstruct how heterogeneous mean/median, fold-change and maximum-observation statements were normalised into this single field.
- No new interpretation of that calculation is asserted here. The record-level values remain heterogeneous descriptive summaries.

## Discrepancies

### Existing Figure 6 derived-file comparison

| file | exact match | rounding difference | discrepancy |
| --- | --- | --- | --- |
| fig6_lifespan_magnitude_data_availability.csv | 171 | 9 | 0 |
| fig6_lifespan_magnitude_extremes.csv | 795 | 5 | 0 |
| fig6_lifespan_magnitude_summary.csv | 117 | 18 | 0 |

All counts, categories, record ordering and exact-valued statistics match. Non-exact cells are explained by the six-significant-digit export format used for the existing files; these include percentages, non-integer statistics, two magnitude strings and three large Entrez-ID strings in the old extremes table. **No Figure 6 sample-size, magnitude, row-inclusion, maximum-label or plotting discrepancy was found.** The three existing Figure 6 derived files were read only and were not overwritten.

### Every non-exact (rounding-only) comparison

| file | field | key | independent | existing | classification |
| --- | --- | --- | --- | --- | --- |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Caenorhabditis elegans / Increase | 74.506579 | 74.5066 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Caenorhabditis elegans / Decrease | 66.260163 | 66.2602 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Mus musculus / Increase | 82.142857 | 82.1429 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Mus musculus / Decrease | 36.764706 | 36.7647 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Saccharomyces cerevisiae / Increase | 72.818792 | 72.8188 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Saccharomyces cerevisiae / Decrease | 42.248062 | 42.2481 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Drosophila melanogaster / Increase | 86.363636 | 86.3636 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Drosophila melanogaster / Decrease | 58.928571 | 58.9286 | rounding difference |
| fig6_lifespan_magnitude_data_availability.csv | usable_magnitude_percent | Podospora anserina / Increase | 33.333333 | 33.3333 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Saccharomyces cerevisiae / Increase | 35.788018 | 35.788 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Saccharomyces cerevisiae / Increase | 45.344476 | 45.3445 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Saccharomyces cerevisiae / Decrease | 40.018349 | 40.0183 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Saccharomyces cerevisiae / Decrease | 16.194089 | 16.1941 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Caenorhabditis elegans / Increase | 33.913466 | 33.9135 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Caenorhabditis elegans / Increase | 51.605627 | 51.6056 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Caenorhabditis elegans / Decrease | 25.177914 | 25.1779 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | median | Caenorhabditis elegans / Decrease | 22.200001 | 22.2 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Caenorhabditis elegans / Decrease | 14.244978 | 14.245 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Drosophila melanogaster / Increase | 35.514825 | 35.5148 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Drosophila melanogaster / Increase | 24.761633 | 24.7616 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Drosophila melanogaster / Decrease | 32.690909 | 32.6909 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Drosophila melanogaster / Decrease | 17.566114 | 17.5661 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | q3 | Drosophila melanogaster / Decrease | 42.599998 | 42.6 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | iqr | Drosophila melanogaster / Decrease | 22.599998 | 22.6 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | mean | Mus musculus / Increase | 22.841304 | 22.8413 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Mus musculus / Increase | 9.6615981 | 9.6616 | rounding difference |
| fig6_lifespan_magnitude_summary.csv | sd | Mus musculus / Decrease | 19.20434 | 19.2043 | rounding difference |
| fig6_lifespan_magnitude_extremes.csv | entrez gene id | row 37: Caenorhabditis elegans / Decrease / rank 7 | 3565786 | 3.56579e+06 | rounding difference |
| fig6_lifespan_magnitude_extremes.csv | entrez gene id | row 48: Drosophila melanogaster / Increase / rank 8 | 3772232 | 3.77223e+06 | rounding difference |
| fig6_lifespan_magnitude_extremes.csv | entrez gene id | row 55: Drosophila melanogaster / Decrease / rank 5 | 3772064 | 3.77206e+06 | rounding difference |
| fig6_lifespan_magnitude_extremes.csv | avg lifespan change (max obsv) | row 59: Drosophila melanogaster / Decrease / rank 9 | 42.59999847 | 42.6 | rounding difference |
| fig6_lifespan_magnitude_extremes.csv | avg lifespan change (max obsv) | row 69: Mus musculus / Increase / rank 9 | 33.29999924 | 33.3 | rounding difference |

The old extremes table was exported with a six-significant-digit float format. Consequently, four large Entrez IDs are written in scientific notation; one remains numerically exact and three lose their final one-digit precision in that old display-only table. The verified plotting-data and maximum-record CSVs preserve Entrez IDs as source strings. This is classified as a rounding difference in the pre-existing audit table, not a Figure 6 point or annotation discrepancy, because Entrez ID is not used to position points and all raw row identities, GenAge IDs, symbols and magnitudes were independently matched.

## Verdict

**PASS — Figure 6 data and annotations independently verified.**
