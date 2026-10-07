# Figure 6 lifespan-effect magnitude exploration

This is a descriptive audit and two diagnostic previews, not the final Figure 6. All calculations use GenAge records (rows); repeated genes were not deduplicated.

## Source audit

- Source: `GenAge Project/models_genes/genage_models.csv`
- Raw dataset size: **2,205 rows x 8 raw columns**
- Organism field: `organism`
- Lifespan-effect field: `lifespan effect`
- Lifespan-change field: `avg lifespan change (max obsv)`

Raw columns:

- `GenAge ID`
- `symbol`
- `name`
- `organism`
- `entrez gene id`
- `avg lifespan change (max obsv)`
- `lifespan effect`
- `longevity influence`

Raw data types:

- `GenAge ID`: `int64`
- `lifespan effect`: `str`
- `organism`: `str`
- `entrez gene id`: `float64`
- `symbol`: `str`
- `avg lifespan change (max obsv)`: `float64`
- `name`: `str`
- `longevity influence`: `str`

## Data availability

| lifespan_effect_category | records | usable | missing_magnitude | usable_percent |
| --- | --- | --- | --- | --- |
| Increase | 1,155 | 832 | 323 | 72.03% |
| Decrease | 892 | 440 | 452 | 49.33% |
| Increase and Decrease | 144 | 0 | 144 | 0.00% |
| Missing/Other | 14 | 0 | 14 | 0.00% |

Across the full dataset, **1,272 of 2,205 records (57.69%)** have a usable numeric magnitude. There were **0** non-missing values that failed numeric conversion; the remaining unavailable values were blank. `Increase and Decrease` and missing/other lifespan-effect records have no usable magnitudes.

The full organism-by-category audit is in `results/tables/fig6_lifespan_magnitude_data_availability.csv`.

## Major-organism sample sizes

Usable numeric magnitude records:

| organism | Increase | Decrease | Increase and Decrease | Missing/Other |
| --- | --- | --- | --- | --- |
| Saccharomyces cerevisiae | 217 | 218 | 0 | 0 |
| Caenorhabditis elegans | 453 | 163 | 0 | 0 |
| Drosophila melanogaster | 114 | 33 | 0 | 0 |
| Mus musculus | 46 | 25 | 0 | 0 |
| Schizosaccharomyces pombe | 0 | 0 | 0 | 0 |

The five major organisms are retained in both previews. *Schizosaccharomyces pombe* has no usable magnitude values and is shown as `n=0`, rather than being silently excluded.

## Magnitude distributions and extremes

The requested mean, median, SD, quartiles, IQR, minimum, and maximum for each major organism and effect category are in `results/tables/fig6_lifespan_magnitude_summary.csv`. Category maxima are:

| organism | Increase | Decrease | Increase and Decrease |
| --- | --- | --- | --- |
| Saccharomyces cerevisiae | 600 | 85 | NA |
| Caenorhabditis elegans | 1,000 | 84 | NA |
| Drosophila melanogaster | 150 | 75 | NA |
| Mus musculus | 46 | 83 | NA |
| Schizosaccharomyces pombe | NA | NA | NA |

The largest extension-side values are 1,000 (*C. elegans*, `age-1`), 600 (*S. cerevisiae*, `SIR2`), and 300 (*S. cerevisiae*, `SCH9`). The largest observed decrease magnitude is 85 (*S. cerevisiae*, `RPT5`). These records are retained. The top ten usable records per major organism/effect category are in `results/tables/fig6_lifespan_magnitude_extremes.csv`; groups with no usable values cannot contribute ranked records.

Because the range spans 4 to 1,000, Preview A uses a logarithmic magnitude axis and Preview B uses a symmetric-log signed axis. These transformations expose the dense central distributions while keeping every extreme observation visible.

## Is the magnitude field unsigned?

Yes, for this release it behaves as an unsigned field: among 1,272 usable values, **1,272 are positive, 0 are zero, and 0 are negative**. Positive values occur in both `Increase` and `Decrease` records. The raw column was not changed.

For diagnosis only, `signed_lifespan_change` was reconstructed as positive absolute magnitude for `Increase`, negative absolute magnitude for `Decrease`, and undefined for `Increase and Decrease` or missing/other classifications.

Signed distribution by major organism:

| organism | n | mean | median | sd | q1 | q3 | iqr | minimum | maximum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Saccharomyces cerevisiae | 435 | -2.20 | -10.00 | 50.94 | -40.00 | 28.50 | 68.50 | -85.00 | 600.00 |
| Caenorhabditis elegans | 616 | 18.28 | 19.05 | 51.88 | -8.25 | 32.00 | 40.25 | -84.00 | 1000.00 |
| Drosophila melanogaster | 147 | 20.20 | 24.00 | 36.85 | 8.00 | 39.00 | 31.00 | -75.00 | 150.00 |
| Mus musculus | 71 | 2.88 | 16.00 | 30.49 | -19.00 | 21.50 | 40.50 | -83.00 | 46.00 |
| Schizosaccharomyces pombe | 0 | NA | NA | NA | NA | NA | NA | NA | NA |

## Preview comparison

- **Preview A (extension only)** is the more direct match to Pedro's wording about "life extending effects" and avoids reconstructing direction. Its main limitations are the 1,000 and 600 extremes, uneven coverage, and the complete absence of usable *S. pombe* values.
- **Preview B (signed Increase versus Decrease)** is more informative for understanding the full direction-and-magnitude structure and exposes organism-specific mixtures of extension and shortening. However, its sign is reconstructed from a separate category, mixed-direction records are undefined, and it goes beyond a literal extension-only request.

Preview A is therefore the leading interpretation if Pedro's wording is taken literally; Preview B is the stronger diagnostic companion. Neither should be treated as final until the intended scope is confirmed.

## Ambiguities and review risks

- The column name does not encode units. The values appear to be percentage-like, but the unit should be confirmed from GenAge documentation or observation-level records before final axis wording.
- `avg lifespan change (max obsv)` appears to retain a maximum observed summary rather than all underlying observations; a gene-level row may therefore conceal intervention, sex, dose, tissue, or study-level variation.
- All 144 `Increase and Decrease` records lack a numeric magnitude, so their two directions cannot be reconstructed from this table.
- *S. pombe* has 61 records but zero usable magnitude values, preventing a five-species quantitative comparison without supplementary observation-level data.
- The extreme extension values materially affect linear-scale views and mean estimates. They should be checked against the underlying GenAge observations, not automatically removed.
- The records are not independent effect estimates in a formal meta-analytic sense; no inferential test or biological ranking is justified at this stage.

## Diagnostic figure contract

- Core question: determine whether available GenAge summary values support a defensible cross-organism comparison of reported lifespan-effect magnitude.
- Archetype: quantitative comparison.
- Backend: Python/matplotlib only.
- Evidence: Preview A directly compares extension magnitudes; Preview B tests whether reconstructing direction around zero adds useful information.
- Source data: the three Figure 6 CSV audit tables produced by the companion script.
- Reviewer risk: unsigned source values, uneven missingness, mixed-direction records without magnitudes, non-independent records, and a small number of extreme observations.
