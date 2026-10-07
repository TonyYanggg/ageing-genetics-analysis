# Basic GenAge model-organism figures

This exploratory stage used the **2,205 rows/GenAge records** in `GenAge Project/models_genes/genage_models.csv`. The raw CSV was read without modification. Record count is the primary measure; repeated gene records were not merged.

Reproducible script: [`scripts/01_basic_descriptive_figures.py`](../scripts/01_basic_descriptive_figures.py)

Run from the project root with:

```text
python "task/task 1/scripts/01_basic_descriptive_figures.py"
```

The script validates the expected columns and all requested counts before writing outputs. It stops with an error if the input no longer reproduces the validated 2,205-record inventory.

## Record count versus unique gene identifiers

Unique Entrez counts exclude missing IDs.

| Organism | Records | Unique symbols | Unique Entrez IDs |
|---|---:|---:|---:|
| *Saccharomyces cerevisiae* | 911 | 911 | 911 |
| *Caenorhabditis elegans* | 889 | 887 | 886 |
| *Drosophila melanogaster* | 202 | 201 | 200 |
| *Mus musculus* | 136 | 136 | 136 |
| *Schizosaccharomyces pombe* | 61 | 61 | 61 |
| *Podospora anserina* | 3 | 2 | 2 |
| *Caenorhabditis briggsae* | 1 | 1 | 1 |
| *Danio rerio* | 1 | 1 | 1 |
| *Mesocricetus auratus* | 1 | 1 | 1 |
| **Total / sum across organisms** | **2,205** | **2,201** | **2,199** |

The small differences confirm that a few organism-specific records reuse a symbol or Entrez ID. They remain separate records in all three figures. Source table: [`results/tables/organism_record_gene_counts.csv`](../results/tables/organism_record_gene_counts.csv).

## Figure 1 — GenAge model-organism records by species

What was counted: rows grouped by `organism`.  
Denominator: 2,205 model-organism records.

| Organism | Count | Percentage |
|---|---:|---:|
| *Saccharomyces cerevisiae* | 911 | 41.32% |
| *Caenorhabditis elegans* | 889 | 40.32% |
| *Drosophila melanogaster* | 202 | 9.16% |
| *Mus musculus* | 136 | 6.17% |
| *Schizosaccharomyces pombe* | 61 | 2.77% |
| *Podospora anserina* | 3 | 0.14% |
| *Caenorhabditis briggsae* | 1 | 0.05% |
| *Danio rerio* | 1 | 0.05% |
| *Mesocricetus auratus* | 1 | 0.05% |

Data issue: organism is complete, but record counts are slightly larger than unique-symbol/Entrez counts for three organisms. Percentages are rounded independently and therefore sum to 100.03%.

Outputs:

- [`figures/exploratory/fig1_records_by_species.png`](../figures/exploratory/fig1_records_by_species.png)
- [`figures/exploratory/fig1_records_by_species.svg`](../figures/exploratory/fig1_records_by_species.svg)
- [`results/tables/fig1_records_by_species.csv`](../results/tables/fig1_records_by_species.csv)

## Figure 2 — Longevity influence

What was counted: rows grouped by `longevity influence`.  
Denominator: 2,205 model-organism records.

| Classification | Count | Percentage |
|---|---:|---:|
| Anti-Longevity | 1,101 | 49.93% |
| Pro-Longevity | 546 | 24.76% |
| Necessary for fitness | 497 | 22.54% |
| Unannotated | 34 | 1.54% |
| Unclear | 27 | 1.22% |

Data issue: no values are missing. Percentages are rounded independently and sum to 99.99%.

Outputs:

- [`figures/exploratory/fig2_longevity_influence.png`](../figures/exploratory/fig2_longevity_influence.png)
- [`figures/exploratory/fig2_longevity_influence.svg`](../figures/exploratory/fig2_longevity_influence.svg)
- [`results/tables/fig2_longevity_influence.csv`](../results/tables/fig2_longevity_influence.csv)

## Figure 3 — Reported lifespan-effect direction

What was counted: rows grouped by `lifespan effect`, with blank values retained as `Missing`. This field is treated as a summary of reported directions across associated genetic manipulations, not as a knockout-only response.  
Denominator: 2,205 model-organism records, including missing values.

| Lifespan-effect direction | Count | Percentage |
|---|---:|---:|
| Increase | 1,155 | 52.38% |
| Decrease | 892 | 40.45% |
| Increase and Decrease | 144 | 6.53% |
| Missing | 14 | 0.63% |

Data issue: 14 records have no lifespan-effect direction. They are displayed explicitly rather than excluded. Percentages are rounded independently and sum to 99.99%.

Outputs:

- [`figures/exploratory/fig3_lifespan_effect.png`](../figures/exploratory/fig3_lifespan_effect.png)
- [`figures/exploratory/fig3_lifespan_effect.svg`](../figures/exploratory/fig3_lifespan_effect.svg)
- [`results/tables/fig3_lifespan_effect.csv`](../results/tables/fig3_lifespan_effect.csv)

## Validation summary

- All organism counts match the previously inspected inventory.
- Longevity-influence counts sum to 2,205 and match the expected categories.
- Lifespan-effect counts, including 14 missing values, sum to 2,205 and match the expected categories.
- No deduplication, biological interpretation, lifespan-magnitude analysis, hypothesis testing, or additional analysis was performed.
