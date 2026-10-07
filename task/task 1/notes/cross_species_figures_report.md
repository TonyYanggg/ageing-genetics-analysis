# Cross-species descriptive figures

This stage compares `longevity influence` and `lifespan effect` across GenAge model organisms using the immutable `genage_models.csv` source. All calculations use **rows/GenAge records**; repeated gene records were not deduplicated. The full tables retain all nine organisms, while Figures 4A–5B show only the five previously validated major organisms.

Reproducible script: [`scripts/02_cross_species_descriptive_figures.py`](../scripts/02_cross_species_descriptive_figures.py)

Run from the project root with:

```text
python "task/task 1/scripts/02_cross_species_descriptive_figures.py"
```

## Longevity influence by organism

### Counts

| Organism | Anti-Longevity | Pro-Longevity | Necessary for fitness | Unannotated | Unclear | Total |
|---|---:|---:|---:|---:|---:|---:|
| *Saccharomyces cerevisiae* | 349 | 51 | 492 | 10 | 9 | 911 |
| *Caenorhabditis elegans* | 578 | 283 | 0 | 15 | 13 | 889 |
| *Drosophila melanogaster* | 73 | 125 | 0 | 1 | 3 | 202 |
| *Mus musculus* | 45 | 85 | 0 | 4 | 2 | 136 |
| *Schizosaccharomyces pombe* | 54 | 1 | 5 | 1 | 0 | 61 |
| *Podospora anserina* | 1 | 0 | 0 | 2 | 0 | 3 |
| *Caenorhabditis briggsae* | 1 | 0 | 0 | 0 | 0 | 1 |
| *Danio rerio* | 0 | 1 | 0 | 0 | 0 | 1 |
| *Mesocricetus auratus* | 0 | 0 | 0 | 1 | 0 | 1 |
| **All records** | **1,101** | **546** | **497** | **34** | **27** | **2,205** |

Source: [`results/tables/longevity_influence_by_organism_counts.csv`](../results/tables/longevity_influence_by_organism_counts.csv)

### Percentage within organism

| Organism | Anti-Longevity | Pro-Longevity | Necessary for fitness | Unannotated | Unclear |
|---|---:|---:|---:|---:|---:|
| *Saccharomyces cerevisiae* | 38.31% | 5.60% | 54.01% | 1.10% | 0.99% |
| *Caenorhabditis elegans* | 65.02% | 31.83% | 0.00% | 1.69% | 1.46% |
| *Drosophila melanogaster* | 36.14% | 61.88% | 0.00% | 0.50% | 1.49% |
| *Mus musculus* | 33.09% | 62.50% | 0.00% | 2.94% | 1.47% |
| *Schizosaccharomyces pombe* | 88.52% | 1.64% | 8.20% | 1.64% | 0.00% |
| *Podospora anserina* | 33.33% | 0.00% | 0.00% | 66.67% | 0.00% |
| *Caenorhabditis briggsae* | 100.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| *Danio rerio* | 0.00% | 100.00% | 0.00% | 0.00% | 0.00% |
| *Mesocricetus auratus* | 0.00% | 0.00% | 0.00% | 100.00% | 0.00% |

Source: [`results/tables/longevity_influence_by_organism_percent.csv`](../results/tables/longevity_influence_by_organism_percent.csv)

Data notes:

- No longevity-influence value is missing.
- Each count row equals the corresponding organism's validated record total.
- Unrounded percentages equal 100% within every organism; displayed percentages are rounded.
- The four minor organisms remain in the tables but are excluded from the figures because each has only 1–3 records.

Figures:

- Figure 4A: [`PNG`](../figures/exploratory/fig4a_longevity_influence_counts_by_species.png) · [`SVG`](../figures/exploratory/fig4a_longevity_influence_counts_by_species.svg)
- Figure 4B: [`PNG`](../figures/exploratory/fig4b_longevity_influence_percent_by_species.png) · [`SVG`](../figures/exploratory/fig4b_longevity_influence_percent_by_species.svg)

## Lifespan effect by organism

`lifespan effect` is retained as a summary of directions reported across associated genetic manipulations; it is not treated as a knockout-only variable. Blank values are represented explicitly as `Missing`.

### Counts

| Organism | Increase | Decrease | Increase and Decrease | Missing | Total |
|---|---:|---:|---:|---:|---:|
| *Saccharomyces cerevisiae* | 298 | 516 | 90 | 7 | 911 |
| *Caenorhabditis elegans* | 608 | 246 | 31 | 4 | 889 |
| *Drosophila melanogaster* | 132 | 56 | 13 | 1 | 202 |
| *Mus musculus* | 56 | 68 | 10 | 2 | 136 |
| *Schizosaccharomyces pombe* | 56 | 5 | 0 | 0 | 61 |
| *Podospora anserina* | 3 | 0 | 0 | 0 | 3 |
| *Caenorhabditis briggsae* | 1 | 0 | 0 | 0 | 1 |
| *Danio rerio* | 0 | 1 | 0 | 0 | 1 |
| *Mesocricetus auratus* | 1 | 0 | 0 | 0 | 1 |
| **All records** | **1,155** | **892** | **144** | **14** | **2,205** |

Source: [`results/tables/lifespan_effect_by_organism_counts.csv`](../results/tables/lifespan_effect_by_organism_counts.csv)

### Percentage within organism

| Organism | Increase | Decrease | Increase and Decrease | Missing |
|---|---:|---:|---:|---:|
| *Saccharomyces cerevisiae* | 32.71% | 56.64% | 9.88% | 0.77% |
| *Caenorhabditis elegans* | 68.39% | 27.67% | 3.49% | 0.45% |
| *Drosophila melanogaster* | 65.35% | 27.72% | 6.44% | 0.50% |
| *Mus musculus* | 41.18% | 50.00% | 7.35% | 1.47% |
| *Schizosaccharomyces pombe* | 91.80% | 8.20% | 0.00% | 0.00% |
| *Podospora anserina* | 100.00% | 0.00% | 0.00% | 0.00% |
| *Caenorhabditis briggsae* | 100.00% | 0.00% | 0.00% | 0.00% |
| *Danio rerio* | 0.00% | 100.00% | 0.00% | 0.00% |
| *Mesocricetus auratus* | 100.00% | 0.00% | 0.00% | 0.00% |

Source: [`results/tables/lifespan_effect_by_organism_percent.csv`](../results/tables/lifespan_effect_by_organism_percent.csv)

Data notes:

- All 14 missing lifespan-effect values remain in the tables and figures: 7 in *S. cerevisiae*, 4 in *C. elegans*, 1 in *D. melanogaster*, and 2 in *M. musculus*.
- Each count row equals the corresponding organism's validated record total.
- Unrounded percentages equal 100% within every organism; displayed percentages are rounded.
- No gene records were deduplicated.

Figures:

- Figure 5A: [`PNG`](../figures/exploratory/fig5a_lifespan_effect_counts_by_species.png) · [`SVG`](../figures/exploratory/fig5a_lifespan_effect_counts_by_species.svg)
- Figure 5B: [`PNG`](../figures/exploratory/fig5b_lifespan_effect_percent_by_species.png) · [`SVG`](../figures/exploratory/fig5b_lifespan_effect_percent_by_species.svg)

## Observed descriptive patterns

- `Necessary for fitness` accounts for 54.01% of *S. cerevisiae* records and 8.20% of *S. pombe* records; it is 0% in the other three major organisms. Of all 497 records in this category, 492 (98.99%) are *S. cerevisiae* records.
- `Anti-Longevity` represents 88.52% of *S. pombe* records and 65.02% of *C. elegans* records.
- `Pro-Longevity` represents 62.50% of *M. musculus* records and 61.88% of *D. melanogaster* records.
- `Increase` represents 91.80% of *S. pombe*, 68.39% of *C. elegans*, and 65.35% of *D. melanogaster* records.
- `Decrease` represents 56.64% of *S. cerevisiae* records and 50.00% of *M. musculus* records.
- `Increase and Decrease` is present in four major organisms, ranging from 3.49% in *C. elegans* to 9.88% in *S. cerevisiae*; it is absent from *S. pombe*.
- Missing lifespan-effect values constitute at most 1.47% of any major-organism total.

No hypothesis testing, biological explanation, mechanistic interpretation, or additional analysis was performed.
