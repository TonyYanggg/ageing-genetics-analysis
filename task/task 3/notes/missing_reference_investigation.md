# Task 3 — investigation of the two entries without GenAge references

## Decision

Keep the current timeline as the **database-derived analysis** and keep both entries undated there. Store the two primary papers below as an explicit external-curation layer. If the curated layer is shown, label it as a sensitivity analysis rather than silently mixing it into the GenAge/LibAge extraction.

This preserves the distinction between what the database supplies and what was recovered by literature searching.

## Recovered papers

| GenAge entry | Externally recovered primary paper | Match assessment |
|---|---|---|
| `ZC395.10` / `daf-41` (ID 2081) | Horikawa M, Sural S, Hsu AL, Antebi A. *Co-chaperone p23 regulates C. elegans Lifespan in Response to Temperature.* PLoS Genetics (2015). PMID: [25830239](https://pubmed.ncbi.nlm.nih.gov/25830239/); DOI: [10.1371/journal.pgen.1005023](https://doi.org/10.1371/journal.pgen.1005023) | High confidence. The paper explicitly identifies `daf-41/ZC395.10` and reports longer lifespan at 25 °C and shorter lifespan at 15 °C. |
| `Mir29a` (ID 2173) | Takeda T, Tanabe H. *Lifespan and reproduction in brain-specific miR-29-knockdown mouse.* Biochemical and Biophysical Research Communications (2016). PMID: [26902119](https://pubmed.ncbi.nlm.nih.gov/26902119/); DOI: [10.1016/j.bbrc.2016.02.055](https://doi.org/10.1016/j.bbrc.2016.02.055) | High confidence that this is the source record: its table exactly matches the GenAge summary (male mean +14.2%, maximum +7.1%; female mean −11.1%, maximum −5.6%). Caveat: the experiment targets brain miR-29, so treating it as uniquely specific to `Mir29a` is a GenAge curation choice. |

The structured record is `../results/verified_timeline/external_reference_candidates.csv`. Neither paper occurs in the collected LibAge bibliography, and neither entry has an observation/reference block on its current GenAge page or in the pinned HAGR-R observation export.

## Sensitivity effect

If these papers are used to supply 2015 and 2016 respectively:

| Quantity | Database-derived | With external curation |
|---|---:|---:|
| Dated entries | 2,203 / 2,205 | 2,205 / 2,205 |
| Distinct publications | 708 | 710 |
| C. elegans entries first referenced in 2015 | 28 | 29 |
| Mouse entries first referenced in 2016 | 3 | 4 |
| Final C. elegans cumulative count | 888 | 889 |
| Final mouse cumulative count | 135 | 136 |

Only 2 of 2,205 entries (0.09%) change. The overall interpretation and visible shape of the curves do not change; only the organism endpoints and the 2015/2016 annual counts increase by one.

## Recommended presentation

Use the database-derived figure as the main result because it is fully reproducible from the declared GenAge/LibAge sources. Report `2,203 of 2,205 entries dated (99.91%)` and name the two missing database links. If complete coverage is useful, include the externally curated result in a caption, note, or sensitivity panel and cite both primary papers.
