# GenAge data inventory

Inspection date: 2026-08-18  
Project scope: GenAge human genes and GenAge model-organism genes only

## Scope and method

This report inventories the downloaded files without modifying the raw sources. It is an inspection-stage document only: no biological inference, publication figure, modelling, or data cleaning has been performed.

The directory was scanned recursively. The two in-scope CSV files were parsed with a standards-compliant CSV reader and profiled for dimensions, inferred types, missingness, distinct values, duplicated records, whitespace, and encoding. Row counts below exclude the header row. Percentages use the number of data rows in the relevant file as the denominator.

The local GenAge release notes identify the data as **GenAge Build 21**, released **28 August 2023**, with **307 human genes** and **2,205 model-organism genes**. The two copies of `release.html` under the human and model folders are byte-identical. The observed CSV row counts match the release notes exactly.

Primary local documentation:

- [`GenAge Project/human_genes/release.html`](../../../GenAge%20Project/human_genes/release.html)
- [`GenAge Project/models_genes/release.html`](../../../GenAge%20Project/models_genes/release.html)
- HAGR's linked [GenAge help page](https://genomics.senescence.info/help.html), used only where the local release note does not define a category

## Overview of project files

The project contained **31 files in 15 directories** before this report was created.

### In-scope GenAge datasets

| Dataset | Path | Format | Size | Data rows | Columns | Role |
|---|---|---:|---:|---:|---:|---|
| Human genes | [`genage_human.csv`](../../../GenAge%20Project/human_genes/genage_human.csv) | CSV | 22,778 bytes (22.24 KiB) | 307 | 6 | Curated human ageing-related gene catalogue |
| Model-organism genes | [`genage_models.csv`](../../../GenAge%20Project/models_genes/genage_models.csv) | CSV | 214,711 bytes (209.68 KiB) | 2,205 | 8 | Curated genes associated with ageing or longevity in model organisms |

### GenAge-associated but outside the present scope

| File | Classification | Reason not profiled here |
|---|---|---|
| [`ageing_signatures.xlsx`](../../../GenAge%20Project/signatures_supplement/ageing_signatures.xlsx), 52,781 bytes (51.54 KiB) | GenAge Build 21 ageing-expression signature supplement. It contains a source sheet and global overexpressed/underexpressed gene tables. | The requested scope is the human-gene and model-organism-gene catalogues. This supplement should be inventoried separately if expression signatures are added to the analysis plan. |

### Files belonging to other HAGR resources

These files are not inputs for the present GenAge analysis:

| Resource | Files/folders identified | Basis for classification |
|---|---|---|
| CellAge | `Other downloads from HAGR/cellAge/cellage3.tsv`; its `release.html` | The release note identifies CellAge Build 3 and describes curated cellular-senescence genes. |
| CellAge signatures | `Other downloads from HAGR/cellSignatures/signatures1.csv`; its `release.html` | The release note identifies a CellAge cellular-senescence expression signature, not a GenAge signature. |
| DrugAge | `Other downloads from HAGR/dataset (1)_DrugAge/drugage.csv`; its `release.html` | The release note identifies DrugAge Build 5. |
| AnAge | `Other downloads from HAGR/dataset (2)_AnAge/anage_data.txt`; its `release.html` | The release note identifies AnAge Build 15 and the text file contains species life-history fields. |
| GenDR | `Other downloads from HAGR/dataset_GenDR/gendr_manipulations.csv`; its `release.html`; `Other downloads from HAGR/DR gene expression signature.xls` | The release note and filenames identify dietary-restriction gene manipulations and a DR expression signature. |
| LongevityMap | `Other downloads from HAGR/longevity_genes/longevity.csv`; its `release.html` | The release note identifies LongevityMap Build 3 and the CSV contains association, population, variant, and gene fields. |
| ARCT software | `Other downloads from HAGR/ARCT-0.9/` | Thirteen Perl source/example/build files. The README identifies ARCT as an ageing-research computational toolkit, not a downloaded biological table. |

## Dataset 1: GenAge human genes

### File summary

- **Filename:** `genage_human.csv`
- **Format:** comma-delimited text; valid UTF-8 without a byte-order mark; CRLF line endings
- **Size:** 22,778 bytes (22.24 KiB)
- **Shape:** 307 rows × 6 columns
- **Column names:** `GenAge ID`, `symbol`, `name`, `entrez gene id`, `uniprot`, `why`
- **Complete duplicated rows:** 0
- **Observed table grain:** one human GenAge entry per row

### Column inventory and interpretation

| Column | Inferred data type | Missing, n (%) | Unique non-missing | Representative values / range | Plain-English interpretation | Cleaning or interpretation notes |
|---|---|---:|---:|---|---|---|
| `GenAge ID` | Integer identifier (`int64` on import) | 0 (0.00%) | 307 | 1–309 | Internal identifier assigned to the human GenAge entry. | Unique in this file, but not contiguous. Treat as an identifier, not a measurement or row number. HAGR documentation indicates that HAGR IDs may change between releases, so it should not be the sole long-term join key. |
| `symbol` | Text identifier | 0 (0.00%) | 307 | `GHR`, `TP53`, `TERC`, `TERT` | Human gene symbol, apparently following HGNC nomenclature for this release. | Unique here. Preserve as text and record the GenAge build because symbols can be updated between database releases. |
| `name` | Text | 0 (0.00%) | 307 | `growth hormone receptor`; `tumor protein p53`; `telomerase RNA component` | Descriptive human gene name. | Useful for labels and checking identifiers, but not a stable key. |
| `entrez gene id` | Integer identifier (`int64` on import) | 0 (0.00%) | 307 | 2–493,856 | NCBI Entrez Gene identifier corresponding to the human gene. | Unique here. Import as text or integer identifier rather than a continuous numeric variable. |
| `uniprot` | Text identifier | 1 (0.33%) | 306 | `GHR_HUMAN`, `P53_HUMAN`, `TERT_HUMAN` | UniProt entry mnemonic associated with the gene/protein. | The only missing value is for `TERC` (GenAge ID 7). Its gene name indicates an RNA component, but whether the missing UniProt value is intentionally inapplicable should be confirmed rather than assumed. |
| `why` | Multi-label categorical text | 0 (0.00%) | 34 exact combinations | `putative`; `functional`; `mammal`; `downstream`; `functional,putative`; `mammal,model,cell` | One or more evidence/reason categories explaining why the human gene was included in GenAge. | Multiple labels are stored in one comma-separated cell. Parse to a long table or separate Boolean flags in a derived dataset. The value `human_link` uses an underscore, whereas prose documentation uses “human link”. |

### `why` categories

The most frequent **exact cell values** are `putative` (53), `functional` (51), `mammal` (43), `downstream` (43), `upstream` (26), `model` (18), and `functional,putative` (13). There are 34 distinct combinations.

After splitting comma-separated cells, nine atomic labels occur:

| Label | Rows containing label | Meaning given or supported by HAGR documentation |
|---|---:|---|
| `human` | 3 | Direct evidence linking the gene product to ageing in humans |
| `mammal` | 62 | Direct evidence in a mammalian model organism |
| `model` | 33 | Direct evidence in a non-mammalian model organism |
| `cell` | 25 | Direct evidence in a cellular model system |
| `human_link` | 5 | Evidence linking the gene/product to human longevity or multiple age-related phenotypes |
| `upstream` | 35 | Evidence that the gene product regulates genes previously linked to ageing |
| `functional` | 81 | Evidence linking the gene product to an ageing-related pathway or mechanism |
| `downstream` | 56 | Evidence that the product acts downstream of an ageing-related pathway, mechanism, or gene product |
| `putative` | 80 | Indirect or inconclusive evidence linking the product to ageing |

These categories describe the stated basis for inclusion; they should not be treated as equivalent evidence strengths without an explicit, documented rule.

## Dataset 2: GenAge model-organism genes

### File summary

- **Filename:** `genage_models.csv`
- **Format:** comma-delimited text; valid UTF-8 without a byte-order mark; CRLF line endings
- **Size:** 214,711 bytes (209.68 KiB)
- **Shape:** 2,205 rows × 8 columns
- **Column names:** `GenAge ID`, `symbol`, `name`, `organism`, `entrez gene id`, `avg lifespan change (max obsv)`, `lifespan effect`, `longevity influence`
- **Complete duplicated rows:** 0
- **Observed table grain:** a GenAge model-organism record per row. The file is not strictly one row per unique symbol or Entrez ID.

### Column inventory and interpretation

| Column | Inferred data type | Missing, n (%) | Unique non-missing | Representative values / range | Plain-English interpretation | Cleaning or interpretation notes |
|---|---|---:|---:|---|---|---|
| `GenAge ID` | Integer identifier (`int64`) | 0 (0.00%) | 2,205 | 1–2,271 | Internal identifier for the GenAge model-organism record. | Unique but not contiguous. Do not treat it as a measurement or assume that IDs are stable across releases. |
| `symbol` | Text identifier | 0 (0.00%) | 2,194 | `aak-2`, `Adcy5`, `SIR2`, `CG3776` | Organism-specific gene symbol or locus name. | Not unique: 11 repeated occurrences beyond the first. A symbol alone is not a safe key across species, and four symbol–organism pairs also repeat. |
| `name` | Text | 2 (0.09%) | 1,715 | `AMP-Activated Kinase`; `adenylate cyclase 5`; `Transaldolase` | Gene name or descriptive annotation. | Missing for `sams-2` (ID 737) and `Y46G5A.6` (ID 750). Twenty-one non-missing values have leading or trailing whitespace; clean only in a derived copy. Names are not unique and sometimes vary in capitalization or description. |
| `organism` | Categorical text | 0 (0.00%) | 9 | *Caenorhabditis elegans*; *Saccharomyces cerevisiae*; *Drosophila melanogaster*; *Mus musculus* | Scientific name of the model organism associated with the record. | Suitable for grouping. Preserve exact strings in the raw data; a derived taxonomic mapping could be added later if needed. |
| `entrez gene id` | Nullable integer identifier (imported as `float64` because of missing cells) | 3 (0.14%) | 2,199 | 2,539–101,830,758 | NCBI Entrez Gene identifier for the organism-specific gene. | Missing for `eat-1` (ID 202), `hsp70` (ID 352), and `SNCA` (ID 618). Use a nullable integer or text type to avoid `.0` suffixes/scientific notation. Three repeated non-missing IDs occur. |
| `avg lifespan change (max obsv)` | Numeric (`float64`) | 933 (42.31%) | 167 | 4–1,000 | Appears to record the maximum observed average-lifespan change for the entry; values appear to be percentages because HAGR displays this field as a percent. | The abbreviated header does not fully define the calculation, denominator, or whether “average” can include median results. The sign is not stored here; direction is in `lifespan effect`. High missingness is a major limitation. Seventy-eight values contain floating-point-looking precision such as `33.29999924`; round only in a derived display field after the intended precision is confirmed. The value 1,000 is extreme and should be source-checked, not automatically removed. |
| `lifespan effect` | Categorical text | 14 (0.63%) | 3 | `Increase`; `Decrease`; `Increase and Decrease` | Direction(s) of the reported lifespan effect across observations associated with the gene. | Do not infer a direction from `longevity influence` when this field is missing. The relationship between this summary and individual interventions is not available in the CSV. |
| `longevity influence` | Categorical text | 0 (0.00%) | 5 | `Anti-Longevity`; `Pro-Longevity`; `Necessary for fitness`; `Unclear`; `Unannotated` | Curated interpretation of whether the normal gene activity tends to promote or oppose longevity, or whether evidence is insufficient/ambiguous. | This is conceptually different from `lifespan effect`: an intervention that increases lifespan can imply an anti-longevity role for the unmanipulated gene. The CSV does not include the manipulation needed to reconstruct that logic. |

### Representative categorical distributions

#### Organism

| Organism | Rows |
|---|---:|
| *Saccharomyces cerevisiae* | 911 |
| *Caenorhabditis elegans* | 889 |
| *Drosophila melanogaster* | 202 |
| *Mus musculus* | 136 |
| *Schizosaccharomyces pombe* | 61 |
| *Podospora anserina* | 3 |
| *Mesocricetus auratus* | 1 |
| *Danio rerio* | 1 |
| *Caenorhabditis briggsae* | 1 |

#### Lifespan effect

| Value | Rows | Percentage of all rows |
|---|---:|---:|
| `Increase` | 1,155 | 52.38% |
| `Decrease` | 892 | 40.45% |
| `Increase and Decrease` | 144 | 6.53% |
| Missing | 14 | 0.63% |

#### Longevity influence

| Value | Rows | Percentage of all rows |
|---|---:|---:|
| `Anti-Longevity` | 1,101 | 49.93% |
| `Pro-Longevity` | 546 | 24.76% |
| `Necessary for fitness` | 497 | 22.54% |
| `Unannotated` | 34 | 1.54% |
| `Unclear` | 27 | 1.22% |

HAGR documentation describes pro-longevity genes as genes whose decreased expression reduces lifespan and/or whose overexpression extends lifespan; anti-longevity genes show the converse pattern. `Necessary for fitness` is used particularly for yeast results where shortened lifespan is not sufficient to establish accelerated ageing. `Unclear` and `Unannotated` denote conflicting or insufficient evidence, respectively.

## Data-quality observations

### File integrity and encoding

- Both CSVs are valid UTF-8, have no byte-order mark, use CRLF line endings, and contain no non-ASCII characters.
- Every parsed row has the expected number of columns.
- Both files contain values with commas, including multi-label evidence and gene descriptions. They should always be read with a CSV parser rather than split manually on commas.
- No complete duplicated row was detected in either CSV.

### Missingness

- The human table is nearly complete: only one UniProt value is missing.
- The model table has limited missingness in `name` (2 rows), `entrez gene id` (3 rows), and `lifespan effect` (14 rows).
- `avg lifespan change (max obsv)` is missing in 933 of 2,205 rows (42.31%). Any later summary of lifespan-change magnitude will therefore describe a selected subset and should report its denominator.

### Repeated identifiers and potentially conflicting records

The model file has no identical full rows, but it is not a simple unique-gene table:

- `GenAge ID` is unique.
- Eleven rows repeat a symbol used earlier in the file.
- Four symbol–organism pairs have a repeated occurrence.
- Three non-missing Entrez IDs repeat an earlier Entrez ID.

Important repeated groups include:

| Gene records | Observation |
|---|---|
| `CG3776`, *D. melanogaster*, IDs 857 and 858, Entrez 37996 | Same symbol/organism/Entrez, but different lifespan-change values and opposite `Anti-Longevity` versus `Pro-Longevity` classifications. This may represent different manipulations or evidence, but the CSV lacks the detail needed to resolve it. |
| `cox1`, *P. anserina*, IDs 137 and 146, Entrez 802462 | Same symbol/organism/Entrez with slightly different gene-name text. This is a candidate duplicate or alternate annotation requiring source-page review. |
| `eat-18`, *C. elegans*, IDs 203 and 2066 | Same symbol and organism but different Entrez IDs and descriptions. This may reflect identifier/annotation changes or distinct loci; do not merge automatically. |
| `mir-238`, *C. elegans*, IDs 2168 and 2169 | Same symbol and organism but different Entrez IDs, effect directions, and longevity classifications. Do not merge without provenance. |
| `tald-1` and `Y24D9A.8`, *C. elegans*, Entrez 177202 | Different symbols sharing one Entrez ID, consistent with possible synonym or annotation changes. |

### Fields requiring parsing or careful typing

- Parse human `why` as a multi-label field. Retain the original value and create derived atomic flags or a long evidence table.
- Import Entrez IDs as strings or nullable integers, not ordinary floating-point measurements.
- Treat `GenAge ID`, gene symbols, UniProt values, and Entrez IDs as identifiers.
- Strip leading/trailing whitespace from the 21 affected model-gene `name` values only in a derived cleaned table.
- Preserve the raw lifespan-change values, but decide and document a display/analysis precision before rounding the 78 floating-point-like values.
- Keep `lifespan effect` separate from `longevity influence`; they summarize different concepts.

### Information not present in the two CSV exports

The following scientifically useful details are **not explicit columns** in the local human/model gene CSVs:

- intervention or manipulation type, such as knockout, knockdown, mutation, or overexpression;
- detailed phenotype or experimental condition;
- strain, sex, tissue, dose, or assay type;
- publication identifiers or reference lists;
- evidence descriptions at the individual-study level;
- human homologues of model-organism genes;
- explicit orthology mappings or cross-species gene families.

The human `why` column provides a high-level evidence category, and the model table provides summarized lifespan direction and longevity influence, but neither replaces study-level provenance. These missing dimensions may exist on individual HAGR web entries or in a richer export; that should be checked before designing analyses that require them.

## Potential analysis dimensions

The following variables are available for later descriptive statistics and figures. This is only an availability map, not a scientific analysis plan.

| Dimension | Human genes | Model-organism genes | Availability note |
|---|---|---|---|
| Gene identity | `symbol`, `name`, `entrez gene id`, `uniprot`, `GenAge ID` | `symbol`, `name`, `entrez gene id`, `GenAge ID` | Strong basic identifier coverage, with the duplicate/missing-ID caveats above. |
| Organism/species | Implicitly human; no organism column | `organism` with 9 species | Supports counts and cross-tabulation by organism. |
| Evidence/reason for inclusion | `why` multi-label categories | Not provided as an evidence field | Human evidence labels can be counted after parsing. |
| Pro- versus anti-longevity role | Not present | `longevity influence` | Includes pro-, anti-, necessary-for-fitness, unclear, and unannotated categories. |
| Lifespan-effect direction | Not present | `lifespan effect` | Increase, decrease, both, or missing. |
| Lifespan-change magnitude | Not present | `avg lifespan change (max obsv)` | Numeric but 42.31% missing; definition/precision need confirmation. |
| Intervention/manipulation type | Not present | Not present | Cannot be analysed from these CSVs alone. |
| Detailed phenotype | Not present | Not present beyond summarized lifespan fields | A richer source would be required. |
| References/evidence provenance | Only broad `why` labels | Not present | Study-level claims cannot be traced from these exports alone. |
| Human homologues | Not applicable within human table | Not present | Requires a separate orthology/homology source or richer GenAge export. |
| Cross-species information | Could be joined cautiously by curated identifiers | Organism field is available, but no explicit orthology | Gene-symbol matching alone is insufficient for homologous-gene comparisons. |
| Identifier completeness/coverage | Nearly complete | Can summarize missing Entrez IDs and repeated mappings | Useful as quality-control reporting rather than biological inference. |

## Questions / uncertainties

1. What is the exact definition and unit of `avg lifespan change (max obsv)`? In particular, does “average” include median lifespan, what is the denominator, and how is the maximum chosen when several interventions or studies exist?
2. Should multiple model-organism records for the same symbol or Entrez ID be retained as separate evidence records, or does HAGR provide a canonical gene-level grouping elsewhere?
3. What manipulation or experimental observation produced each `lifespan effect` and `longevity influence` classification? The CSV does not contain enough detail to audit apparently conflicting cases such as `CG3776`.
4. Are `Unclear`, `Unannotated`, and `Necessary for fitness` intended to be included in the same summaries as pro- and anti-longevity genes, or reported separately?
5. Are the comma-separated human `why` labels all co-equal evidence categories, or is one intended to be the primary reason for inclusion?
6. Which identifier releases were used for HGNC symbols, Entrez IDs, and UniProt mnemonics in Build 21? Identifier versions matter for later joins.
7. Is the missing UniProt value for `TERC` intentionally not applicable, and how should non-protein-coding genes be represented in protein-based summaries?
8. Should the apparent whitespace and floating-point artifacts be normalized in a separately generated analysis-ready table? The raw CSVs should remain untouched.
9. Are detailed references, interventions, phenotypes, and human homologues downloadable from HAGR, or would they have to be extracted from individual GenAge entry pages?
10. Will the GenAge ageing-expression signature supplement be part of a later phase? It is GenAge-associated, but it is deliberately outside the present human/model gene catalogue audit.

## Inspection-stage conclusion

The two core files for the stated project scope are `genage_human.csv` and `genage_models.csv`. Their row counts agree with GenAge Build 21 documentation, their encodings and basic tabular structures are sound, and neither contains exact duplicated rows. Before scientific analysis, the main preparation issues are multi-label parsing in `why`, identifier typing and duplicate handling, model-gene name whitespace, ambiguous floating-point precision, and substantial missingness plus incomplete documentation for the lifespan-change field. No raw source file has been altered.
