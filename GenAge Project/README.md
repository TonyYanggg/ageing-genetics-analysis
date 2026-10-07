# GenAge source data

This folder contains the downloaded GenAge datasets used as shared source material across the project. The local release notes identify them as **GenAge Build 21**, released on **28 August 2023**.

## Contents

| Folder | Downloaded files | Description |
|---|---|---|
| `human_genes/` | `genage_human.csv`, `release.html` | Catalogue of 307 curated human ageing-related genes with identifiers and inclusion-reason categories |
| `models_genes/` | `genage_models.csv`, `release.html` | Catalogue of 2,205 model-organism GenAge records, including organism, lifespan-effect direction, longevity influence, and a partly populated lifespan-change field |
| `models_genes/` | `genage_models_orthologs_export.tsv` | Existing frozen human-targeted website export: 1,055 mapping rows; mapping version and retrieval date unverified; Task 4 audits and quarantines four conflicting-identity rows |
| `signatures_supplement/` | `ageing_signatures.xlsx` | Build 21 ageing-expression signature supplement; it has not yet been incorporated into the current analyses |

The detailed local audit is in the [Task 1 data inventory](../task/task%201/notes/data_inventory.md).

## Modification policy

The downloaded CSV, XLSX, and HTML files are immutable source files.

- Do not edit, normalize, reformat, rename, move, or overwrite them.
- Do not deduplicate or correct records in place.
- Read them directly for reproducible analysis.
- Write cleaned data, mappings, extracted records, and other derived files into the relevant numbered task's `results/` folder.
- Record the GenAge build and input path in scripts and reports.
- Updating this README is allowed; adding a new source download requires documenting its origin and release/build information here.

## Important interpretation notes

- The model-organism CSV represents GenAge records and is not strictly one row per unique gene.
- `lifespan effect` and `longevity influence` describe different concepts and should remain separate.
- The lifespan-change field is missing for many records and does not include observation-level experimental details.
- The stable Build 21 CSV files do not include full reference lists, manipulation details, phenotypes, or explicit cross-species orthology mappings. The separate website TSV supplies human-to-model mappings with limited provenance. Its Windows download metadata identifies `https://genomics.senescence.info/genes/search.php?organism=Homo%20sapiens&show=5&export=true&homologs=true`; this does not establish a mapping release or confidence level.

## Current use

- Task 1 reads `models_genes/genage_models.csv` for descriptive Figures 1–5.
- Task 2 reads the same file for the lifespan-extension comparison.
- Task 4 will use the complete 2,205-record model CSV across all nine organisms plus a versioned model-to-model orthology source. Its earlier use of the human CSV and frozen human-targeted TSV is retained only as a documented exploratory screen.
- Future timeline, homology, and human-relevance tasks may combine these downloads with separately documented web-derived or external data.
