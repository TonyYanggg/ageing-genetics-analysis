# Resolving the orthology-group problem

Historical direct-ID baseline: this note explains the change to explicit pairs. Current counts, stable-ID recovery, experimental review and captions are in [the refined review note](refined_review.md). Do not use the 219-pair baseline below as the current result.

## Method decision

Use **explicit model-to-model ortholog pairs**, not all cross-species combinations within ancient OrthoDB groups, for the current primary evidence table. OMA release All.May2026 provides the main pairwise calls. This is a change of analysis unit, not a claim that OMA is infallible or that all OrthoDB groups are erroneous. Do not transitively merge these pairs into new orthogroups.

Exact NCBI Gene cross-references plus matching reference-organism taxonomy identify proteins. Retain multiple isoforms; collapse repeated protein pairs to a unique unordered pair of species+Entrez gene keys. Exclude a gene if one of its mapped proteins links to multiple GenAge gene keys. Preserve all excluded records in the mapping audit. Relation cardinality is oriented from species1 to species2 in the derived tables; raw API rows retain original orientation. No human genes or anchors are used.

All nine organisms remain in the inventory. OMA represents seven; Podospora and hamster are source gaps, not negative results. Ensembl checks are separately labelled supplementary evidence and are not used to merge groups or silently inflate OMA coverage. The source-specific pair table is not a completed all-nine survey.

## Results and interpretation

- Exact, unambiguous mappings: 2,073 distinct GenAge gene keys.
- 219 explicit gene pairs have qualifying lifespan evidence on both sides; 136 remain when both sides require Pro/Anti annotations.
- Broad-pair influence summaries: 70 Anti, 33 Pro, 33 discordant/mixed and 83 unresolved. These are annotation comparisons, not intervention-matched conserved effects.
- C. briggsae Cbr-daf-2 is recovered, with pairs to C. elegans daf-2, fly InR, and mouse Igf1r and Insr.
- The Coq2/HSP104 pair is not called by OMA. Other supported relationships from that old group are retained as explicit pairs, rather than discarding the whole family.
- Of 181 old recurrent source groups, 74 have all enumerated lifespan-evidence pairs assessable and called by OMA; 99 contain at least one assessable pair without an OMA call; 14 have at least one unassessable pair. Categories overlap except the fully-supported category. No OMA call is not proof of non-orthology.

Independent verification reconstructs the 219 broad and 136 directional pairs from raw cached API responses, covering all 21 genome pairs and 124,274 source protein pairs. It checks source hashes, pagination totals, duplicate collapse and record provenance. This verifies computation, not biological ground truth.

## Figure contract and caption

Core conclusion: explicit OMA ortholog pairs link lifespan-associated genes across represented model organisms, while coverage limitations remain visible across the complete nine-organism inventory. Quantitative-grid archetype: panel a gives exact-identifier mapping coverage; panel b is the main species-pair evidence matrix. Continue the existing Python workflow. Export at 183 x 145 mm, editable SVG/PDF plus PNG and 600-dpi TIFF. Counts are descriptive catalogue counts, not independent biological replicates; no hypothesis test, error bar or enrichment claim applies.

**Figure | Explicit ortholog-pair evidence and source coverage.** a, Unique GenAge gene keys with unambiguous exact-ID mappings divided by identifiable gene keys in each organism. NA means the organism is absent from this OMA release, not zero conservation. b, Number of unique gene pairs with source-supported orthology and a recorded lifespan effect for both genes. Each pair appears twice in the symmetric matrix, but is counted once in the total of 219; genes can participate in multiple pairs. Grey NA cells are unassessed source gaps. Numerical zero means no qualifying pair recovered under this mapping procedure, not absence of an evolutionary or longevity relationship. GenAge Build 21 and OMA All.May2026. Ensembl supplementary checks are excluded from this source-specific figure. Sources: coverage.csv, longevity_pairs.csv and figure_pair_counts.csv in the paired-output directory.

## Remaining work

Ensembl 116 supplementary checks establish hamster Csnk1e (Entrez 101830758) to mouse 27373 and zebrafish 100006858 as one-to-one ortholog calls, but neither target has a GenAge lifespan entry here. Podospora Entrez 6194500 and 802462 yield no direct Ensembl cross-reference hits. The observed OrthoDB identity bridge 6194500 → UniProt B2B296 also returns no Ensembl hit; no identity or orthology is inferred from these missing mappings. OMA absence was checked against names, taxon IDs and relevant synonyms, not merely one species-name spelling.

Finish exact-identifier/orthology recovery for Podospora and hamster beyond the supplementary comparisons; audit other unmapped and ambiguous genes. Experimental intervention review and any defensible multi-species group presentation follow the pairwise evidence, not unrestricted graph connectivity. The old OrthoDB atlas remains a historical source audit and is not the current main result.

## Method sources

- [OMA orthology definitions](https://omabrowser.org/oma/homologs/)
- [OMA data access and SPARQL](https://omabrowser.org/oma/uses/)
- [Ensembl orthology API](https://rest.ensembl.org/documentation/info/homology_species_gene_id)
