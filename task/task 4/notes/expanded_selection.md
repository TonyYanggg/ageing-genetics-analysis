# Selection of the expanded biological illustration — 25 September 2026

Historical selection rationale. On 26 September the user requested the same style for all genes; the [all-gene atlas](all_gene_atlas_legend.md) now includes every qualifying gene and pair. The selection below applies only to the retained six-family companion, not atlas inclusion.

The user approved expanding Pedro's three-row template with more ortholog families in the same four organisms. Pedro's earlier approval supports the format; he has not approved these new selections. The primary review deliverable is now the expanded biological illustration. The comprehensive 206-pair catalogue and overview support selection and source auditing.

## Figure contract

Claim: selected orthologous genes have recurring lifespan evidence across major models, with effects that depend on intervention and context. The evidence grid links each cell to a primary experiment and each row to explicitly enumerated OMA pairs. Archetype: qualitative evidence grid, with the gene cards as the central evidence and edge/caveat strips underneath. Use the existing Python workflow, 183 mm width, editable SVG/PDF, PNG and 600-dpi TIFF. Six rows fit at 247 mm height, allowing one continuous comparison rather than two separate sheets. Colour indicates an experiment's outcome, never a gene's inherent good/bad role. The figure has no new inferential statistics or pooled effect sizes.

## Reproducible candidate screen

All 206 frozen four-organism pairs were screened. Connected networks partition their 320 gene nodes into 133 candidate sets, of which 22 span at least three species. The network operation is solely a way to find candidates; it is not orthogroup inference, and it adds no direct relationships. Multiple GenAge IDs for one source gene are retained together. Every component, including all 111 two-species candidates, is in `results/expanded/candidate_screen.csv`.

Preferred candidates have lifespan evidence in at least three species, explicit direct relationships among displayed genes and interpretable primary experiments. Selection then balances continuity with the approved figure, mechanistic breadth and informative caveats. Directional agreement is not an eligibility filter. No numerical score or 'top six' claim is made.

| Selected family | Why included | Cells / direct pairs |
|---|---|---|
| TOR | Approved example; spans all four organisms | 4 / 6 |
| S6 kinase | Approved example; preserves the sex-specific mouse result | 3 / 2 |
| Insulin/IGF receptor | Approved example; preserves receptor paralogy and replication conflict | 3 / 2 |
| AMPK | Adds energy sensing and a contrast between selected yeast and animal outcomes | 3 / 2 |
| ATG7 | Adds autophagy with traceable gene-level experiments across yeast, worm and fly | 3 / 3 |
| SOD2 | Adds oxidative stress and all-four coverage with context-dependent direction | 4 / 6 |

Six families, 20 displayed genes and 21 explicit pairs are a bounded editorial selection. Six is within the proposed 6–10 range; no extra rows were added merely to meet a larger target. These do not represent six independent pathways.

## Eligible alternatives retained

The other 16 networks spanning three or four species contain thioredoxins; peroxiredoxins; 14-3-3 proteins; AKT; INDY-related transporters; adenine nucleotide transporters; the KSS1/mpk-1/fly MAPK network; alpha-mannosidases; QCR7-related proteins; the RIM11/jnk-1/bsk kinase network; ATG8; SOD1; ATG5; Ku70; COQ7; and catalase. These are candidates, not rejected biological hypotheses or validated uniform families. Their displayed membership would require further checks of paralogy, intervention identity and paper-level outcomes. In particular, kinase connectivity alone must not define an orthogroup.

ATG5 mouse overexpression has a verified lifespan-extension abstract (PMID 23939249), but the retrieved worm autophagy abstract (19469880) and yeast screen abstract (20421943) do not alone establish precise ATG5-specific assay claims. ATG7 supplies a clearer first autophagy row. Thioredoxin/peroxiredoxin alternatives would require additional construct and gene-member checks; SOD2 supplies a fully traced oxidative-stress comparison without needing those substitutions. Unselected entries have not all been given primary-paper review and must not be described as weak or absent evidence.

## Search and evidence review

Workflow: targeted multi-source primary-paper search, using existing GenAge reference mappings to resolve PMIDs, Europe PMC's structured PubMed records for metadata/abstracts, and PubMed/PMC/publisher pages for targeted checks. The academic-search MCP server was unavailable. PubMed web opens sometimes returned empty pages; structured Europe PMC records supplied the abstracts. The Europe PMC full-text endpoint failed for the yeast AMPK paper (HTTP 500); exact AMPK experiment statements were checked in indexed primary PMC full-text passages and the review depth explicitly records this. No inaccessible full text is claimed as downloaded or fully reviewed. Cache collector treats full-text failures as nonfatal and reports them.

New primary-paper claims and their scope are in `expanded_additional_evidence.tsv` and the expanded legend. Existing three-example evidence is reused without changing its original fields. Source annotations and experiment outcomes remain separate, especially worm atg-7 (GenAge Increase/Unclear versus the selected eat-2 longevity-requirement experiment). No raw inputs are modified.

## Verification

Builder checks all 117 text elements for assigned-box containment and pairwise overlap. Independent verification reconstructs candidates using species/Entrez keys, reconciles all 21 selected source rows, all 20 gene records, original evidence preservation, cached new PMIDs, SVG labels, PDF bounds and physical dimensions, raster validity and 600-dpi TIFF. PNG and a Python-rendered PDF preview are visually inspected. Existing complete-analysis and focused-figure verifiers must also pass.
