# Ageing genetics: cross-species evidence and human disease associations

Reproducible computational analyses supporting a review of the genetics of ageing. The project integrates GenAge/HAGR catalogues, publication records, explicit orthology relationships and human gene–disease associations to produce traceable figures and research tables.

**Research and analysis: Zhou Yang.** Prepared in support of João Pedro de Magalhães's ageing-genetics review. This repository is a research archive; review-ready outputs are not a claim of publication or approval of every figure.

## Results at a glance

| Analysis | Coverage | Main deliverable |
|---|---|---|
| [GenAge inventory](task/task%201/README.md) | 2,205 model-organism records across nine organisms | Descriptive counts and lifespan-annotation figures |
| [Lifespan-extension comparison](task/task%202/README.md) | 830 qualifying records across four organisms | Verified effect-magnitude figure |
| [Reference timeline](task/task%203/README.md) | 2,203/2,205 entries dated; 708 mapped references | Overall and organism-specific timelines |
| [Cross-species orthology](task/task%204/README.md) | 320 genes; 206 direct pairs across four organisms | 36-page atlas, searchable index and experimental examples |
| [Human gene associations](task/task%205/README.md) | All 307 GenAge human genes retained; 306 queryable targets | 422,250 unique gene–term relations and evidence scores |

The first two workstreams are complete. Timeline figures, the expanded orthology atlas and human disease-association package are ready for review. The earlier three-example orthology illustration received positive feedback; this does not imply approval of the expanded atlas. Functional enrichment was considered but not pursued.

## Explore the figures

### Lifespan extension across model organisms

![Reported lifespan-extension magnitudes across four model organisms](task/task%202/figures/exploratory/fig6_lifespan_extension_magnitude.png)

Individual source records are shown alongside distribution summaries. The analysis retains source-reported extremes and uses explicit inclusion rules. [Caption and interpretation](task/task%202/notes/fig6_caption.md).

### Orthologs and experimental context

![Selected orthologous genes and reported lifespan interventions](task/task%204/figures/expanded/expanded_ortholog_examples.png)

The selected TOR, S6 kinase, insulin/IGF receptor, AMPK, ATG7 and SOD2 examples retain intervention-specific qualifications. [Complete atlas PDF](task/task%204/figures/all_gene_atlas/all_gene_ortholog_atlas.pdf) · [Atlas legend](task/task%204/notes/all_gene_atlas_legend.md) · [Complete pair table](task/task%204/results/complete_four_species/all_qualifying_pairs.csv).

After cloning, open `task/task 4/figures/all_gene_atlas/index.html` locally to search the atlas by gene. GitHub displays HTML source rather than running this local index.

### Additional review outputs

- [Overall reference timeline](task/task%203/figures/review/timeline_review.png) and [organism-specific timelines](task/task%203/figures/review/timeline_by_organism_review.png).
- [Human gene–disease association figure](task/task%205/figures/task5_human_gene_disease_review.png) and [caption](task/task%205/notes/figure_1_legend.md).
- [Task index and detailed methods](task/README.md).

## Reproduce and verify

The publication copy was prepared on 7 October 2026. The original analyses use frozen source releases; retrieving a newer database release is a separate analysis.

Install Git LFS before cloning so large tables and exports are retrieved as actual files:

```sh
git lfs install
git clone https://github.com/TonyYanggg/ageing-genetics-analysis.git
cd ageing-genetics-analysis
git lfs pull
python -m venv .venv
```

Activate the environment with `.venv\Scripts\Activate.ps1` on Windows PowerShell, or `source .venv/bin/activate` on macOS/Linux, then install the recorded dependencies:

```sh
python -m pip install -r requirements.txt
```

Python 3.13 on Windows is the recorded verification environment. Cross-platform execution has not been tested; some historical manifests preserve Windows-relative paths. UTF-8 mode is recommended. The full archive is approximately 1 GB before Git metadata; nine files totaling approximately 353 MB use LFS. Allow additional space for the checkout and LFS object cache.

Run the saved verification workflows from the repository root, without refreshing external databases:

```sh
python -X utf8 "task/task 2/scripts/verify_fig6_data.py"
python -X utf8 "task/task 4/scripts/verify_complete_four_species.py"
python -X utf8 "task/task 4/scripts/verify_focused_examples.py"
python -X utf8 "task/task 4/scripts/verify_expanded_examples.py"
python -X utf8 "task/task 4/scripts/verify_all_gene_atlas.py"
python -X utf8 "task/task 5/scripts/verify_opentargets_outputs.py"
```

Verifiers may refresh their task-local reports. Task READMEs provide figure-generation commands and distinguish optional network collectors from analysis using cached inputs. See [publication verification](PUBLICATION_VERIFICATION.md) for the upload checks.

## Evidence boundaries

- GenAge records are not always unique biological genes. Counts retain the documented unit of analysis.
- Timeline dates are earliest associated references captured by the source, not definitive historical discovery dates.
- Direct ortholog pairs do not establish conserved causal longevity mechanisms. Missing mappings are not biological negatives.
- The atlas distinguishes 20 reviewed experimental cards from 300 catalogue-annotation cards; it does not imply full-paper review for all 320 genes.
- Gene–term relations include diseases, phenotypes and other classes. Disease association does not establish a human-lifespan effect; database scores are not probabilities or causal effect sizes.

## Archive structure and provenance

`GenAge Project/` and `Other downloads from HAGR/` contain unchanged source material. `task/task 1` through `task/task 5` retain scripts, notes, cached sources, derived tables, figures and verification records. Historical alternatives are retained and labelled in the task READMEs. No Task 6 analysis is included.

The public copy omits temporary Python caches and internal working instructions. Documentation has been adapted for readers, while research code, source datasets and figure files were preserved during preparation. The original pre-publication archive had no Git history; dated references to that state in historical reports describe the analysis period, not the current repository.

## Sources and attribution

The source databases and publications remain the work of their respective authors. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for source releases, links, attribution and reuse terms. No blanket licence is imposed on mixed third-party materials, and this archive does not imply endorsement by database providers.
