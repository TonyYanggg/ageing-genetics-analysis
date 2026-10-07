# Model-organism analysis and figure contract

2026-09-18. Historical implementation contract for the nine-organism source audit, not the current focused review figure. **The common-level grouping is provisional, not an accepted biological result.** A Coq2/HSP104 cross-resource disagreement requires resolution before any group-level claim; see `source_validation.md` and the current Task 4 README.

## Evidence and mapping

Use all 2,205 Build 21 records as the inventory. Qualifying lifespan evidence means a nonempty `lifespan effect` of Increase, Decrease, or Increase and Decrease. This deliberately includes records without Pro/Anti labels and records annotated Necessary for fitness, which remain separately identified. A narrower Pro/Anti sensitivity count accompanies the broad lifespan-evidence count. Records without a lifespan-effect category stay in the inventory and audit, but do not establish the headline recurrence.

Join only explicit OrthoDB NCBI-gene cross-references in the matching organism. Use a single documented representative assembly per organism. S288C yeast, Bristol N2 worm, reference fly, C57BL/6J mouse, 972h- fission yeast, S mat+ Podospora and the selected briggsae/zebrafish/hamster assemblies are recorded in source_species.csv. Genome representation does not guarantee mapping of every GenAge gene or experimental strain.

Retain all returned OrthoDB gene IDs. Multiple IDs linked to one Entrez ID may be collapsed at a taxonomic level only if they assign it to the same group at that level; otherwise mark that mapping ambiguous and exclude it at that level. Never accept the first of multiple search hits or join on symbol alone. Missing IDs and missing cross-references remain unresolved, not negative orthology results.

## Group resolution

The primary common-level summary uses source-defined Eukaryota-level groups (taxon 2759), allowing the same unit across fungi and animals. These are ancient ancestral groups and may encompass multiple paralogous subfamilies; membership is not a claim of one-to-one orthology between every pair. Report model-gene multiplicity and do not name a broad group after a single selected gene.

Retain finer-level groups in a separate table with counts by taxonomic level. Do not sum nested groups into a purported number of independent conserved genes. The broad-versus-fine resolution limitation must be visible in the caption and final interpretation. Source group annotations are descriptions, not independent functional validation.

## Figure contract

Provisional question: which source-defined ancestral groups contain recorded lifespan effects in multiple model organisms, and how do their GenAge influence annotations compare? Numerical conclusions follow the data.

Archetype: quantitative grid. Python is the chosen backend. Panel a: mapping coverage for all nine organisms, with mapped/identifiable-gene denominators and missing-ID counts in source data. Panel b: evidence matrix for a transparent selection of recurrent groups, using source group IDs and model-gene symbols; distinguish Pro, Anti, Mixed, other lifespan evidence, a mapped ortholog without qualifying GenAge evidence, and no recovered group member. Panel c: all common-level groups counted by number of species with lifespan evidence, with direction agreement/discordance/unresolved categories. A complete all-group, all-nine table accompanies the displayed selection.

Selection for the compact matrix is by descending species recurrence, then fewer member genes (to avoid giving the broadest families automatic visual priority), then stable group ID; disclose any separate illustrative subset. Do not select rows from the historical human-anchor result.

Draft review output: 183 mm wide, legible 7 pt text, editable SVG/PDF, PNG preview, 600 dpi TIFF. No significance tests or error bars: genes/groups are catalogue units, not independent experimental replicates. A missing mapped member is unresolved in this analysis, not demonstrated biological absence. Full orthology resolution and selected experimental evidence need review before strong conservation claims.

## Checks

Freeze REST release and RDF source version, public query text, raw responses, retrieval times and checksums. Reject partial/time-limited query responses. Reconcile queried IDs with all input IDs. Audit taxonomy mismatches, ambiguous group assignments, missing IDs, multi-ID mappings and coverage. Verify every derived recurrence count from cached raw memberships independently, maintain raw-file integrity, and inspect figure exports.
