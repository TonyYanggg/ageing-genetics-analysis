# Orthology source audit — 18 September 2026

**Update:** the [pairwise replacement](pairwise_resolution.md) addresses the operational grouping problem by requiring explicit ortholog calls. This note preserves the original reason for rejecting the broad-group interpretation. It is not a claim that every source conflict or coverage gap is biologically resolved.

## Decision

The all-nine-organism pipeline is implemented, but its OrthoDB Eukaryota-level groups are **source-audit outputs only**. Do not describe the 181 recurrent source groups as validated longevity orthogroups or conserved longevity genes. No final biological group count has been accepted.

## Concrete concern

OrthoDB 12.2 group `3459473at2759` includes fly Coq2 (Entrez 40988; OrthoDB `7227_0:000d2a`) and yeast HSP104 (YLL026W), alongside other members. REST queries reproduce the SPARQL assignment, so this is not merely a CSV join discrepancy. Coq2 is a coenzyme-Q biosynthesis protein; the broad group's label is chaperone.

An independent OMA lookup assigns fly `DROME19540` / `COQ2_DROME` to root HOG 848365 and yeast `P31539` / `HS104_YEAST` to root HOG 848368. This is a disagreement in group resolution/assignment, not proof that either database is universally correct. Removing one suspicious row would not validate the remaining groups.

Reproduce with `python "task/task 4/scripts/audit_orthology_source.py"`. Raw JSON responses and request metadata/checksums are retained in `results/model_orthology/source_cache/`; the summary is `source_validation.json`. OMA checks use live API responses frozen at retrieval, not an asserted release number.

## Next scientific step

Assess an independent model-to-model orthology resource and its coverage across all nine organisms, then choose defensible group resolution or explicit pairwise relationships. Resolve one-to-many relations and cross-source disagreements before accepting recurrent groups. Preserve the current audit to make the source-selection decision traceable. Experimental-context review follows, not precedes, orthology validation.

The coverage panel measures matching to this particular source and level; it is not biological conservation. Missing exact identifiers remain unresolved. In particular, zero mapped C. briggsae records does not mean absence of a longevity ortholog. Non-Pro/Anti lifespan evidence remains eligible for the broad screen. The narrower Pro/Anti sensitivity analysis is still pending.
