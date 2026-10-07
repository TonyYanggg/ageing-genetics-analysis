# Task 3 figure contract

## Core conclusion

Show the distribution of earliest dated associated references for current GenAge Build 21 entries, and how cumulative entry counts differ across organisms. These dates do not establish historical first discoveries.

## Evidence chain

1. Extract every GenAge observation-to-primary-reference link for the selected entries.
2. Preserve the full reference table before deriving dates, so genes with several observations or publications remain auditable.
3. Derive an earliest-reference year per GenAge entry for a gene-discovery view.
4. Retain all unique references for a complementary research-activity view.
5. Compare organism-specific coverage and missingness before deciding which view becomes the main panel.

## Figure archetype

Quantitative timeline. The likely hero panel is cumulative genes by earliest associated reference year, stratified by major model organism. Annual unique-reference counts are a possible supporting panel if they add distinct evidence.

## Backend and export contract

- Backend: Python with matplotlib only.
- White background and restrained organism palette.
- Current exploratory exports: PNG and SVG.
- Submission exports, after the design is approved: editable SVG/PDF and 600-dpi TIFF, with a source-data CSV.

## Review risks

- A reference associated with a current GenAge entry is not automatically the historical first discovery of that gene.
- One gene can have several observations and references; one reference can support several genes.
- Live website contents can differ from stable Build 21 downloads.
- Publication year is a coarse proxy for discovery date.
- Missing or dynamically loaded references could bias early-year assignments.
- Organisms differ greatly in record count and research history.

For these reasons, the first output is a validated extraction pilot rather than a publication figure.
