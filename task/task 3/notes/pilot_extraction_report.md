# Task 3 pilot extraction report

Extraction time: `2026-09-17T13:28:24+00:00`  
Requested GenAge IDs: `1, 6, 7, 8, 21, 128, 137, 508, 606, 1792, 2025, 2083`

## Scope

This is a small parser and coverage pilot, not the full timeline dataset. It spans all nine organisms in the Build 21 download and includes several well-known multi-reference entries.

## Results

- Entry pages parsed: **12**
- Cached gene-detail pages reused: **12**
- Cached LibAge reference fragments reused: **25**
- Observation rows: **25**
- Observation-reference rows: **25**
- Unique reference records in the pilot: **22**
- Reference rows without a parsed year: **0**
- Observed publication-year range: **1988–2017**

| Organism | Reference rows | Unique references |
|---|---:|---:|
| Caenorhabditis briggsae | 1 | 1 |
| Caenorhabditis elegans | 12 | 9 |
| Danio rerio | 1 | 1 |
| Drosophila melanogaster | 3 | 3 |
| Mesocricetus auratus | 1 | 1 |
| Mus musculus | 1 | 1 |
| Podospora anserina | 1 | 1 |
| Saccharomyces cerevisiae | 4 | 4 |
| Schizosaccharomyces pombe | 1 | 1 |

## Validation

- Every live page symbol and organism matched the corresponding immutable Build 21 CSV row.
- The `age-1` page yielded the expected reference years 1988, 1990, 1995, 2005, and 2008.
- Raw HTML is cached under `results/pilot/raw_html/`; each derived row records its source URL, extraction time, and page SHA-256.
- Observation rows and reference rows are stored separately, preventing repeated references or multiple observations from being silently collapsed.

## Interpretation rule for the next stage

Extract all references first. Then derive both:

1. the earliest associated reference per GenAge entry, for a candidate gene-discovery timeline; and
2. all unique associated references, for a literature-activity timeline.

The main review figure should be chosen only after comparing coverage, duplicates, and organism-specific distributions. An associated reference year must not be described automatically as the definitive historical discovery year.
