"""Extract a small, auditable pilot of GenAge observation references.

The script intentionally defaults to a fixed pilot rather than the full GenAge
catalogue. HAGR asks users to make contact before issuing large query volumes.
Each downloaded HTML page is cached so validation and parser development do not
repeat requests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

from bs4 import BeautifulSoup, Tag
import pandas as pd
import requests


TASK_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
GENAGE_FILE = PROJECT_ROOT / "GenAge Project" / "models_genes" / "genage_models.csv"
PILOT_ROOT = TASK_ROOT / "results" / "pilot"
RAW_HTML_DIR = PILOT_ROOT / "raw_html"
RAW_REFERENCE_DIR = PILOT_ROOT / "raw_reference_html"
OBSERVATIONS_FILE = PILOT_ROOT / "genage_observations_pilot.csv"
REFERENCES_FILE = PILOT_ROOT / "genage_references_pilot.csv"
REPORT_FILE = TASK_ROOT / "notes" / "pilot_extraction_report.md"

DETAIL_URL = "https://genomics.senescence.info/genes/details.php?id={genage_id}"
REFERENCE_API_URL = "https://libage.ageing-map.org/api/get/genage_model/{source_id}/"
USER_AGENT = "GenAgeTimelineResearch/0.1 (academic, cached, rate-limited pilot)"
PILOT_IDS = [1, 6, 7, 8, 21, 128, 137, 508, 606, 1792, 2025, 2083]
EXPECTED_AGE1_YEARS = {1988, 1990, 1995, 2005, 2008}


def clean_text(value: str) -> str:
    """Collapse HTML whitespace without changing words or punctuation."""

    return " ".join(value.split())


def direct_definition_pairs(block: Tag) -> dict[str, Tag]:
    """Map direct-child dt labels to their following dd elements."""

    pairs: dict[str, Tag] = {}
    for term in block.find_all("dt", recursive=False):
        definition = term.find_next_sibling("dd")
        if definition is not None:
            pairs[clean_text(term.get_text(" ", strip=True))] = definition
    return pairs


def extract_year(citation: str) -> int | None:
    """Extract the first parenthesized four-digit publication year."""

    match = re.search(r"\(((?:18|19|20)\d{2})\)", citation)
    return int(match.group(1)) if match else None


def extract_identifier(href: str | None, pattern: str) -> str:
    """Extract one identifier from a URL or return an empty string."""

    if not href:
        return ""
    match = re.search(pattern, href)
    return match.group(1) if match else ""


def parse_detail_page(
    html: str,
    *,
    requested_id: int,
    source_url: str,
    source_sha256: str,
    extracted_at_utc: str,
) -> tuple[dict[str, str], list[dict[str, object]], list[dict[str, object]]]:
    """Parse gene information, observations, and primary references."""

    soup = BeautifulSoup(html, "html.parser")
    heading = soup.find(lambda tag: tag.name in {"h1", "h2"} and "Gene information" in tag.get_text(" "))
    if heading is None:
        raise ValueError(f"GenAge gene-information heading not found for ID {requested_id}")

    heading_match = re.search(r"ID:\s*(\d+)", heading.get_text(" ", strip=True))
    page_id = int(heading_match.group(1)) if heading_match else requested_id
    if page_id != requested_id:
        raise ValueError(f"Requested GenAge ID {requested_id}, received ID {page_id}")

    info_block = soup.select_one("dl.section-entry:not(.boxed)")
    if info_block is None:
        raise ValueError(f"Gene-information block not found for ID {requested_id}")
    info = direct_definition_pairs(info_block)
    organism_definition = info.get("Organism")
    organism_tag = organism_definition.find("i") if organism_definition else None
    organism = clean_text(
        (organism_tag or organism_definition or Tag(name="dd")).get_text(" ", strip=True)
    )
    gene = {
        "genage_id": str(page_id),
        "symbol": clean_text(info.get("Symbol", Tag(name="dd")).get_text(" ", strip=True)),
        "name": clean_text(info.get("Name", Tag(name="dd")).get_text(" ", strip=True)),
        "organism": organism,
    }

    observation_rows: list[dict[str, object]] = []
    reference_rows: list[dict[str, object]] = []
    for observation_index, block in enumerate(soup.select("dl.section-entry.boxed"), start=1):
        fields = direct_definition_pairs(block)
        reference_definition = fields.get("Primary reference")
        observation_source_id = ""
        if reference_definition:
            reference_list = reference_definition.select_one("ul.references-list")
            observation_source_id = extract_identifier(
                reference_list.get("id") if reference_list else None,
                r"ref-(\d+)",
            )

        observation = {
            **gene,
            "observation_index": observation_index,
            "observation_source_id": observation_source_id,
            "observation": clean_text(fields.get("Observations", Tag(name="dd")).get_text(" ", strip=True)),
            "lifespan_effect": clean_text(fields.get("Lifespan Effect", Tag(name="dd")).get_text(" ", strip=True)),
            "genetic_manipulation": clean_text(fields.get("Genetic Manipulation", Tag(name="dd")).get_text(" ", strip=True)),
            "longevity_category": clean_text(fields.get("Longevity Category", Tag(name="dd")).get_text(" ", strip=True)),
            "reference_count": 0,
            "source_url": source_url,
            "source_sha256": source_sha256,
            "extracted_at_utc": extracted_at_utc,
        }
        observation_rows.append(observation)

    return gene, observation_rows, reference_rows


def parse_reference_fragment(
    html: str,
    *,
    observation: dict[str, object],
    reference_source_url: str,
    reference_source_sha256: str,
    extracted_at_utc: str,
) -> list[dict[str, object]]:
    """Parse the HTML fragment returned by the LibAge observation endpoint."""

    try:
        payload = json.loads(html)
    except json.JSONDecodeError:
        payload = None

    if isinstance(payload, list):
        items: list[tuple[Tag, str]] = []
        for record in payload:
            if not isinstance(record, dict):
                continue
            citation_soup = BeautifulSoup(str(record.get("reference", "")), "html.parser")
            wrapper_soup = BeautifulSoup("<li></li>", "html.parser")
            item = wrapper_soup.li
            if item is None:
                continue
            libage_link = wrapper_soup.new_tag("a", href=str(record.get("url", "")))
            libage_link.string = str(record.get("id", ""))
            item.append(libage_link)
            citation_span = wrapper_soup.new_tag("span", attrs={"class": "libage-reference"})
            for child in list(citation_soup.contents):
                citation_span.append(child)
            item.append(citation_span)
            items.append((item, str(record.get("url", ""))))
    else:
        soup = BeautifulSoup(html, "html.parser")
        items = [(item, "") for item in soup.select("li")]

    reference_rows: list[dict[str, object]] = []
    for reference_index, (item, api_libage_path) in enumerate(items, start=1):
        citation_span = item.select_one("span.libage-reference")
        citation = clean_text((citation_span or item).get_text(" ", strip=True))
        title_tag = citation_span.find("b") if citation_span else None
        pubmed_link = item.select_one("a.libage-pubmed")
        libage_link = item.find(
            "a", href=re.compile(r"libage\.ageing-map\.org/entries/entry/\d+")
        )
        libage_href = urljoin(
            reference_source_url,
            libage_link.get("href", "") if libage_link else api_libage_path,
        )
        pubmed_href = urljoin(reference_source_url, pubmed_link.get("href", "")) if pubmed_link else ""
        reference_rows.append(
            {
                "genage_id": observation["genage_id"],
                "symbol": observation["symbol"],
                "name": observation["name"],
                "organism": observation["organism"],
                "observation_index": observation["observation_index"],
                "observation_source_id": observation["observation_source_id"],
                "reference_index": reference_index,
                "libage_id": extract_identifier(libage_href, r"/entry/(\d+)/?"),
                "pubmed_id": extract_identifier(pubmed_href, r"/pubmed/(\d+)")
                or extract_identifier(pubmed_href, r"/([0-9]+)$"),
                "publication_year": extract_year(citation),
                "title": clean_text(title_tag.get_text(" ", strip=True)) if title_tag else "",
                "citation": citation,
                "libage_url": libage_href,
                "pubmed_url": pubmed_href,
                "detail_source_url": observation["source_url"],
                "detail_source_sha256": observation["source_sha256"],
                "reference_source_url": reference_source_url,
                "reference_source_sha256": reference_source_sha256,
                "extracted_at_utc": extracted_at_utc,
            }
        )
    return reference_rows


def fetch_or_read_cached(
    session: requests.Session,
    genage_id: int,
    *,
    delay_seconds: float,
    refresh: bool,
) -> tuple[str, str, bool]:
    """Read cached HTML or make one polite request and cache the response."""

    RAW_HTML_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = RAW_HTML_DIR / f"genage_{genage_id}.html"
    url = DETAIL_URL.format(genage_id=genage_id)
    if cache_file.exists() and not refresh:
        return cache_file.read_text(encoding="utf-8"), url, True

    response = session.get(url, timeout=30)
    response.raise_for_status()
    response.encoding = response.apparent_encoding or "utf-8"
    html = response.text
    cache_file.write_text(html, encoding="utf-8", newline="\n")
    time.sleep(delay_seconds)
    return html, str(response.url), False


def fetch_or_read_reference_cached(
    session: requests.Session,
    source_id: str,
    *,
    delay_seconds: float,
    refresh: bool,
) -> tuple[str, str, bool]:
    """Read or cache one small LibAge reference HTML fragment."""

    RAW_REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = RAW_REFERENCE_DIR / f"genage_model_{source_id}.html"
    url = REFERENCE_API_URL.format(source_id=source_id)
    if cache_file.exists() and not refresh:
        return cache_file.read_text(encoding="utf-8"), url, True

    response = session.get(url, timeout=30)
    response.raise_for_status()
    response.encoding = response.apparent_encoding or "utf-8"
    html = response.text
    cache_file.write_text(html, encoding="utf-8", newline="\n")
    time.sleep(delay_seconds)
    return html, str(response.url), False


def validate_against_download(gene: dict[str, str], downloaded: pd.DataFrame) -> None:
    """Require the live page identity to match the stable downloaded row."""

    genage_id = int(gene["genage_id"])
    matched = downloaded.loc[downloaded["GenAge ID"].eq(genage_id)]
    if len(matched) != 1:
        raise ValueError(f"Expected one downloaded row for GenAge ID {genage_id}; found {len(matched)}")
    row = matched.iloc[0]
    if gene["symbol"] != str(row["symbol"]):
        raise ValueError(f"Symbol mismatch for GenAge ID {genage_id}: page={gene['symbol']!r}, CSV={row['symbol']!r}")
    if gene["organism"] != str(row["organism"]):
        raise ValueError(f"Organism mismatch for GenAge ID {genage_id}: page={gene['organism']!r}, CSV={row['organism']!r}")


def write_report(
    ids: list[int],
    observations: pd.DataFrame,
    references: pd.DataFrame,
    *,
    extracted_at_utc: str,
    page_cache_hits: int,
    reference_cache_hits: int,
) -> None:
    """Write a compact pilot coverage and validation report."""

    unique_references = references.drop_duplicates(subset=["libage_id", "pubmed_id", "citation"])
    missing_years = int(references["publication_year"].isna().sum())
    age1_years = set(
        references.loc[references["genage_id"].astype(str).eq("7"), "publication_year"]
        .dropna()
        .astype(int)
    )
    missing_expected_age1 = sorted(EXPECTED_AGE1_YEARS - age1_years)
    if missing_expected_age1:
        raise ValueError(f"age-1 pilot validation is missing expected years: {missing_expected_age1}")

    organism_counts = (
        references.groupby("organism", dropna=False)
        .agg(reference_rows=("citation", "size"), unique_references=("citation", "nunique"))
        .reset_index()
    )
    organism_lines = [
        "| Organism | Reference rows | Unique references |",
        "|---|---:|---:|",
    ]
    organism_lines.extend(
        f"| {row.organism} | {row.reference_rows} | {row.unique_references} |"
        for row in organism_counts.itertuples(index=False)
    )
    organism_markdown = "\n".join(organism_lines)
    year_min = int(references["publication_year"].min()) if not references.empty else None
    year_max = int(references["publication_year"].max()) if not references.empty else None

    report = f"""# Task 3 pilot extraction report

Extraction time: `{extracted_at_utc}`  
Requested GenAge IDs: `{', '.join(map(str, ids))}`

## Scope

This is a small parser and coverage pilot, not the full timeline dataset. It spans all nine organisms in the Build 21 download and includes several well-known multi-reference entries.

## Results

- Entry pages parsed: **{len(ids)}**
- Cached gene-detail pages reused: **{page_cache_hits}**
- Cached LibAge reference fragments reused: **{reference_cache_hits}**
- Observation rows: **{len(observations)}**
- Observation-reference rows: **{len(references)}**
- Unique reference records in the pilot: **{len(unique_references)}**
- Reference rows without a parsed year: **{missing_years}**
- Observed publication-year range: **{year_min}–{year_max}**

{organism_markdown}

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
"""
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(report, encoding="utf-8", newline="\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ids",
        nargs="+",
        type=int,
        default=PILOT_IDS,
        help="GenAge IDs to extract. Defaults to the fixed 12-entry pilot.",
    )
    parser.add_argument(
        "--delay-seconds",
        type=float,
        default=2.0,
        help="Delay after each network request; cached pages do not wait.",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Redownload existing cached pages.",
    )
    args = parser.parse_args()
    ids = sorted(set(args.ids))
    if len(ids) > 25:
        raise ValueError("This pilot script is capped at 25 IDs to avoid unapproved bulk querying.")
    if args.delay_seconds < 1.0:
        raise ValueError("Use a delay of at least one second between network requests.")

    downloaded = pd.read_csv(GENAGE_FILE)
    known_ids = set(downloaded["GenAge ID"].astype(int))
    unknown_ids = sorted(set(ids) - known_ids)
    if unknown_ids:
        raise ValueError(f"IDs not present in the Build 21 download: {unknown_ids}")

    extracted_at_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "text/html, */*;q=0.8"})
    all_observations: list[dict[str, object]] = []
    all_references: list[dict[str, object]] = []
    page_cache_hits = 0
    reference_cache_hits = 0

    for genage_id in ids:
        html, source_url, from_cache = fetch_or_read_cached(
            session,
            genage_id,
            delay_seconds=args.delay_seconds,
            refresh=args.refresh,
        )
        page_cache_hits += int(from_cache)
        source_sha256 = hashlib.sha256(html.encode("utf-8")).hexdigest()
        gene, observations, references = parse_detail_page(
            html,
            requested_id=genage_id,
            source_url=source_url,
            source_sha256=source_sha256,
            extracted_at_utc=extracted_at_utc,
        )
        validate_against_download(gene, downloaded)
        if len(all_observations) + len(observations) > 100:
            raise ValueError("Pilot exceeded 100 observations; narrow --ids before querying references.")
        for observation in observations:
            source_id = str(observation["observation_source_id"])
            if not source_id:
                raise ValueError(
                    f"Missing LibAge observation source ID for GenAge ID {genage_id}, "
                    f"observation {observation['observation_index']}"
                )
            reference_html, reference_url, reference_from_cache = fetch_or_read_reference_cached(
                session,
                source_id,
                delay_seconds=args.delay_seconds,
                refresh=args.refresh,
            )
            reference_cache_hits += int(reference_from_cache)
            reference_sha256 = hashlib.sha256(reference_html.encode("utf-8")).hexdigest()
            observation_references = parse_reference_fragment(
                reference_html,
                observation=observation,
                reference_source_url=reference_url,
                reference_source_sha256=reference_sha256,
                extracted_at_utc=extracted_at_utc,
            )
            observation["reference_count"] = len(observation_references)
            references.extend(observation_references)
        all_observations.extend(observations)
        all_references.extend(references)
        print(
            f"GenAge ID {genage_id}: {gene['symbol']} ({gene['organism']}), "
            f"observations={len(observations)}, references={len(references)}, "
            f"source={'cache' if from_cache else 'network'}"
        )

    observations_frame = pd.DataFrame(all_observations)
    references_frame = pd.DataFrame(all_references)
    if observations_frame.empty:
        raise ValueError("No observation rows were extracted")
    if references_frame.empty:
        raise ValueError("No reference rows were extracted")

    PILOT_ROOT.mkdir(parents=True, exist_ok=True)
    observations_frame.to_csv(OBSERVATIONS_FILE, index=False)
    references_frame.to_csv(REFERENCES_FILE, index=False)
    write_report(
        ids,
        observations_frame,
        references_frame,
        extracted_at_utc=extracted_at_utc,
        page_cache_hits=page_cache_hits,
        reference_cache_hits=reference_cache_hits,
    )
    print(f"Wrote {OBSERVATIONS_FILE}")
    print(f"Wrote {REFERENCES_FILE}")
    print(f"Wrote {REPORT_FILE}")


if __name__ == "__main__":
    main()
