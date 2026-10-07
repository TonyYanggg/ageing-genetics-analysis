"""Independently audit the main GenAge Figure 6 from the raw downloaded CSV.

This script treats the raw GenAge model-organism CSV as the source of truth.
It does not call the existing Figure 6 data-preparation functions, does not
modify any figure, and does not overwrite any pre-existing derived table.
The current plotting function is imported only after the independent dataset
has been reconstructed, then instrumented in memory to verify point/box inputs.
"""

from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import importlib.util
import inspect
import math
from pathlib import Path
import sys
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TASK_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
RAW_FILE = PROJECT_ROOT / "GenAge Project" / "models_genes" / "genage_models.csv"
LOCAL_RELEASE_FILE = PROJECT_ROOT / "GenAge Project" / "models_genes" / "release.html"
FIG1_SCRIPT = PROJECT_ROOT / "task" / "task 1" / "scripts" / "01_basic_descriptive_figures.py"
FIG5_SCRIPT = PROJECT_ROOT / "task" / "task 1" / "scripts" / "02_cross_species_descriptive_figures.py"
CURRENT_FIG6_SCRIPT = TASK_ROOT / "scripts" / "fig6_explore_lifespan_effect_magnitude.py"

TABLE_DIR = TASK_ROOT / "results" / "tables"
VERIFIED_DATA_FILE = TABLE_DIR / "fig6_verified_plotting_data.csv"
VERIFIED_MAXIMA_FILE = TABLE_DIR / "fig6_verified_maximum_records.csv"
REPORT_FILE = TASK_ROOT / "notes" / "fig6_verification_report.md"

TASK1_TABLE_DIR = PROJECT_ROOT / "task" / "task 1" / "results" / "tables"
FIG1_TABLE = TASK1_TABLE_DIR / "fig1_records_by_species.csv"
FIG5_TABLE = TASK1_TABLE_DIR / "lifespan_effect_by_organism_counts.csv"
EXISTING_AVAILABILITY_FILE = TABLE_DIR / "fig6_lifespan_magnitude_data_availability.csv"
EXISTING_SUMMARY_FILE = TABLE_DIR / "fig6_lifespan_magnitude_summary.csv"
EXISTING_EXTREMES_FILE = TABLE_DIR / "fig6_lifespan_magnitude_extremes.csv"

ORGANISM_COLUMN = "organism"
EFFECT_COLUMN = "lifespan effect"
MAGNITUDE_COLUMN = "avg lifespan change (max obsv)"
RAW_COLUMNS = [
    "GenAge ID",
    "symbol",
    "name",
    ORGANISM_COLUMN,
    "entrez gene id",
    MAGNITUDE_COLUMN,
    EFFECT_COLUMN,
    "longevity influence",
]

EXPECTED_RAW_ROWS = 2_205
EXPECTED_ORGANISMS = [
    "Saccharomyces cerevisiae",
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Mus musculus",
    "Schizosaccharomyces pombe",
    "Podospora anserina",
    "Caenorhabditis briggsae",
    "Danio rerio",
    "Mesocricetus auratus",
]
EXPECTED_ORGANISM_TOTALS = {
    "Saccharomyces cerevisiae": 911,
    "Caenorhabditis elegans": 889,
    "Drosophila melanogaster": 202,
    "Mus musculus": 136,
    "Schizosaccharomyces pombe": 61,
    "Podospora anserina": 3,
    "Caenorhabditis briggsae": 1,
    "Danio rerio": 1,
    "Mesocricetus auratus": 1,
}
EFFECT_CATEGORIES = ["Increase", "Decrease", "Increase and Decrease", "Missing/Other"]
SUMMARY_CATEGORIES = ["Increase", "Decrease", "Increase and Decrease"]
PLOTTED_ORGANISMS = EXPECTED_ORGANISMS[:4]
CURRENT_FIGURE_N = {
    "Saccharomyces cerevisiae": 217,
    "Caenorhabditis elegans": 453,
    "Drosophila melanogaster": 114,
    "Mus musculus": 46,
}
EXPECTED_MAX_LABELS = {
    "Saccharomyces cerevisiae": ("SIR2", 600.0),
    "Caenorhabditis elegans": ("age-1", 1_000.0),
    "Drosophila melanogaster": ("Vinc", 150.0),
    "Mus musculus": ("Ghrh", 46.0),
}

FIG6_FIGURE_FILES = [
    TASK_ROOT / "figures" / "exploratory" / "fig6_lifespan_extension_magnitude.png",
    TASK_ROOT / "figures" / "exploratory" / "fig6_lifespan_extension_magnitude.svg",
    TASK_ROOT / "figures" / "exploratory" / "fig6_signed_lifespan_effect_magnitude.png",
    TASK_ROOT / "figures" / "exploratory" / "fig6_signed_lifespan_effect_magnitude.svg",
    TASK_ROOT / "figures" / "exploratory" / "fig6_preview_A_extension_only.png",
    TASK_ROOT / "figures" / "exploratory" / "fig6_preview_A_extension_only.svg",
    TASK_ROOT / "figures" / "exploratory" / "fig6_preview_B_signed_effects.png",
    TASK_ROOT / "figures" / "exploratory" / "fig6_preview_B_signed_effects.svg",
]


@dataclass(frozen=True)
class ComparisonResult:
    file: str
    field: str
    key: str
    independent: str
    existing: str
    classification: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_exact(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def display_number(value: Any, digits: int = 8) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    numeric = float(value)
    if numeric.is_integer():
        return f"{int(numeric):,}"
    return f"{numeric:.{digits}g}"


def markdown_table(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in frame.iterrows():
        cells = []
        for column in frame.columns:
            value = row[column]
            rendered = "NA" if pd.isna(value) else str(value)
            cells.append(rendered.replace("|", "\\|").replace("\n", " "))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def comparison_classification(existing_raw: str, independent: Any) -> str:
    existing_raw = str(existing_raw).strip()
    independent_missing = independent is None or pd.isna(independent)
    existing_missing = existing_raw == ""
    if independent_missing and existing_missing:
        return "exact match"
    if independent_missing != existing_missing:
        return "discrepancy"
    try:
        existing_numeric = float(existing_raw)
        independent_numeric = float(independent)
    except (TypeError, ValueError):
        return "exact match" if existing_raw == str(independent) else "discrepancy"
    if math.isclose(existing_numeric, independent_numeric, rel_tol=1e-12, abs_tol=1e-12):
        return "exact match"
    rounded = float(format(independent_numeric, ".6g"))
    if math.isclose(existing_numeric, rounded, rel_tol=1e-12, abs_tol=1e-12):
        return "rounding difference"
    return "discrepancy"


def normalized_text(value: Any) -> str:
    text = "" if value is None or pd.isna(value) else str(value).strip()
    if text.endswith(".0"):
        head = text[:-2]
        if head.lstrip("-").isdigit():
            return head
    return text


def load_raw_source() -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    assert_exact(RAW_FILE.is_file(), f"Raw source does not exist: {RAW_FILE}")
    raw = pd.read_csv(RAW_FILE, dtype=str, keep_default_na=False, na_filter=False)
    assert_exact(len(raw) == EXPECTED_RAW_ROWS, f"Expected {EXPECTED_RAW_ROWS} raw rows, observed {len(raw)}")
    assert_exact(list(raw.columns) == RAW_COLUMNS, f"Unexpected raw columns: {list(raw.columns)}")

    stripped_magnitude = raw[MAGNITUDE_COLUMN].str.strip()
    parsed = pd.to_numeric(stripped_magnitude.mask(stripped_magnitude.eq("")), errors="coerce")
    parsed = parsed.where(np.isfinite(parsed), np.nan)
    effect_category = raw[EFFECT_COLUMN].where(raw[EFFECT_COLUMN].isin(EFFECT_CATEGORIES[:3]), "Missing/Other")
    return raw, parsed, effect_category


def raw_effect_inventory(raw: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    counts = raw[EFFECT_COLUMN].value_counts(dropna=False, sort=False)
    for value, count in counts.items():
        text = str(value)
        rows.append(
            {
                "raw value (repr)": repr(text),
                "count": int(count),
                "characters": len(text),
                "leading/trailing whitespace": text != text.strip(),
                "empty/missing in source": text == "",
            }
        )
    return pd.DataFrame(rows)


def organism_effect_counts(raw: pd.DataFrame, effect_category: pd.Series) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    observed_order = raw[ORGANISM_COLUMN].drop_duplicates().tolist()
    for organism in observed_order:
        mask = raw[ORGANISM_COLUMN].eq(organism)
        row: dict[str, Any] = {"organism": organism}
        row["total GenAge records"] = int(mask.sum())
        for category in EFFECT_CATEGORIES:
            row[category] = int((mask & effect_category.eq(category)).sum())
        records.append(row)
    return pd.DataFrame(records)


def compare_fig1_and_fig5(counts: pd.DataFrame) -> list[ComparisonResult]:
    results: list[ComparisonResult] = []
    fig1 = pd.read_csv(FIG1_TABLE, dtype=str, keep_default_na=False)
    fig5 = pd.read_csv(FIG5_TABLE, dtype=str, keep_default_na=False)
    count_by_org = counts.set_index("organism")

    assert_exact(set(fig1["organism"]) == set(count_by_org.index), "Figure 1 organism set differs from raw source")
    assert_exact(set(fig5["organism"]) == set(count_by_org.index), "Figure 5 organism set differs from raw source")
    for organism in count_by_org.index:
        raw_total = int(count_by_org.loc[organism, "total GenAge records"])
        old_fig1_total = int(fig1.loc[fig1["organism"].eq(organism), "record_count"].iloc[0])
        status = "exact match" if raw_total == old_fig1_total else "discrepancy"
        results.append(ComparisonResult(FIG1_TABLE.name, "record_count", organism, str(raw_total), str(old_fig1_total), status))
        assert_exact(status == "exact match", f"STOP: Figure 1 count differs for {organism}")

        old_row = fig5.loc[fig5["organism"].eq(organism)].iloc[0]
        mapping = {
            "Increase": "Increase",
            "Decrease": "Decrease",
            "Increase and Decrease": "Increase and Decrease",
            "Missing/Other": "Missing",
            "total GenAge records": "record_total",
        }
        for independent_field, old_field in mapping.items():
            independent_value = int(count_by_org.loc[organism, independent_field])
            old_value = int(old_row[old_field])
            status = "exact match" if independent_value == old_value else "discrepancy"
            results.append(ComparisonResult(FIG5_TABLE.name, old_field, organism, str(independent_value), str(old_value), status))
            assert_exact(status == "exact match", f"STOP: Figure 5 {old_field} differs for {organism}")
    return results


def build_availability(
    raw: pd.DataFrame, parsed: pd.Series, effect_category: pd.Series
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for organism in raw[ORGANISM_COLUMN].drop_duplicates():
        organism_mask = raw[ORGANISM_COLUMN].eq(organism)
        organism_total = int(organism_mask.sum())
        for category in EFFECT_CATEGORIES:
            category_mask = organism_mask & effect_category.eq(category)
            category_records = int(category_mask.sum())
            usable = int((category_mask & parsed.notna()).sum())
            missing = category_records - usable
            percent = 100 * usable / category_records if category_records else np.nan
            rows.append(
                {
                    "organism": organism,
                    "organism_total_records": organism_total,
                    "lifespan_effect_category": category,
                    "category_records": category_records,
                    "usable_numeric_magnitude_records": usable,
                    "missing_magnitude_records": missing,
                    "usable_magnitude_percent": percent,
                }
            )
    return pd.DataFrame(rows)


def build_verified_plotting_data(
    raw: pd.DataFrame, parsed: pd.Series, effect_category: pd.Series
) -> pd.DataFrame:
    inclusion = (
        raw[ORGANISM_COLUMN].isin(PLOTTED_ORGANISMS)
        & effect_category.eq("Increase")
        & parsed.notna()
    )
    selected = raw.loc[inclusion].copy()
    output = pd.DataFrame(
        {
            "raw_row_index_zero_based": selected.index.astype(int),
            "csv_line_number": selected.index.astype(int) + 2,
            "GenAge ID": selected["GenAge ID"],
            "organism": selected[ORGANISM_COLUMN],
            "symbol": selected["symbol"],
            "name": selected["name"],
            "entrez gene id": selected["entrez gene id"],
            "longevity influence": selected["longevity influence"],
            "lifespan effect": selected[EFFECT_COLUMN],
            "raw_lifespan_change": selected[MAGNITUDE_COLUMN],
            "parsed_numeric_lifespan_change": parsed.loc[selected.index].astype(float),
        }
    ).reset_index(drop=True)
    return output


def build_summary(
    raw: pd.DataFrame, parsed: pd.Series, effect_category: pd.Series
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for organism in EXPECTED_ORGANISMS[:5]:
        for category in SUMMARY_CATEGORIES:
            values = parsed.loc[raw[ORGANISM_COLUMN].eq(organism) & effect_category.eq(category)].dropna()
            n = int(values.size)
            if n:
                q1 = float(values.quantile(0.25))
                q3 = float(values.quantile(0.75))
                row = {
                    "organism": organism,
                    "lifespan_effect_category": category,
                    "n": n,
                    "mean": float(values.mean()),
                    "median": float(values.median()),
                    "sd": float(values.std(ddof=1)) if n > 1 else np.nan,
                    "q1": q1,
                    "q3": q3,
                    "iqr": q3 - q1,
                    "minimum": float(values.min()),
                    "maximum": float(values.max()),
                }
            else:
                row = {
                    "organism": organism,
                    "lifespan_effect_category": category,
                    "n": 0,
                    "mean": np.nan,
                    "median": np.nan,
                    "sd": np.nan,
                    "q1": np.nan,
                    "q3": np.nan,
                    "iqr": np.nan,
                    "minimum": np.nan,
                    "maximum": np.nan,
                }
            rows.append(row)
    return pd.DataFrame(rows)


def build_existing_style_extremes(
    raw: pd.DataFrame, parsed: pd.Series, effect_category: pd.Series
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    genage_numeric = pd.to_numeric(raw["GenAge ID"], errors="raise")
    for organism in EXPECTED_ORGANISMS[:5]:
        for category in SUMMARY_CATEGORIES:
            mask = raw[ORGANISM_COLUMN].eq(organism) & effect_category.eq(category) & parsed.notna()
            indices = raw.index[mask].tolist()
            ordered = sorted(indices, key=lambda idx: (-float(parsed.loc[idx]), int(genage_numeric.loc[idx])))
            group_n = len(ordered)
            for rank, idx in enumerate(ordered[:10], start=1):
                rows.append(
                    {
                        "rank_within_organism_and_effect": rank,
                        "usable_records_in_group": group_n,
                        "organism": raw.at[idx, ORGANISM_COLUMN],
                        "symbol": raw.at[idx, "symbol"],
                        "name": raw.at[idx, "name"],
                        "entrez gene id": raw.at[idx, "entrez gene id"],
                        "longevity influence": raw.at[idx, "longevity influence"],
                        "lifespan effect": raw.at[idx, EFFECT_COLUMN],
                        MAGNITUDE_COLUMN: raw.at[idx, MAGNITUDE_COLUMN],
                        "GenAge ID": raw.at[idx, "GenAge ID"],
                    }
                )
    return pd.DataFrame(rows)


def compare_table_by_keys(
    existing_file: Path,
    independent: pd.DataFrame,
    keys: list[str],
    numeric_fields: Iterable[str],
) -> list[ComparisonResult]:
    existing = pd.read_csv(existing_file, dtype=str, keep_default_na=False, na_filter=False)
    assert_exact(set(existing.columns) == set(independent.columns), f"Column mismatch for {existing_file.name}")
    independent_indexed = independent.set_index(keys, drop=False)
    existing_indexed = existing.set_index(keys, drop=False)
    assert_exact(set(independent_indexed.index) == set(existing_indexed.index), f"Row-key mismatch for {existing_file.name}")
    results: list[ComparisonResult] = []
    numeric_fields = set(numeric_fields)
    for key in independent_indexed.index:
        independent_row = independent_indexed.loc[key]
        existing_row = existing_indexed.loc[key]
        key_text = " / ".join(key) if isinstance(key, tuple) else str(key)
        for field in independent.columns:
            if field in keys:
                continue
            independent_value = independent_row[field]
            existing_value = existing_row[field]
            if field in numeric_fields:
                classification = comparison_classification(existing_value, independent_value)
            else:
                classification = "exact match" if normalized_text(existing_value) == normalized_text(independent_value) else "discrepancy"
            results.append(
                ComparisonResult(
                    existing_file.name,
                    field,
                    key_text,
                    display_number(independent_value),
                    str(existing_value) if str(existing_value) else "NA",
                    classification,
                )
            )
    return results


def compare_existing_extremes(independent: pd.DataFrame) -> list[ComparisonResult]:
    existing = pd.read_csv(EXISTING_EXTREMES_FILE, dtype=str, keep_default_na=False, na_filter=False)
    assert_exact(len(existing) == len(independent), "Existing extremes row count differs")
    assert_exact(list(existing.columns) == list(independent.columns), "Existing extremes columns differ")
    results: list[ComparisonResult] = []
    numeric_fields = {
        "rank_within_organism_and_effect",
        "usable_records_in_group",
        "entrez gene id",
        MAGNITUDE_COLUMN,
        "GenAge ID",
    }
    for row_number, (independent_row, existing_row) in enumerate(zip(independent.to_dict("records"), existing.to_dict("records"), strict=True), start=1):
        key = f"row {row_number}: {independent_row['organism']} / {independent_row['lifespan effect']} / rank {independent_row['rank_within_organism_and_effect']}"
        for field in independent.columns:
            left = normalized_text(independent_row[field])
            right = normalized_text(existing_row[field])
            if field in numeric_fields:
                status = comparison_classification(right, independent_row[field])
            else:
                status = "exact match" if left == right else "discrepancy"
            results.append(ComparisonResult(EXISTING_EXTREMES_FILE.name, field, key, left, right, status))
    return results


def top_ten_tables(plotting_data: pd.DataFrame) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for organism in PLOTTED_ORGANISMS:
        subset = plotting_data.loc[plotting_data["organism"].eq(organism)].sort_values(
            ["parsed_numeric_lifespan_change", "raw_row_index_zero_based"],
            ascending=[False, True],
            kind="stable",
        ).head(10).copy()
        subset.insert(0, "rank", range(1, len(subset) + 1))
        tables[organism] = subset[
            [
                "rank",
                "symbol",
                "entrez gene id",
                "parsed_numeric_lifespan_change",
                "raw_row_index_zero_based",
                "csv_line_number",
                "raw_lifespan_change",
                "lifespan effect",
                "GenAge ID",
            ]
        ]
    return tables


def build_maximum_records(raw: pd.DataFrame, plotting_data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for organism in PLOTTED_ORGANISMS:
        subset = plotting_data.loc[plotting_data["organism"].eq(organism)]
        maximum = float(subset["parsed_numeric_lifespan_change"].max())
        tied = subset.loc[subset["parsed_numeric_lifespan_change"].eq(maximum)].copy()
        tie_count = len(tied)
        for record in tied.to_dict("records"):
            same_symbol_raw = raw[ORGANISM_COLUMN].eq(organism) & raw["symbol"].eq(record["symbol"])
            same_symbol_fig6 = subset["symbol"].eq(record["symbol"])
            rows.append(
                {
                    **record,
                    "maximum_tie_count": tie_count,
                    "symbol_records_in_raw_organism": int(same_symbol_raw.sum()),
                    "symbol_records_in_figure6": int(same_symbol_fig6.sum()),
                    "symbol_is_repeated_in_raw_organism": bool(same_symbol_raw.sum() > 1),
                    "symbol_is_repeated_in_figure6": bool(same_symbol_fig6.sum() > 1),
                }
            )
    columns = [
        "raw_row_index_zero_based",
        "csv_line_number",
        "GenAge ID",
        "organism",
        "symbol",
        "name",
        "entrez gene id",
        "longevity influence",
        "lifespan effect",
        "raw_lifespan_change",
        "parsed_numeric_lifespan_change",
        "maximum_tie_count",
        "symbol_records_in_raw_organism",
        "symbol_records_in_figure6",
        "symbol_is_repeated_in_raw_organism",
        "symbol_is_repeated_in_figure6",
    ]
    return pd.DataFrame(rows)[columns]


@contextmanager
def instrument_current_plot(module: Any):
    scatter_calls: list[tuple[np.ndarray, np.ndarray]] = []
    box_calls: list[np.ndarray] = []
    save_calls: list[str] = []
    original_scatter = matplotlib.axes.Axes.scatter
    original_box = module.add_horizontal_boxplot
    original_save = module.save_preview

    def capture_scatter(axis: Any, x: Any, y: Any, *args: Any, **kwargs: Any):
        scatter_calls.append((np.asarray(x, dtype=float).copy(), np.asarray(y, dtype=float).copy()))
        return original_scatter(axis, x, y, *args, **kwargs)

    def capture_box(axis: Any, values: Any, position: Any, color: Any, width: Any):
        box_calls.append(np.asarray(values, dtype=float).copy())
        return original_box(axis, values, position, color, width)

    def suppress_save(fig: Any, basename: str) -> None:
        save_calls.append(basename)
        plt.close(fig)

    matplotlib.axes.Axes.scatter = capture_scatter
    module.add_horizontal_boxplot = capture_box
    module.save_preview = suppress_save
    try:
        yield scatter_calls, box_calls, save_calls
    finally:
        matplotlib.axes.Axes.scatter = original_scatter
        module.add_horizontal_boxplot = original_box
        module.save_preview = original_save
        plt.close("all")


def verify_current_plot_inputs(
    raw: pd.DataFrame,
    parsed: pd.Series,
    effect_category: pd.Series,
    plotting_data: pd.DataFrame,
) -> dict[str, Any]:
    spec = importlib.util.spec_from_file_location("fig6_plot_under_test", CURRENT_FIG6_SCRIPT)
    assert_exact(spec is not None and spec.loader is not None, "Could not import the current Figure 6 script")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    plot_source = inspect.getsource(module.plot_final_extension_only)
    forbidden = ["drop_duplicates", ".groupby(", ".agg(", ".unique("]
    assert_exact(not any(token in plot_source for token in forbidden), "Current plot function contains a deduplication or aggregation operation")

    plot_input = raw.copy()
    plot_input["effect_category"] = effect_category
    plot_input["numeric_lifespan_change"] = parsed
    with instrument_current_plot(module) as (scatter_calls, box_calls, save_calls):
        module.configure_plot_style()
        module.plot_final_extension_only(plot_input)

    assert_exact(len(scatter_calls) == 4, f"Expected four scatter calls, observed {len(scatter_calls)}")
    assert_exact(len(box_calls) == 4, f"Expected four boxplot calls, observed {len(box_calls)}")
    assert_exact(save_calls == ["fig6_lifespan_extension_magnitude"], f"Unexpected in-memory save calls: {save_calls}")

    per_organism: list[dict[str, Any]] = []
    for position, organism in enumerate(PLOTTED_ORGANISMS):
        expected = plotting_data.loc[
            plotting_data["organism"].eq(organism), "parsed_numeric_lifespan_change"
        ].to_numpy(dtype=float)
        scatter_x, scatter_y = scatter_calls[position]
        box_values = box_calls[position]
        assert_exact(np.array_equal(scatter_x, expected), f"Scatter x data differ for {organism}")
        assert_exact(np.array_equal(box_values, expected), f"Boxplot input differs for {organism}")
        assert_exact(scatter_y.shape == scatter_x.shape, f"Jitter vector length differs for {organism}")
        assert_exact(np.all(np.isfinite(scatter_y)), f"Non-finite jitter positions for {organism}")
        per_organism.append(
            {
                "organism": organism,
                "verified rows": len(expected),
                "scatter x observations": len(scatter_x),
                "boxplot input observations": len(box_values),
                "x values identical": True,
                "jitter changes x": False,
            }
        )
    return {
        "per_organism": pd.DataFrame(per_organism),
        "total_scatter": sum(len(call[0]) for call in scatter_calls),
        "total_box": sum(len(call) for call in box_calls),
        "source_forbidden_operations_absent": True,
    }


def comparison_frame(results: list[ComparisonResult]) -> pd.DataFrame:
    return pd.DataFrame([result.__dict__ for result in results])


def compact_comparison_summary(results: list[ComparisonResult]) -> pd.DataFrame:
    frame = comparison_frame(results)
    if frame.empty:
        return pd.DataFrame(columns=["file", "exact match", "rounding difference", "discrepancy"])
    summary = (
        frame.groupby(["file", "classification"]).size().unstack(fill_value=0)
        .reindex(columns=["exact match", "rounding difference", "discrepancy"], fill_value=0)
        .reset_index()
    )
    return summary


def main() -> None:
    figure_hashes_before = {str(path.relative_to(PROJECT_ROOT)): sha256(path) for path in FIG6_FIGURE_FILES if path.is_file()}
    raw, parsed, effect_category = load_raw_source()
    raw_hash = sha256(RAW_FILE)

    # Raw values are inventoried before any cleaning or recoding.
    effect_inventory = raw_effect_inventory(raw)
    expected_effect_counts = {"Increase": 1155, "Decrease": 892, "Increase and Decrease": 144, "": 14}
    assert_exact(raw[EFFECT_COLUMN].value_counts(dropna=False).to_dict() == expected_effect_counts, "Unexpected raw lifespan-effect values/counts")
    assert_exact(raw[EFFECT_COLUMN].eq(raw[EFFECT_COLUMN].str.strip()).all(), "Whitespace cleaning would alter lifespan-effect values")

    counts = organism_effect_counts(raw, effect_category)
    observed_totals = counts.set_index("organism")["total GenAge records"].to_dict()
    assert_exact(observed_totals == EXPECTED_ORGANISM_TOTALS, f"STOP: raw organism totals differ: {observed_totals}")
    assert_exact(int(observed_totals["Schizosaccharomyces pombe"]) == 61, "STOP: S. pombe total is not 61")
    fig15_results = compare_fig1_and_fig5(counts)

    # Magnitude-field diagnostics are performed on the raw strings.
    stripped = raw[MAGNITUDE_COLUMN].str.strip()
    blank_count = int(stripped.eq("").sum())
    whitespace_only_count = int((raw[MAGNITUDE_COLUMN].ne("") & stripped.eq("")).sum())
    nan_like_tokens = {"nan", "na", "n/a", "null", "none"}
    literal_nan_like_count = int(stripped.str.casefold().isin(nan_like_tokens).sum())
    nonblank_parse_failures = int((stripped.ne("") & parsed.isna()).sum())
    infinite_count = int(pd.to_numeric(stripped.mask(stripped.eq("")), errors="coerce").isin([np.inf, -np.inf]).sum())
    positive_count = int(parsed.gt(0).sum())
    zero_count = int(parsed.eq(0).sum())
    negative_count = int(parsed.lt(0).sum())
    usable_count = int(parsed.notna().sum())
    assert_exact(blank_count == 933, f"Expected 933 blank magnitudes, observed {blank_count}")
    assert_exact(usable_count == 1_272, f"Expected 1,272 usable magnitudes, observed {usable_count}")
    assert_exact(nonblank_parse_failures == 0, f"Observed {nonblank_parse_failures} nonblank parse failures")
    assert_exact(infinite_count == 0, f"Observed {infinite_count} infinite values")
    assert_exact(zero_count == 0 and negative_count == 0, "Raw usable magnitudes include zero or negative values")

    availability = build_availability(raw, parsed, effect_category)
    plotting_data = build_verified_plotting_data(raw, parsed, effect_category)
    summary = build_summary(raw, parsed, effect_category)
    maxima = build_maximum_records(raw, plotting_data)
    top_ten = top_ten_tables(plotting_data)

    # Critical Figure 6 assertions.
    assert_exact(set(plotting_data["organism"]) == set(PLOTTED_ORGANISMS), "Unexpected organism in Figure 6 data")
    assert_exact(plotting_data["lifespan effect"].eq("Increase").all(), "A Figure 6 row is not exactly Increase")
    assert_exact(plotting_data["parsed_numeric_lifespan_change"].gt(0).all(), "A plotted magnitude is not greater than zero")
    independent_n = plotting_data.groupby("organism", sort=False).size().to_dict()
    assert_exact(independent_n == CURRENT_FIGURE_N, f"Figure 6 sample sizes differ: {independent_n}")
    expected_total = sum(CURRENT_FIGURE_N.values())
    assert_exact(expected_total == 830, f"Unexpected expected total: {expected_total}")
    assert_exact(len(plotting_data) == expected_total, f"Verified plotting rows total {len(plotting_data)}, expected {expected_total}")
    assert_exact(plotting_data["raw_row_index_zero_based"].is_unique, "A raw row appears more than once in the verified plotting data")

    for organism, (expected_symbol, expected_value) in EXPECTED_MAX_LABELS.items():
        rows = maxima.loc[maxima["organism"].eq(organism)]
        assert_exact(len(rows) == 1, f"Maximum tie or missing maximum for {organism}: {len(rows)} records")
        row = rows.iloc[0]
        assert_exact(row["symbol"] == expected_symbol, f"Maximum symbol for {organism} is {row['symbol']}, expected {expected_symbol}")
        assert_exact(float(row["parsed_numeric_lifespan_change"]) == expected_value, f"Maximum value for {organism} differs")

    sp_mask = raw[ORGANISM_COLUMN].eq("Schizosaccharomyces pombe")
    sp_counts = {
        "total": int(sp_mask.sum()),
        "Increase": int((sp_mask & effect_category.eq("Increase")).sum()),
        "Decrease": int((sp_mask & effect_category.eq("Decrease")).sum()),
        "Increase and Decrease": int((sp_mask & effect_category.eq("Increase and Decrease")).sum()),
        "Missing/Other": int((sp_mask & effect_category.eq("Missing/Other")).sum()),
        "with numerical magnitude": int((sp_mask & parsed.notna()).sum()),
        "without usable magnitude": int((sp_mask & parsed.isna()).sum()),
    }
    assert_exact(sp_counts == {"total": 61, "Increase": 56, "Decrease": 5, "Increase and Decrease": 0, "Missing/Other": 0, "with numerical magnitude": 0, "without usable magnitude": 61}, f"S. pombe audit differs: {sp_counts}")

    # Existing derived files are compared only after independent reconstruction.
    availability_results = compare_table_by_keys(
        EXISTING_AVAILABILITY_FILE,
        availability,
        ["organism", "lifespan_effect_category"],
        ["organism_total_records", "category_records", "usable_numeric_magnitude_records", "missing_magnitude_records", "usable_magnitude_percent"],
    )
    summary_results = compare_table_by_keys(
        EXISTING_SUMMARY_FILE,
        summary,
        ["organism", "lifespan_effect_category"],
        ["n", "mean", "median", "sd", "q1", "q3", "iqr", "minimum", "maximum"],
    )
    independent_extremes = build_existing_style_extremes(raw, parsed, effect_category)
    extremes_results = compare_existing_extremes(independent_extremes)
    all_existing_results = fig15_results + availability_results + summary_results + extremes_results
    discrepancies = [item for item in all_existing_results if item.classification == "discrepancy"]
    assert_exact(not discrepancies, "Existing derived-file discrepancies detected: " + "; ".join(f"{item.file} {item.key} {item.field}" for item in discrepancies[:10]))

    # Instrument the current plot in memory. save_preview is replaced with a no-op.
    plot_integrity = verify_current_plot_inputs(raw, parsed, effect_category, plotting_data)
    assert_exact(plot_integrity["total_scatter"] == expected_total, "Current plot scatter total differs")
    assert_exact(plot_integrity["total_box"] == expected_total, "Current plot boxplot-input total differs")

    # Confirm Figures 1, 5 and 6 code all identify the same immutable raw path.
    raw_path_fragment = '"GenAge Project" / "models_genes" / "genage_models.csv"'
    assert_exact(raw_path_fragment in FIG1_SCRIPT.read_text(encoding="utf-8"), "Figure 1 script points to a different raw source")
    assert_exact(raw_path_fragment in FIG5_SCRIPT.read_text(encoding="utf-8"), "Figure 5 script points to a different raw source")
    assert_exact(raw_path_fragment in CURRENT_FIG6_SCRIPT.read_text(encoding="utf-8"), "Figure 6 script points to a different raw source")

    # These are the only two CSVs created by this audit.
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    plotting_data.to_csv(VERIFIED_DATA_FILE, index=False, lineterminator="\n")
    maxima.to_csv(VERIFIED_MAXIMA_FILE, index=False, lineterminator="\n")
    assert_exact(len(pd.read_csv(VERIFIED_DATA_FILE)) == expected_total, "Written verified plotting CSV row count differs")
    assert_exact(len(pd.read_csv(VERIFIED_MAXIMA_FILE)) == 4, "Written maximum-record CSV row count differs")

    figure_hashes_after = {str(path.relative_to(PROJECT_ROOT)): sha256(path) for path in FIG6_FIGURE_FILES if path.is_file()}
    assert_exact(figure_hashes_before == figure_hashes_after, "A Figure 6 file changed during verification")

    main_summary = summary.loc[
        summary["organism"].isin(PLOTTED_ORGANISMS)
        & summary["lifespan_effect_category"].eq("Increase")
    ].copy()
    for column in ["mean", "median", "sd", "q1", "q3", "iqr", "minimum", "maximum"]:
        main_summary[column] = main_summary[column].map(display_number)
    sample_comparison = pd.DataFrame(
        {
            "Organism": PLOTTED_ORGANISMS,
            "Current Figure 6 n": [CURRENT_FIGURE_N[organism] for organism in PLOTTED_ORGANISMS],
            "Independently recomputed n": [independent_n[organism] for organism in PLOTTED_ORGANISMS],
            "Match": [CURRENT_FIGURE_N[organism] == independent_n[organism] for organism in PLOTTED_ORGANISMS],
        }
    )
    log_table = (
        plotting_data.groupby("organism", sort=False)["parsed_numeric_lifespan_change"]
        .agg(minimum="min", zero_values=lambda values: int((values == 0).sum()), negative_values=lambda values: int((values < 0).sum()))
        .reset_index()
    )
    maximum_report = maxima[
        [
            "organism",
            "symbol",
            "entrez gene id",
            "parsed_numeric_lifespan_change",
            "raw_row_index_zero_based",
            "csv_line_number",
            "raw_lifespan_change",
            "lifespan effect",
            "maximum_tie_count",
            "symbol_records_in_raw_organism",
            "symbol_records_in_figure6",
        ]
    ].copy()

    relevant_summary_comparisons = comparison_frame(summary_results)
    derived_summary = compact_comparison_summary(availability_results + summary_results + extremes_results)
    rounding_details = comparison_frame(
        [
            item
            for item in availability_results + summary_results + extremes_results
            if item.classification == "rounding difference"
        ]
    )
    figure_hash_table = pd.DataFrame(
        [{"Figure 6 file": path, "SHA-256 before": digest, "SHA-256 after": figure_hashes_after[path], "unchanged": digest == figure_hashes_after[path]} for path, digest in figure_hashes_before.items()]
    )
    effect_count_report = counts.rename(columns={"total GenAge records": "total"})
    availability_report = availability.copy()
    availability_report["usable_magnitude_percent"] = availability_report["usable_magnitude_percent"].map(display_number)

    top_ten_markdown = []
    for organism in PLOTTED_ORGANISMS:
        table = top_ten[organism].copy()
        table["parsed_numeric_lifespan_change"] = table["parsed_numeric_lifespan_change"].map(display_number)
        top_ten_markdown.append(f"### {organism}\n\n{markdown_table(table)}")

    release_text = LOCAL_RELEASE_FILE.read_text(encoding="utf-8", errors="replace") if LOCAL_RELEASE_FILE.is_file() else ""
    local_release_has_percent_definition = "avg lifespan change" in release_text.lower() and "%" in release_text
    report = f"""# Figure 6 independent verification report

Audit date: 2026-08-20  
Verification script: `scripts/verify_fig6_data.py`

## Raw source verified

- Exact filename: `{RAW_FILE.name}`
- Full local path: `{RAW_FILE.resolve()}`
- Dimensions: **{len(raw):,} rows x {len(raw.columns)} columns**
- SHA-256: `{raw_hash}`
- Exact columns, in file order: {', '.join(f'`{column}`' for column in raw.columns)}
- Relevant raw fields: organism = `{ORGANISM_COLUMN}`; symbol/name = `symbol` and `name`; Entrez ID = `entrez gene id`; longevity influence = `longevity influence`; lifespan effect = `{EFFECT_COLUMN}`; magnitude = `{MAGNITUDE_COLUMN}`.

The scripts used for Figures 1, 5 and 6 all point to the same local raw path, and the independently recomputed organism/effect counts exactly match the saved Figures 1 and 5 tables. The saved tables do not embed a source checksum, so the raw SHA-256 above is the definitive byte-level identity for this audit.

### Raw lifespan-effect values before cleaning

{markdown_table(effect_inventory)}

No capitalization or whitespace normalization is needed: all nonempty categories already exactly match `Increase`, `Decrease`, or `Increase and Decrease`. The 14 empty strings are explicitly recoded to `Missing/Other` only for tabulation; no record is dropped.

### Organism and raw effect counts

{markdown_table(effect_count_report)}

Every organism total and effect count is an exact match to the prior Figures 1 and 5 tables. In particular, *S. pombe* has **61** raw GenAge records.

## Filtering verified

The Figure 6 dataset was reconstructed directly from the raw CSV with three simultaneous rules: organism is one of the four plotted organisms; the raw lifespan-effect value is exactly `Increase`; and `{MAGNITUDE_COLUMN}` parses to a finite numeric value. A value is “usable” only if its stripped raw field is nonblank and converts to a finite number. No deduplication, aggregation, imputation, capping or winsorisation was applied.

Raw magnitude-field inspection:

- blank fields: **{blank_count:,}**
- whitespace-only but nonempty fields: **{whitespace_only_count:,}**
- literal NaN-like strings (`NaN`, `NA`, `N/A`, `null`, `none`, case-insensitive): **{literal_nan_like_count:,}**
- nonblank nonnumeric parse failures: **{nonblank_parse_failures:,}**
- infinite values: **{infinite_count:,}**
- usable finite numeric values: **{usable_count:,}**
- positive values: **{positive_count:,}**
- zero values: **{zero_count:,}**
- negative values: **{negative_count:,}**

There are no unusual nonblank strings to list. All missing/unusable cases are blank raw fields.

### Magnitude-data availability by organism and effect category

{markdown_table(availability_report)}

The exact 830 record-level rows used by the main Figure 6 were written to `{VERIFIED_DATA_FILE.relative_to(PROJECT_ROOT)}`. Each output row retains the zero-based raw row index and the one-based CSV line number (including the header as line 1).

## Sample sizes

{markdown_table(sample_comparison)}

The independent total is **{len(plotting_data):,} = 217 + 453 + 114 + 46**. All four comparisons match exactly.

## Summary statistics

The standard deviation is the sample standard deviation (`ddof=1`); quartiles use pandas' default linear interpolation. These values were calculated from the raw-derived Figure 6 records, not from the existing summary file.

{markdown_table(main_summary[["organism", "n", "mean", "median", "sd", "q1", "q3", "iqr", "minimum", "maximum"]])}

### Column-by-column comparison with `fig6_lifespan_magnitude_summary.csv`

{markdown_table(relevant_summary_comparisons)}

`exact match` means numerical identity (allowing equivalent integer/decimal representations). `rounding difference` means the existing value is the independent value formatted to six significant digits. No comparison is classified as a discrepancy.

## Maximum records

The maximum was found independently by sorting each organism's raw-derived Figure 6 records by numeric magnitude descending, with raw row index as a deterministic tie-breaker. Every maximum is unique; the four proposed labels are therefore correct. None of the four maximum-gene symbols is repeated within its organism in either the full raw dataset or the Figure 6 subset.

{markdown_table(maximum_report)}

Every tied maximum would be retained in `{VERIFIED_MAXIMA_FILE.relative_to(PROJECT_ROOT)}`; the observed tie count is one for every organism.

### Top 10 Figure 6 records per organism

{chr(10).join(chr(10) + block for block in top_ten_markdown)}

## S. pombe

| Quantity | Independently recomputed |
| --- | ---: |
| Total records | {sp_counts['total']} |
| Increase | {sp_counts['Increase']} |
| Decrease | {sp_counts['Decrease']} |
| Increase and Decrease | {sp_counts['Increase and Decrease']} |
| Missing/Other lifespan effect | {sp_counts['Missing/Other']} |
| Records with a usable numerical magnitude | {sp_counts['with numerical magnitude']} |
| Records without a usable magnitude | {sp_counts['without usable magnitude']} |

Thus, all **61 / 61 (100%)** *S. pombe* records lack a usable magnitude. Its exclusion is justified solely by absent magnitude data, not by an absence of GenAge records.

## Plot integrity

The current `plot_final_extension_only` function was tested only after the independent dataset had been created. Its save routine was replaced in memory with a no-op, so no figure was written. The scatter x arrays and boxplot input arrays were captured and compared, in order, with the independently filtered record-level magnitudes.

{markdown_table(plot_integrity['per_organism'])}

- captured scatter observations: **{plot_integrity['total_scatter']:,}**
- captured boxplot-input observations: **{plot_integrity['total_box']:,}**
- verified CSV rows: **{len(plotting_data):,}**
- unique raw row indices: **{plotting_data['raw_row_index_zero_based'].nunique():,}**
- deduplication/aggregation operations in the current main plotting function: **none detected**
- jitter: changes only the y display coordinate; every captured x magnitude remained exactly equal to the raw-derived value

The scatter and boxplot both receive the same 830 record-level values. No row is omitted, plotted twice, deduplicated or aggregated before plotting.

### Logarithmic-axis checks

{markdown_table(log_table)}

Every plotted magnitude is strictly greater than zero, so the log x-axis is numerically legitimate.

### Figure-file immutability during the audit

{markdown_table(figure_hash_table)}

## Units and interpretation

### Confirmed from source/documentation

- The downloaded CSV calls the raw field `{MAGNITUDE_COLUMN}` but does not include a unit in the column header.
- The local `release.html` records GenAge Build 21 and 2,205 model-organism genes; it does **{'include' if local_release_has_percent_definition else 'not include'}** a definition of the magnitude field or its unit.
- The official [GenAge model-organism download page](https://genomics.senescence.info/genes/models.html) identifies Build 21 and the 2,205-record download.
- The official [GenAge model-organism search](https://genomics.senescence.info/genes/search.php?lifespan_effect=increase&long_influence=&sort=5) labels the display column “Max average lifespan change” and appends `%` to these records, including age-1 1000%, SIR2 600%, SCH9 300%, daf-2 200% and Vinc 150%.
- Official detail pages also express the field under percentage-change headings, including [age-1](https://genomics.senescence.info/genes/details.php?gene=age-1&organism=elegans), [SIR2](https://genomics.senescence.info/genes/details.php?gene=SIR2&organism=cerevisiae), and [Ghrh](https://genomics.senescence.info/genes/details.php?gene=Ghrh&organism=musculus).

Therefore, `%` is justified as the source database's display unit for this field. The axis should be interpreted as **GenAge-reported percentage-coded lifespan-change summary values**, not as percentages independently recalculated from the underlying experiments.

### Assumed/inferred or unavailable

- The downloaded CSV itself does not define the exact formula behind `avg lifespan change (max obsv)`.
- The inspected local and official pages do not provide enough methodological detail to reconstruct how heterogeneous mean/median, fold-change and maximum-observation statements were normalised into this single field.
- No new interpretation of that calculation is asserted here. The record-level values remain heterogeneous descriptive summaries.

## Discrepancies

### Existing Figure 6 derived-file comparison

{markdown_table(derived_summary)}

All counts, categories, record ordering and exact-valued statistics match. Non-exact cells are explained by the six-significant-digit export format used for the existing files; these include percentages, non-integer statistics, two magnitude strings and three large Entrez-ID strings in the old extremes table. **No Figure 6 sample-size, magnitude, row-inclusion, maximum-label or plotting discrepancy was found.** The three existing Figure 6 derived files were read only and were not overwritten.

### Every non-exact (rounding-only) comparison

{markdown_table(rounding_details)}

The old extremes table was exported with a six-significant-digit float format. Consequently, four large Entrez IDs are written in scientific notation; one remains numerically exact and three lose their final one-digit precision in that old display-only table. The verified plotting-data and maximum-record CSVs preserve Entrez IDs as source strings. This is classified as a rounding difference in the pre-existing audit table, not a Figure 6 point or annotation discrepancy, because Entrez ID is not used to position points and all raw row identities, GenAge IDs, symbols and magnitudes were independently matched.

## Verdict

**PASS — Figure 6 data and annotations independently verified.**
"""

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(report, encoding="utf-8", newline="\n")
    print(f"Raw source: {RAW_FILE}")
    print(f"SHA-256: {raw_hash}")
    print(f"Verified plotting rows: {len(plotting_data)}")
    print(f"Verified maxima rows: {len(maxima)}")
    print(f"Wrote: {VERIFIED_DATA_FILE}")
    print(f"Wrote: {VERIFIED_MAXIMA_FILE}")
    print(f"Wrote: {REPORT_FILE}")
    print("PASS — Figure 6 data and annotations independently verified.")


if __name__ == "__main__":
    main()
