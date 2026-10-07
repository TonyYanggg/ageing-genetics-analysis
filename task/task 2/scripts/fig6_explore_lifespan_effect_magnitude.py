"""Audit GenAge lifespan-effect magnitudes and create Figure 6 previews.

This is an exploratory, record-level analysis. It does not deduplicate genes,
alter the source column, modify Figures 1--5b, or create a final Figure 6.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


TASK_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = PROJECT_ROOT / "GenAge Project" / "models_genes" / "genage_models.csv"
TABLE_DIR = TASK_ROOT / "results" / "tables"
FIGURE_DIR = TASK_ROOT / "figures" / "exploratory"
NOTES_FILE = TASK_ROOT / "notes" / "fig6_exploration_notes.md"
FIGURE_NOTE_FILE = TASK_ROOT / "notes" / "fig6_figure_note.md"
FIGURE_CAPTION_FILE = TASK_ROOT / "notes" / "fig6_caption.md"

ORGANISM_COLUMN = "organism"
EFFECT_COLUMN = "lifespan effect"
MAGNITUDE_COLUMN = "avg lifespan change (max obsv)"

EXPECTED_COLUMNS = {
    "GenAge ID",
    "symbol",
    "name",
    ORGANISM_COLUMN,
    "entrez gene id",
    MAGNITUDE_COLUMN,
    EFFECT_COLUMN,
    "longevity influence",
}
EXPECTED_TOTAL_RECORDS = 2_205
FULL_ORGANISM_ORDER = [
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
EXPECTED_ORGANISM_COUNTS = {
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
MAJOR_ORGANISMS = FULL_ORGANISM_ORDER[:5]
FINAL_PLOTTED_ORGANISMS = MAJOR_ORGANISMS[:4]
EXPECTED_FINAL_USABLE_COUNTS = {
    "Saccharomyces cerevisiae": {"Increase": 217, "Decrease": 218},
    "Caenorhabditis elegans": {"Increase": 453, "Decrease": 163},
    "Drosophila melanogaster": {"Increase": 114, "Decrease": 33},
    "Mus musculus": {"Increase": 46, "Decrease": 25},
}
SUMMARY_EFFECT_CATEGORIES = ["Increase", "Decrease", "Increase and Decrease"]
AVAILABILITY_EFFECT_CATEGORIES = [
    "Increase",
    "Decrease",
    "Increase and Decrease",
    "Missing/Other",
]

INCREASE_COLOR = "#59A14F"
DECREASE_COLOR = "#E15759"
POINT_EDGE = "#FFFFFF"
NEUTRAL_COLOR = "#4E79A7"
RANDOM_SEED = 20260820


def configure_plot_style() -> None:
    """Apply the restrained style used by the existing GenAge figures."""

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.axisbelow": True,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def load_and_validate_data() -> pd.DataFrame:
    """Load the immutable source and add diagnostic columns to a copy."""

    if not INPUT_FILE.is_file():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    data = pd.read_csv(INPUT_FILE)
    missing_columns = EXPECTED_COLUMNS.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Input is missing expected columns: {sorted(missing_columns)}")
    if len(data) != EXPECTED_TOTAL_RECORDS:
        raise ValueError(
            f"Expected {EXPECTED_TOTAL_RECORDS:,} records; observed {len(data):,}."
        )

    observed_counts = data[ORGANISM_COLUMN].value_counts().to_dict()
    if observed_counts != EXPECTED_ORGANISM_COUNTS:
        raise ValueError(
            "Organism counts differ from the validated inventory. "
            f"Observed: {observed_counts}"
        )

    effect_clean = data[EFFECT_COLUMN].astype("string").str.strip()
    data = data.copy()
    data["effect_category"] = effect_clean.where(
        effect_clean.isin(SUMMARY_EFFECT_CATEGORIES), "Missing/Other"
    )
    data["numeric_lifespan_change"] = pd.to_numeric(
        data[MAGNITUDE_COLUMN], errors="coerce"
    ).replace([np.inf, -np.inf], np.nan)

    # The raw field is preserved. Sign is reconstructed only in this new column.
    magnitude = data["numeric_lifespan_change"].abs()
    data["signed_lifespan_change"] = np.select(
        [data["effect_category"].eq("Increase"), data["effect_category"].eq("Decrease")],
        [magnitude, -magnitude],
        default=np.nan,
    )
    return data


def print_source_audit(data: pd.DataFrame) -> None:
    """Print the source path, schema, and exact analysis field names."""

    print("Step 1 -- source data audit")
    print(f"File path used: {INPUT_FILE.resolve()}")
    print(f"Number of rows: {len(data):,}")
    print("Column names:")
    for column in data.columns[: len(data.columns) - 3]:
        print(f"  - {column}")
    print("Data types (raw columns):")
    raw_columns = [column for column in data.columns if column not in {
        "effect_category", "numeric_lifespan_change", "signed_lifespan_change"
    }]
    for column, dtype in data[raw_columns].dtypes.items():
        print(f"  - {column}: {dtype}")
    print("Exact fields used:")
    print(f"  - organism: {ORGANISM_COLUMN}")
    print(f"  - lifespan effect: {EFFECT_COLUMN}")
    print(f"  - magnitude: {MAGNITUDE_COLUMN}")
    print("No gene deduplication was performed.\n")


def build_availability_table(data: pd.DataFrame) -> pd.DataFrame:
    """Count magnitude availability for every organism and effect category."""

    rows: list[dict[str, object]] = []
    for organism in FULL_ORGANISM_ORDER:
        organism_data = data.loc[data[ORGANISM_COLUMN].eq(organism)]
        total_records = len(organism_data)
        for category in AVAILABILITY_EFFECT_CATEGORIES:
            category_data = organism_data.loc[
                organism_data["effect_category"].eq(category)
            ]
            usable = int(category_data["numeric_lifespan_change"].notna().sum())
            category_records = len(category_data)
            rows.append(
                {
                    "organism": organism,
                    "organism_total_records": total_records,
                    "lifespan_effect_category": category,
                    "category_records": category_records,
                    "usable_numeric_magnitude_records": usable,
                    "missing_magnitude_records": category_records - usable,
                    "usable_magnitude_percent": (
                        100 * usable / category_records if category_records else np.nan
                    ),
                }
            )

    availability = pd.DataFrame(rows)
    totals = availability.groupby("organism", sort=False)["category_records"].sum()
    expected = pd.Series(EXPECTED_ORGANISM_COUNTS).reindex(FULL_ORGANISM_ORDER)
    if not totals.equals(expected):
        raise ValueError("Availability-table category counts do not match organism totals.")
    return availability


def describe_values(values: pd.Series) -> dict[str, float | int]:
    """Return the requested distribution statistics for numeric values."""

    clean = values.dropna().astype(float)
    if clean.empty:
        return {
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
    q1 = float(clean.quantile(0.25))
    q3 = float(clean.quantile(0.75))
    return {
        "n": int(clean.size),
        "mean": float(clean.mean()),
        "median": float(clean.median()),
        "sd": float(clean.std(ddof=1)) if clean.size > 1 else np.nan,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "minimum": float(clean.min()),
        "maximum": float(clean.max()),
    }


def build_magnitude_summary(data: pd.DataFrame) -> pd.DataFrame:
    """Summarize unsigned magnitudes by major organism and effect category."""

    rows: list[dict[str, object]] = []
    for organism in MAJOR_ORGANISMS:
        for category in SUMMARY_EFFECT_CATEGORIES:
            values = data.loc[
                data[ORGANISM_COLUMN].eq(organism)
                & data["effect_category"].eq(category),
                "numeric_lifespan_change",
            ]
            rows.append(
                {
                    "organism": organism,
                    "lifespan_effect_category": category,
                    **describe_values(values),
                }
            )
    return pd.DataFrame(rows)


def build_extremes_table(data: pd.DataFrame) -> pd.DataFrame:
    """Collect the top ten usable records in each major-organism/effect group."""

    identifying_columns = [
        ORGANISM_COLUMN,
        "symbol",
        "name",
        "entrez gene id",
        "longevity influence",
        EFFECT_COLUMN,
        MAGNITUDE_COLUMN,
    ]
    groups: list[pd.DataFrame] = []
    for organism in MAJOR_ORGANISMS:
        for category in SUMMARY_EFFECT_CATEGORIES:
            group = data.loc[
                data[ORGANISM_COLUMN].eq(organism)
                & data["effect_category"].eq(category)
                & data["numeric_lifespan_change"].notna()
            ].sort_values(
                ["numeric_lifespan_change", "GenAge ID"], ascending=[False, True]
            )
            if group.empty:
                continue
            top = group.head(10)[identifying_columns + ["GenAge ID"]].copy()
            top.insert(0, "rank_within_organism_and_effect", range(1, len(top) + 1))
            top.insert(1, "usable_records_in_group", len(group))
            groups.append(top)

    if not groups:
        return pd.DataFrame(
            columns=[
                "rank_within_organism_and_effect",
                "usable_records_in_group",
                *identifying_columns,
                "GenAge ID",
            ]
        )
    return pd.concat(groups, ignore_index=True)


def build_signed_summary(data: pd.DataFrame) -> pd.DataFrame:
    """Describe the diagnostic signed variable by major organism."""

    rows: list[dict[str, object]] = []
    for organism in MAJOR_ORGANISMS:
        values = data.loc[
            data[ORGANISM_COLUMN].eq(organism), "signed_lifespan_change"
        ]
        rows.append({"organism": organism, **describe_values(values)})
    return pd.DataFrame(rows)


def style_species_axis(
    ax: plt.Axes,
    organisms: list[str] = MAJOR_ORGANISMS,
) -> None:
    """Apply shared species-axis styling for the requested organism set."""

    positions = np.arange(len(organisms))
    ax.set_yticks(positions, labels=organisms)
    ax.set_ylim(-0.65, len(organisms) - 0.35)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#777777")
    for label in ax.get_yticklabels():
        label.set_fontstyle("italic")


def add_horizontal_boxplot(
    ax: plt.Axes,
    values: np.ndarray,
    position: float,
    color: str,
    width: float,
) -> None:
    """Draw a compact horizontal boxplot; all observations are plotted separately."""

    if values.size == 0:
        return
    ax.boxplot(
        [values],
        positions=[position],
        widths=width,
        orientation="horizontal",
        patch_artist=True,
        manage_ticks=False,
        showfliers=False,
        boxprops={"facecolor": color, "alpha": 0.26, "edgecolor": color, "linewidth": 1.0},
        whiskerprops={"color": color, "linewidth": 1.0},
        capprops={"color": color, "linewidth": 1.0},
        medianprops={"color": "#202020", "linewidth": 1.25},
    )


def save_preview(fig: plt.Figure, basename: str) -> None:
    """Save the diagnostic preview as PNG plus editable-text SVG."""

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / f"{basename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{basename}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_extension_only(data: pd.DataFrame) -> None:
    """Preview A: Increase records only, with every usable observation visible."""

    rng = np.random.default_rng(RANDOM_SEED)
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    for position, organism in enumerate(MAJOR_ORGANISMS):
        values = data.loc[
            data[ORGANISM_COLUMN].eq(organism)
            & data["effect_category"].eq("Increase")
            & data["numeric_lifespan_change"].notna(),
            "numeric_lifespan_change",
        ].to_numpy(dtype=float)
        add_horizontal_boxplot(ax, values, position, INCREASE_COLOR, width=0.46)
        if values.size:
            jitter = np.clip(rng.normal(0, 0.075, values.size), -0.20, 0.20)
            ax.scatter(
                values,
                position + jitter,
                s=12,
                color=INCREASE_COLOR,
                alpha=0.34,
                edgecolor=POINT_EDGE,
                linewidth=0.20,
                rasterized=False,
                zorder=3,
            )
        ax.text(
            1.015,
            position,
            f"n={values.size:,}",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.8,
            clip_on=False,
        )

    ax.set_xscale("log")
    ax.set_xlim(3.5, 1_200)
    ticks = [5, 10, 20, 50, 100, 200, 500, 1_000]
    ax.set_xticks(ticks, labels=[f"{tick:,}" for tick in ticks])
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Reported lifespan-change magnitude (source values; log scale)", labelpad=9)
    ax.set_title(
        "Diagnostic preview A — extension-only magnitude distributions",
        loc="left",
        pad=12,
    )
    style_species_axis(ax)
    fig.text(
        0.32,
        0.035,
        "Increase records only. Boxes show median and IQR; points show every usable record. This is not the final Figure 6.",
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#4A4A4A",
    )
    fig.subplots_adjust(left=0.30, right=0.84, top=0.88, bottom=0.20)
    save_preview(fig, "fig6_preview_A_extension_only")


def plot_signed_effects(data: pd.DataFrame) -> None:
    """Preview B: reconstruct Increase/Decrease signs and retain every observation."""

    rng = np.random.default_rng(RANDOM_SEED + 1)
    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    category_specs = [
        ("Increase", -0.15, INCREASE_COLOR),
        ("Decrease", 0.15, DECREASE_COLOR),
    ]
    for position, organism in enumerate(MAJOR_ORGANISMS):
        counts: dict[str, int] = {}
        for category, offset, color in category_specs:
            values = data.loc[
                data[ORGANISM_COLUMN].eq(organism)
                & data["effect_category"].eq(category)
                & data["signed_lifespan_change"].notna(),
                "signed_lifespan_change",
            ].to_numpy(dtype=float)
            counts[category] = values.size
            add_horizontal_boxplot(ax, values, position + offset, color, width=0.22)
            if values.size:
                jitter = np.clip(rng.normal(0, 0.035, values.size), -0.085, 0.085)
                ax.scatter(
                    values,
                    position + offset + jitter,
                    s=11,
                    color=color,
                    alpha=0.36,
                    edgecolor=POINT_EDGE,
                    linewidth=0.18,
                    rasterized=False,
                    zorder=3,
                )
        ax.text(
            1.015,
            position,
            f"Increase n={counts['Increase']:,}; Decrease n={counts['Decrease']:,}",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.0,
            clip_on=False,
        )

    ax.axvline(0, color="#222222", linewidth=1.0, zorder=2)
    ax.set_xscale("symlog", linthresh=10, linscale=0.8, base=10)
    ax.set_xlim(-1_200, 1_200)
    ticks = [-1_000, -100, -10, 0, 10, 100, 1_000]
    ax.set_xticks(ticks, labels=[f"{tick:,}" for tick in ticks])
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(
        "Exploratory signed lifespan change (source values; symmetric log scale)",
        labelpad=9,
    )
    ax.set_title(
        "Diagnostic preview B — signed Increase versus Decrease effects",
        loc="left",
        pad=12,
    )
    style_species_axis(ax)
    ax.legend(
        handles=[
            Patch(facecolor=INCREASE_COLOR, alpha=0.65, label="Increase (positive)"),
            Patch(facecolor=DECREASE_COLOR, alpha=0.65, label="Decrease (negative)"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.27),
        ncol=2,
        fontsize=8.8,
    )
    fig.text(
        0.32,
        0.035,
        "Sign reconstructed from the categorical lifespan-effect field; mixed-direction records are undefined. This is not the final Figure 6.",
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#4A4A4A",
    )
    fig.subplots_adjust(left=0.28, right=0.79, top=0.88, bottom=0.25)
    save_preview(fig, "fig6_preview_B_signed_effects")


def plot_final_extension_only(data: pd.DataFrame) -> None:
    """Export the primary Figure 6 using Increase records with usable magnitudes."""

    rng = np.random.default_rng(RANDOM_SEED)
    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    extreme_coordinates: dict[tuple[str, str], tuple[float, float]] = {}
    annotated_symbols = {"SIR2", "age-1", "Vinc", "Ghrh"}

    for position, organism in enumerate(FINAL_PLOTTED_ORGANISMS):
        plot_data = data.loc[
            data[ORGANISM_COLUMN].eq(organism)
            & data["effect_category"].eq("Increase")
            & data["numeric_lifespan_change"].notna(),
            ["symbol", "numeric_lifespan_change"],
        ].reset_index(drop=True)
        values = plot_data["numeric_lifespan_change"].to_numpy(dtype=float)
        expected = EXPECTED_FINAL_USABLE_COUNTS[organism]["Increase"]
        if values.size != expected:
            raise ValueError(
                f"Expected {expected} usable Increase records for {organism}; "
                f"observed {values.size}."
            )

        add_horizontal_boxplot(ax, values, position, INCREASE_COLOR, width=0.46)
        jitter = np.clip(rng.normal(0, 0.075, values.size), -0.20, 0.20)
        point_y = position + jitter
        ax.scatter(
            values,
            point_y,
            s=12,
            color=INCREASE_COLOR,
            alpha=0.34,
            edgecolor=POINT_EDGE,
            linewidth=0.20,
            rasterized=False,
            zorder=3,
        )
        for row_index, row in plot_data.iterrows():
            symbol = str(row["symbol"])
            if symbol in annotated_symbols:
                extreme_coordinates[(organism, symbol)] = (
                    float(row["numeric_lifespan_change"]),
                    float(point_y[row_index]),
                )

        ax.text(
            1.015,
            position,
            f"n = {values.size:,}",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.8,
            clip_on=False,
        )

    annotation_specs = {
        ("Saccharomyces cerevisiae", "SIR2"): ("SIR2 (600%)", 600, -10, -18, "right"),
        ("Caenorhabditis elegans", "age-1"): ("age-1 (1,000%)", 1_000, -10, -18, "right"),
        ("Drosophila melanogaster", "Vinc"): ("Vinc (150%)", 150, 8, 16, "left"),
        ("Mus musculus", "Ghrh"): ("Ghrh (46%)", 46, 9, 16, "left"),
    }
    for key, (label, expected_value, x_offset, y_offset, alignment) in annotation_specs.items():
        if key not in extreme_coordinates:
            raise ValueError(f"Expected extreme record was not found: {key}")
        x_value, y_value = extreme_coordinates[key]
        if not np.isclose(x_value, expected_value):
            raise ValueError(
                f"Expected {key} to have magnitude {expected_value}; observed {x_value}."
            )
        ax.annotate(
            label,
            xy=(x_value, y_value),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            ha=alignment,
            va="center",
            fontsize=7.8,
            color="#303030",
            arrowprops={"arrowstyle": "-", "color": "#666666", "linewidth": 0.65},
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.8},
            zorder=5,
        )

    ax.set_xscale("log")
    ax.set_xlim(3.5, 1_200)
    ticks = [5, 10, 20, 50, 100, 200, 500, 1_000]
    ax.set_xticks(ticks, labels=[f"{tick:,}" for tick in ticks])
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Reported lifespan-extension magnitude (%)", labelpad=9)
    ax.set_title(
        "Reported lifespan-extension magnitude across model organisms",
        loc="left",
        pad=12,
    )
    style_species_axis(ax, FINAL_PLOTTED_ORGANISMS)
    fig.subplots_adjust(left=0.30, right=0.82, top=0.86, bottom=0.17)
    save_preview(fig, "fig6_lifespan_extension_magnitude")


def plot_final_signed_effects(data: pd.DataFrame) -> None:
    """Export the alternative signed Figure 6 view with reconstructed direction."""

    rng = np.random.default_rng(RANDOM_SEED + 1)
    fig, ax = plt.subplots(figsize=(10.2, 5.2))
    category_specs = [
        ("Increase", -0.15, INCREASE_COLOR),
        ("Decrease", 0.15, DECREASE_COLOR),
    ]
    for position, organism in enumerate(FINAL_PLOTTED_ORGANISMS):
        counts: dict[str, int] = {}
        for category, offset, color in category_specs:
            values = data.loc[
                data[ORGANISM_COLUMN].eq(organism)
                & data["effect_category"].eq(category)
                & data["signed_lifespan_change"].notna(),
                "signed_lifespan_change",
            ].to_numpy(dtype=float)
            counts[category] = values.size
            expected = EXPECTED_FINAL_USABLE_COUNTS[organism][category]
            if values.size != expected:
                raise ValueError(
                    f"Expected {expected} usable {category} records for {organism}; "
                    f"observed {values.size}."
                )
            add_horizontal_boxplot(ax, values, position + offset, color, width=0.22)
            jitter = np.clip(rng.normal(0, 0.035, values.size), -0.085, 0.085)
            ax.scatter(
                values,
                position + offset + jitter,
                s=11,
                color=color,
                alpha=0.36,
                edgecolor=POINT_EDGE,
                linewidth=0.18,
                rasterized=False,
                zorder=3,
            )
        ax.text(
            1.015,
            position,
            f"usable n: Increase = {counts['Increase']:,}; Decrease = {counts['Decrease']:,}",
            transform=ax.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.0,
            clip_on=False,
        )

    ax.axvline(0, color="#222222", linewidth=1.0, zorder=2)
    ax.set_xscale("symlog", linthresh=10, linscale=0.8, base=10)
    ax.set_xlim(-1_200, 1_200)
    ticks = [-1_000, -100, -10, 0, 10, 100, 1_000]
    ax.set_xticks(ticks, labels=[f"{tick:,}" for tick in ticks])
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.7)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Reconstructed signed lifespan-effect magnitude (%)", labelpad=9)
    ax.set_title(
        "Signed lifespan-effect magnitude across model organisms",
        loc="left",
        pad=12,
    )
    style_species_axis(ax, FINAL_PLOTTED_ORGANISMS)
    ax.legend(
        handles=[
            Patch(facecolor=INCREASE_COLOR, alpha=0.65, label="Increase (positive)"),
            Patch(facecolor=DECREASE_COLOR, alpha=0.65, label="Decrease (negative)"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.31),
        ncol=2,
        fontsize=8.8,
    )
    fig.text(
        0.28,
        0.050,
        "Sign reconstructed from lifespan effect: Increase → positive; Decrease → negative; Increase and Decrease excluded because direction is not unique.",
        ha="left",
        va="bottom",
        fontsize=8.0,
        color="#4A4A4A",
    )
    fig.text(
        0.28,
        0.020,
        "Usable n is direction-specific because magnitude-data availability differs. S. pombe omitted: none of 61 records had a usable magnitude.",
        ha="left",
        va="bottom",
        fontsize=8.0,
        color="#4A4A4A",
    )
    fig.subplots_adjust(left=0.28, right=0.76, top=0.87, bottom=0.30)
    save_preview(fig, "fig6_signed_lifespan_effect_magnitude")


def write_final_figure_note() -> None:
    """Write the short manuscript-facing note requested for the main Figure 6."""

    note = """# Figure 6 note

**Fig. 6 | Reported lifespan-extension magnitude across model organisms.** Only GenAge records classified as `Increase` and containing a usable numerical lifespan-change magnitude are plotted. Boxes show the median and interquartile range; points represent individual GenAge records, including extreme observations. `Usable n` therefore differs from the total number of GenAge records. *Schizosaccharomyces pombe* was excluded because none of its 61 records contained a usable magnitude in the downloaded dataset. These record-level summary values derive from heterogeneous experiments and should be interpreted descriptively.
"""
    FIGURE_NOTE_FILE.parent.mkdir(parents=True, exist_ok=True)
    FIGURE_NOTE_FILE.write_text(note, encoding="utf-8")


def write_main_figure_caption() -> None:
    """Write the standalone caption for the finalized main Figure 6."""

    caption = """# Figure 6 caption

Figure 6. Reported lifespan-extension magnitude across model organisms. GenAge model-organism records classified as Increase and containing an available numerical lifespan-change magnitude were included. Boxes show the median and interquartile range, and points represent individual GenAge records. The x-axis is logarithmic. n denotes the number of Increase records with an available numerical lifespan-change magnitude. Labels identify the maximum reported value for each plotted organism. Schizosaccharomyces pombe was not plotted because none of its 61 GenAge records contained an available lifespan-change magnitude in the downloaded dataset.
"""
    FIGURE_CAPTION_FILE.parent.mkdir(parents=True, exist_ok=True)
    FIGURE_CAPTION_FILE.write_text(caption, encoding="utf-8")


def export_final_main_figure_only() -> None:
    """Update only the main Figure 6 artwork and its standalone caption."""

    data = load_and_validate_data()
    configure_plot_style()
    plot_final_extension_only(data)
    write_main_figure_caption()
    print("Main Figure 6 updated without writing analysis tables or the signed figure:")
    print(f"  - {FIGURE_DIR / 'fig6_lifespan_extension_magnitude.png'}")
    print(f"  - {FIGURE_DIR / 'fig6_lifespan_extension_magnitude.svg'}")
    print(f"  - {FIGURE_CAPTION_FILE}")


def export_final_figures_only() -> None:
    """Finalize Figure 6 plots without rebuilding or writing analysis tables."""

    data = load_and_validate_data()
    configure_plot_style()
    plot_final_extension_only(data)
    plot_final_signed_effects(data)
    write_final_figure_note()
    print("Final Figure 6 exports written without regenerating CSV tables:")
    print(f"  - {FIGURE_DIR / 'fig6_lifespan_extension_magnitude.png'}")
    print(f"  - {FIGURE_DIR / 'fig6_lifespan_extension_magnitude.svg'}")
    print(f"  - {FIGURE_DIR / 'fig6_signed_lifespan_effect_magnitude.png'}")
    print(f"  - {FIGURE_DIR / 'fig6_signed_lifespan_effect_magnitude.svg'}")
    print(f"  - {FIGURE_NOTE_FILE}")


def markdown_table(table: pd.DataFrame, formats: dict[str, str] | None = None) -> str:
    """Render a small DataFrame as Markdown without an extra dependency."""

    formats = formats or {}
    headers = [str(column) for column in table.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in table.iterrows():
        cells: list[str] = []
        for column in table.columns:
            value = row[column]
            if pd.isna(value):
                rendered = "NA"
            elif column in formats:
                rendered = formats[column].format(value)
            else:
                rendered = str(value)
            cells.append(rendered.replace("|", "\\|"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_notes(
    data: pd.DataFrame,
    availability: pd.DataFrame,
    magnitude_summary: pd.DataFrame,
    signed_summary: pd.DataFrame,
) -> None:
    """Write the requested descriptive exploration notes."""

    overall = (
        availability.groupby("lifespan_effect_category", sort=False)
        .agg(
            records=("category_records", "sum"),
            usable=("usable_numeric_magnitude_records", "sum"),
            missing_magnitude=("missing_magnitude_records", "sum"),
        )
        .reindex(AVAILABILITY_EFFECT_CATEGORIES)
        .reset_index()
    )
    overall["usable_percent"] = np.where(
        overall["records"].gt(0), 100 * overall["usable"] / overall["records"], np.nan
    )

    major_counts = availability.loc[
        availability["organism"].isin(MAJOR_ORGANISMS)
    ].pivot(
        index="organism",
        columns="lifespan_effect_category",
        values="usable_numeric_magnitude_records",
    ).reindex(MAJOR_ORGANISMS)
    major_counts = major_counts.reset_index()[
        ["organism", *AVAILABILITY_EFFECT_CATEGORIES]
    ]

    maxima = magnitude_summary.pivot(
        index="organism", columns="lifespan_effect_category", values="maximum"
    ).reindex(MAJOR_ORGANISMS).reset_index()
    maxima = maxima[["organism", *SUMMARY_EFFECT_CATEGORIES]]

    raw_numeric = data["numeric_lifespan_change"]
    usable_total = int(raw_numeric.notna().sum())
    negative_total = int(raw_numeric.lt(0).sum())
    zero_total = int(raw_numeric.eq(0).sum())
    positive_total = int(raw_numeric.gt(0).sum())
    non_numeric_nonmissing = int(
        data[MAGNITUDE_COLUMN].notna().sum() - raw_numeric.notna().sum()
    )

    dtype_lines = "\n".join(
        f"- `{column}`: `{dtype}`"
        for column, dtype in data[list(EXPECTED_COLUMNS)].dtypes.items()
    )
    column_lines = "\n".join(f"- `{column}`" for column in data.columns[:8])

    notes = f"""# Figure 6 lifespan-effect magnitude exploration

This is a descriptive audit and two diagnostic previews, not the final Figure 6. All calculations use GenAge records (rows); repeated genes were not deduplicated.

## Source audit

- Source: `GenAge Project/models_genes/genage_models.csv`
- Raw dataset size: **{len(data):,} rows x {len(data.columns) - 3} raw columns**
- Organism field: `{ORGANISM_COLUMN}`
- Lifespan-effect field: `{EFFECT_COLUMN}`
- Lifespan-change field: `{MAGNITUDE_COLUMN}`

Raw columns:

{column_lines}

Raw data types:

{dtype_lines}

## Data availability

{markdown_table(overall, {'records': '{:,.0f}', 'usable': '{:,.0f}', 'missing_magnitude': '{:,.0f}', 'usable_percent': '{:.2f}%'})}

Across the full dataset, **{usable_total:,} of {len(data):,} records ({100 * usable_total / len(data):.2f}%)** have a usable numeric magnitude. There were **{non_numeric_nonmissing}** non-missing values that failed numeric conversion; the remaining unavailable values were blank. `Increase and Decrease` and missing/other lifespan-effect records have no usable magnitudes.

The full organism-by-category audit is in `results/tables/fig6_lifespan_magnitude_data_availability.csv`.

## Major-organism sample sizes

Usable numeric magnitude records:

{markdown_table(major_counts, {category: '{:,.0f}' for category in AVAILABILITY_EFFECT_CATEGORIES})}

The five major organisms are retained in both previews. *Schizosaccharomyces pombe* has no usable magnitude values and is shown as `n=0`, rather than being silently excluded.

## Magnitude distributions and extremes

The requested mean, median, SD, quartiles, IQR, minimum, and maximum for each major organism and effect category are in `results/tables/fig6_lifespan_magnitude_summary.csv`. Category maxima are:

{markdown_table(maxima, {category: '{:,.6g}' for category in SUMMARY_EFFECT_CATEGORIES})}

The largest extension-side values are 1,000 (*C. elegans*, `age-1`), 600 (*S. cerevisiae*, `SIR2`), and 300 (*S. cerevisiae*, `SCH9`). The largest observed decrease magnitude is 85 (*S. cerevisiae*, `RPT5`). These records are retained. The top ten usable records per major organism/effect category are in `results/tables/fig6_lifespan_magnitude_extremes.csv`; groups with no usable values cannot contribute ranked records.

Because the range spans 4 to 1,000, Preview A uses a logarithmic magnitude axis and Preview B uses a symmetric-log signed axis. These transformations expose the dense central distributions while keeping every extreme observation visible.

## Is the magnitude field unsigned?

Yes, for this release it behaves as an unsigned field: among {usable_total:,} usable values, **{positive_total:,} are positive, {zero_total} are zero, and {negative_total} are negative**. Positive values occur in both `Increase` and `Decrease` records. The raw column was not changed.

For diagnosis only, `signed_lifespan_change` was reconstructed as positive absolute magnitude for `Increase`, negative absolute magnitude for `Decrease`, and undefined for `Increase and Decrease` or missing/other classifications.

Signed distribution by major organism:

{markdown_table(signed_summary, {'n': '{:,.0f}', 'mean': '{:.2f}', 'median': '{:.2f}', 'sd': '{:.2f}', 'q1': '{:.2f}', 'q3': '{:.2f}', 'iqr': '{:.2f}', 'minimum': '{:.2f}', 'maximum': '{:.2f}'})}

## Preview comparison

- **Preview A (extension only)** is the more direct match to Pedro's wording about "life extending effects" and avoids reconstructing direction. Its main limitations are the 1,000 and 600 extremes, uneven coverage, and the complete absence of usable *S. pombe* values.
- **Preview B (signed Increase versus Decrease)** is more informative for understanding the full direction-and-magnitude structure and exposes organism-specific mixtures of extension and shortening. However, its sign is reconstructed from a separate category, mixed-direction records are undefined, and it goes beyond a literal extension-only request.

Preview A is therefore the leading interpretation if Pedro's wording is taken literally; Preview B is the stronger diagnostic companion. Neither should be treated as final until the intended scope is confirmed.

## Ambiguities and review risks

- The column name does not encode units. The values appear to be percentage-like, but the unit should be confirmed from GenAge documentation or observation-level records before final axis wording.
- `avg lifespan change (max obsv)` appears to retain a maximum observed summary rather than all underlying observations; a gene-level row may therefore conceal intervention, sex, dose, tissue, or study-level variation.
- All 144 `Increase and Decrease` records lack a numeric magnitude, so their two directions cannot be reconstructed from this table.
- *S. pombe* has 61 records but zero usable magnitude values, preventing a five-species quantitative comparison without supplementary observation-level data.
- The extreme extension values materially affect linear-scale views and mean estimates. They should be checked against the underlying GenAge observations, not automatically removed.
- The records are not independent effect estimates in a formal meta-analytic sense; no inferential test or biological ranking is justified at this stage.

## Diagnostic figure contract

- Core question: determine whether available GenAge summary values support a defensible cross-organism comparison of reported lifespan-effect magnitude.
- Archetype: quantitative comparison.
- Backend: Python/matplotlib only.
- Evidence: Preview A directly compares extension magnitudes; Preview B tests whether reconstructing direction around zero adds useful information.
- Source data: the three Figure 6 CSV audit tables produced by the companion script.
- Reviewer risk: unsigned source values, uneven missingness, mixed-direction records without magnitudes, non-independent records, and a small number of extreme observations.
"""
    NOTES_FILE.parent.mkdir(parents=True, exist_ok=True)
    NOTES_FILE.write_text(notes, encoding="utf-8")


def write_tables(
    availability: pd.DataFrame,
    magnitude_summary: pd.DataFrame,
    extremes: pd.DataFrame,
) -> None:
    """Write the three requested Figure 6 audit tables."""

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    availability.to_csv(
        TABLE_DIR / "fig6_lifespan_magnitude_data_availability.csv",
        index=False,
        float_format="%.6g",
    )
    magnitude_summary.to_csv(
        TABLE_DIR / "fig6_lifespan_magnitude_summary.csv",
        index=False,
        float_format="%.6g",
    )
    extremes.to_csv(
        TABLE_DIR / "fig6_lifespan_magnitude_extremes.csv",
        index=False,
        float_format="%.6g",
    )


def print_results_summary(
    data: pd.DataFrame,
    availability: pd.DataFrame,
    signed_summary: pd.DataFrame,
) -> None:
    """Print key audit totals and the signed diagnostic distribution."""

    overall = availability.groupby("lifespan_effect_category", sort=False).agg(
        records=("category_records", "sum"),
        usable=("usable_numeric_magnitude_records", "sum"),
        missing=("missing_magnitude_records", "sum"),
    )
    overall["usable_percent"] = 100 * overall["usable"] / overall["records"]
    print("Step 3 -- lifespan-effect categories and magnitude availability")
    print(overall.to_string(float_format=lambda value: f"{value:.2f}"))
    print("\nStep 6 -- exploratory signed-lifespan-change distribution")
    print(signed_summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    raw = data["numeric_lifespan_change"]
    print(
        "\nUnsigned-field check: "
        f"positive={int(raw.gt(0).sum()):,}, zero={int(raw.eq(0).sum()):,}, "
        f"negative={int(raw.lt(0).sum()):,}."
    )


def main() -> None:
    """Run the complete Figure 6 exploratory audit and preview workflow."""

    data = load_and_validate_data()
    print_source_audit(data)

    availability = build_availability_table(data)
    magnitude_summary = build_magnitude_summary(data)
    extremes = build_extremes_table(data)
    signed_summary = build_signed_summary(data)

    write_tables(availability, magnitude_summary, extremes)
    configure_plot_style()
    plot_extension_only(data)
    plot_signed_effects(data)
    write_notes(data, availability, magnitude_summary, signed_summary)
    print_results_summary(data, availability, signed_summary)

    print("\nOutputs written:")
    print(f"  - 3 CSV audit tables: {TABLE_DIR}")
    print(f"  - 2 PNG and 2 SVG diagnostic previews: {FIGURE_DIR}")
    print(f"  - exploration notes: {NOTES_FILE}")
    print("Figures 1--5b were not read or written by this script.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Explore or finalize GenAge Figure 6 lifespan-effect plots."
    )
    parser.add_argument(
        "--final-plots-only",
        action="store_true",
        help="Export the finalized Figure 6 plots and note without writing CSV tables.",
    )
    parser.add_argument(
        "--final-main-only",
        action="store_true",
        help="Update only the main Figure 6 plot and its standalone caption.",
    )
    arguments = parser.parse_args()
    if arguments.final_main_only:
        export_final_main_figure_only()
    elif arguments.final_plots_only:
        export_final_figures_only()
    else:
        main()
