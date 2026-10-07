"""Create basic descriptive tables and figures for GenAge model records.

This script reads the immutable GenAge model-organism CSV and writes only
derived tables and figures under ``results/`` and ``figures/``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import pandas as pd


TASK_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = (
    PROJECT_ROOT / "GenAge Project" / "models_genes" / "genage_models.csv"
)
TABLE_DIR = TASK_ROOT / "results" / "tables"
FIGURE_DIR = TASK_ROOT / "figures" / "exploratory"

EXPECTED_COLUMNS = {
    "GenAge ID",
    "symbol",
    "name",
    "organism",
    "entrez gene id",
    "avg lifespan change (max obsv)",
    "lifespan effect",
    "longevity influence",
}
EXPECTED_TOTAL_RECORDS = 2_205
EXPECTED_ORGANISM_COUNTS = {
    "Saccharomyces cerevisiae": 911,
    "Caenorhabditis elegans": 889,
    "Drosophila melanogaster": 202,
    "Mus musculus": 136,
    "Schizosaccharomyces pombe": 61,
    "Podospora anserina": 3,
    "Mesocricetus auratus": 1,
    "Danio rerio": 1,
    "Caenorhabditis briggsae": 1,
}
EXPECTED_LONGEVITY_COUNTS = {
    "Anti-Longevity": 1_101,
    "Pro-Longevity": 546,
    "Necessary for fitness": 497,
    "Unannotated": 34,
    "Unclear": 27,
}
EXPECTED_LIFESPAN_EFFECT_COUNTS = {
    "Increase": 1_155,
    "Decrease": 892,
    "Increase and Decrease": 144,
    "Missing": 14,
}

BAR_COLOR = "#4C78A8"


def configure_plot_style() -> None:
    """Apply a restrained, consistent exploratory-science style."""

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
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def assert_counts(actual: pd.Series, expected: dict[str, int], label: str) -> None:
    """Stop if a categorical count differs from the validated inventory."""

    actual_dict = {str(key): int(value) for key, value in actual.items()}
    if actual_dict != expected:
        raise ValueError(
            f"{label} counts differ from the validated inventory. "
            f"Expected {expected}; observed {actual_dict}."
        )


def load_and_validate_data() -> pd.DataFrame:
    """Read the raw CSV and validate its required structure and totals."""

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
    if data["organism"].isna().any():
        raise ValueError("The organism column contains missing values.")
    if data["longevity influence"].isna().any():
        raise ValueError("The longevity influence column contains missing values.")

    return data


def build_summary_tables(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate record/gene counts and the three figure source tables."""

    organism_gene_counts = (
        data.groupby("organism", sort=False)
        .agg(
            record_count=("GenAge ID", "size"),
            unique_gene_symbols=("symbol", "nunique"),
            unique_entrez_gene_ids=("entrez gene id", "nunique"),
        )
        .reset_index()
        .sort_values(["record_count", "organism"], ascending=[False, True])
        .reset_index(drop=True)
    )
    count_columns = [
        "record_count",
        "unique_gene_symbols",
        "unique_entrez_gene_ids",
    ]
    organism_gene_counts[count_columns] = organism_gene_counts[count_columns].astype(int)

    fig1_table = organism_gene_counts[["organism", "record_count"]].copy()
    fig1_table["percentage_of_all_records"] = (
        100 * fig1_table["record_count"] / len(data)
    ).round(2)

    longevity_counts = data["longevity influence"].value_counts()
    fig2_table = longevity_counts.rename_axis("longevity_influence").reset_index(
        name="record_count"
    )
    fig2_table["percentage_of_all_records"] = (
        100 * fig2_table["record_count"] / len(data)
    ).round(2)

    lifespan_effect = data["lifespan effect"].fillna("Missing")
    lifespan_counts = lifespan_effect.value_counts()
    fig3_table = lifespan_counts.rename_axis("lifespan_effect").reset_index(
        name="record_count"
    )
    fig3_table["percentage_of_all_records"] = (
        100 * fig3_table["record_count"] / len(data)
    ).round(2)

    assert_counts(
        organism_gene_counts.set_index("organism")["record_count"],
        EXPECTED_ORGANISM_COUNTS,
        "Organism",
    )
    assert_counts(longevity_counts, EXPECTED_LONGEVITY_COUNTS, "Longevity influence")
    assert_counts(
        lifespan_counts,
        EXPECTED_LIFESPAN_EFFECT_COUNTS,
        "Lifespan effect",
    )

    for table_name, table in {
        "organism record/gene": organism_gene_counts,
        "Figure 1": fig1_table,
        "Figure 2": fig2_table,
        "Figure 3": fig3_table,
    }.items():
        if int(table["record_count"].sum()) != EXPECTED_TOTAL_RECORDS:
            raise ValueError(f"{table_name} table does not sum to 2,205 records.")

    return organism_gene_counts, fig1_table, fig2_table, fig3_table


def style_bar_axis(ax: plt.Axes, maximum: int) -> None:
    """Apply shared formatting to a horizontal count bar chart."""

    ax.set_xlim(0, maximum * 1.22)
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax.xaxis.grid(True, color="#D9D9D9", linewidth=0.7)
    ax.yaxis.grid(False)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#777777")


def save_figure(fig: plt.Figure, basename: str) -> None:
    """Save a raster preview and an editable vector version."""

    fig.savefig(FIGURE_DIR / f"{basename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{basename}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_records_by_species(table: pd.DataFrame) -> None:
    """Create Figure 1: GenAge records by species."""

    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    bars = ax.barh(table["organism"], table["record_count"], color=BAR_COLOR)
    ax.invert_yaxis()
    ax.set_title("GenAge model-organism records by species", loc="left", pad=12)
    ax.set_xlabel("Number of GenAge records", labelpad=8)
    ax.set_ylabel("")
    style_bar_axis(ax, int(table["record_count"].max()))

    for label in ax.get_yticklabels():
        label.set_fontstyle("italic")
    for bar, count in zip(bars, table["record_count"], strict=True):
        ax.text(
            bar.get_width() + table["record_count"].max() * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{int(count):,}",
            va="center",
            ha="left",
            fontsize=9,
        )

    fig.subplots_adjust(left=0.34, right=0.97, top=0.88, bottom=0.14)
    save_figure(fig, "fig1_records_by_species")


def plot_categorical_counts(
    table: pd.DataFrame,
    category_column: str,
    title: str,
    basename: str,
    figsize: tuple[float, float],
) -> None:
    """Create a horizontal count chart with count and percentage labels."""

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.barh(table[category_column], table["record_count"], color=BAR_COLOR)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", pad=12)
    ax.set_xlabel("Number of GenAge records", labelpad=8)
    ax.set_ylabel("")
    style_bar_axis(ax, int(table["record_count"].max()))

    maximum = table["record_count"].max()
    for bar, count, percentage in zip(
        bars,
        table["record_count"],
        table["percentage_of_all_records"],
        strict=True,
    ):
        ax.text(
            bar.get_width() + maximum * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{int(count):,} ({percentage:.2f}%)",
            va="center",
            ha="left",
            fontsize=9,
        )

    fig.subplots_adjust(left=0.30, right=0.97, top=0.84, bottom=0.18)
    save_figure(fig, basename)


def write_outputs(
    organism_gene_counts: pd.DataFrame,
    fig1_table: pd.DataFrame,
    fig2_table: pd.DataFrame,
    fig3_table: pd.DataFrame,
) -> None:
    """Write all derived source-data tables and figures."""

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    organism_gene_counts.to_csv(
        TABLE_DIR / "organism_record_gene_counts.csv", index=False
    )
    fig1_table.to_csv(TABLE_DIR / "fig1_records_by_species.csv", index=False)
    fig2_table.to_csv(TABLE_DIR / "fig2_longevity_influence.csv", index=False)
    fig3_table.to_csv(TABLE_DIR / "fig3_lifespan_effect.csv", index=False)

    configure_plot_style()
    plot_records_by_species(fig1_table)
    plot_categorical_counts(
        fig2_table,
        category_column="longevity_influence",
        title="Longevity influence classifications in GenAge model-organism records",
        basename="fig2_longevity_influence",
        figsize=(8.4, 4.2),
    )
    plot_categorical_counts(
        fig3_table,
        category_column="lifespan_effect",
        title="Reported lifespan-effect directions in GenAge model-organism records",
        basename="fig3_lifespan_effect",
        figsize=(8.4, 3.8),
    )


def main() -> None:
    """Run the complete descriptive-table and figure workflow."""

    data = load_and_validate_data()
    tables = build_summary_tables(data)
    write_outputs(*tables)

    print(f"Validated {len(data):,} GenAge model-organism records.")
    print(f"Wrote 4 CSV tables to: {TABLE_DIR}")
    print(f"Wrote 3 PNG and 3 SVG figures to: {FIGURE_DIR}")
    print("All organism, longevity-influence, and lifespan-effect totals matched.")


if __name__ == "__main__":
    main()
