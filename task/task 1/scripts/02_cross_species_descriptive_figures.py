"""Create cross-species GenAge category tables and exploratory figures.

The raw model-organism CSV is treated as immutable. All nine organisms are
retained in the source tables, while figures compare the five major organisms.
No records are deduplicated and no statistical tests are performed.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter, StrMethodFormatter
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

LONGEVITY_CATEGORIES = [
    "Anti-Longevity",
    "Pro-Longevity",
    "Necessary for fitness",
    "Unannotated",
    "Unclear",
]
LIFESPAN_CATEGORIES = [
    "Increase",
    "Decrease",
    "Increase and Decrease",
    "Missing",
]

EXPECTED_LONGEVITY_TOTALS = {
    "Anti-Longevity": 1_101,
    "Pro-Longevity": 546,
    "Necessary for fitness": 497,
    "Unannotated": 34,
    "Unclear": 27,
}
EXPECTED_LIFESPAN_TOTALS = {
    "Increase": 1_155,
    "Decrease": 892,
    "Increase and Decrease": 144,
    "Missing": 14,
}

# Fixed category palettes used consistently between absolute and percentage plots.
LONGEVITY_COLORS = {
    "Anti-Longevity": "#4E79A7",
    "Pro-Longevity": "#59A14F",
    "Necessary for fitness": "#F28E2B",
    "Unannotated": "#BAB0AC",
    "Unclear": "#B07AA1",
}
LIFESPAN_COLORS = {
    "Increase": "#59A14F",
    "Decrease": "#E15759",
    "Increase and Decrease": "#4E79A7",
    "Missing": "#BAB0AC",
}


def configure_plot_style() -> None:
    """Apply the same restrained style as the first exploratory figures."""

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
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def load_and_validate_data() -> pd.DataFrame:
    """Read the raw CSV and validate the structure and known record totals."""

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

    observed_organisms = data["organism"].value_counts().to_dict()
    if observed_organisms != EXPECTED_ORGANISM_COUNTS:
        raise ValueError(
            "Organism record counts differ from the validated inventory. "
            f"Expected {EXPECTED_ORGANISM_COUNTS}; observed {observed_organisms}."
        )

    return data


def make_crosstab(
    data: pd.DataFrame,
    source_column: str,
    categories: list[str],
    missing_label: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return full-organism count and within-organism percentage tables."""

    values = data[source_column]
    if missing_label is not None:
        values = values.fillna(missing_label)

    unexpected = set(values.dropna().unique()).difference(categories)
    if unexpected:
        raise ValueError(
            f"Unexpected values in {source_column!r}: {sorted(unexpected)}"
        )

    counts = pd.crosstab(data["organism"], values, dropna=False)
    counts = counts.reindex(index=FULL_ORGANISM_ORDER, columns=categories, fill_value=0)
    counts.index.name = "organism"
    counts["record_total"] = counts[categories].sum(axis=1)

    expected_totals = pd.Series(EXPECTED_ORGANISM_COUNTS).reindex(FULL_ORGANISM_ORDER)
    if not counts["record_total"].equals(expected_totals.rename("record_total")):
        raise ValueError(
            f"{source_column} row totals do not match the organism record totals."
        )

    unrounded_percent = counts[categories].div(counts["record_total"], axis=0) * 100
    if not unrounded_percent.sum(axis=1).sub(100).abs().lt(1e-10).all():
        raise ValueError(f"{source_column} percentages do not sum to 100%.")

    percentages = unrounded_percent.round(4)
    percentages["percentage_total"] = percentages[categories].sum(axis=1).round(4)
    if not percentages["percentage_total"].sub(100).abs().le(0.01).all():
        raise ValueError(
            f"Rounded {source_column} percentages are not approximately 100%."
        )

    return counts.reset_index(), percentages.reset_index()


def validate_global_category_totals(
    counts: pd.DataFrame,
    categories: list[str],
    expected: dict[str, int],
    label: str,
) -> None:
    """Validate complete-dataset category totals before writing outputs."""

    observed = {category: int(counts[category].sum()) for category in categories}
    if observed != expected:
        raise ValueError(
            f"{label} totals differ from the validated inventory. "
            f"Expected {expected}; observed {observed}."
        )


def major_rows(table: pd.DataFrame) -> pd.DataFrame:
    """Select and order the five major organisms for plotting."""

    return table.set_index("organism").loc[MAJOR_ORGANISMS].reset_index()


def style_species_axis(ax: plt.Axes) -> None:
    """Apply shared species-axis styling."""

    ax.invert_yaxis()
    ax.set_ylabel("")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#777777")
    for label in ax.get_yticklabels():
        label.set_fontstyle("italic")


def add_legend(
    fig: plt.Figure,
    categories: list[str],
    colors: dict[str, str],
) -> None:
    """Add a consistent category legend below the plotting area."""

    handles = [
        Patch(facecolor=colors[category], edgecolor="none", label=category)
        for category in categories
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=len(categories),
        fontsize=8.8,
        handlelength=1.2,
        columnspacing=1.2,
    )


def save_figure(fig: plt.Figure, basename: str) -> None:
    """Save a 300-dpi PNG and editable-text SVG."""

    fig.savefig(FIGURE_DIR / f"{basename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_DIR / f"{basename}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_stacked_counts(
    table: pd.DataFrame,
    categories: list[str],
    colors: dict[str, str],
    title: str,
    basename: str,
) -> None:
    """Plot absolute cross-species category counts as stacked bars."""

    plot_data = major_rows(table)
    fig, ax = plt.subplots(figsize=(8.8, 5.1))
    left = pd.Series(0.0, index=plot_data.index)

    for category in categories:
        values = plot_data[category].astype(float)
        ax.barh(
            plot_data["organism"],
            values,
            left=left,
            color=colors[category],
            edgecolor="white",
            linewidth=0.5,
        )
        left = left + values

    maximum = int(plot_data["record_total"].max())
    ax.set_xlim(0, maximum * 1.13)
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    ax.xaxis.grid(True, color="#D9D9D9", linewidth=0.7)
    ax.set_xlabel("Number of GenAge records", labelpad=8)
    ax.set_title(title, loc="left", pad=12)
    style_species_axis(ax)

    for y_position, total in enumerate(plot_data["record_total"]):
        ax.text(
            float(total) + maximum * 0.012,
            y_position,
            f"n={int(total):,}",
            va="center",
            ha="left",
            fontsize=9,
        )

    add_legend(fig, categories, colors)
    fig.subplots_adjust(left=0.34, right=0.97, top=0.86, bottom=0.25)
    save_figure(fig, basename)


def plot_stacked_percent(
    table: pd.DataFrame,
    count_table: pd.DataFrame,
    categories: list[str],
    colors: dict[str, str],
    title: str,
    basename: str,
) -> None:
    """Plot within-organism composition as 100% stacked horizontal bars."""

    plot_data = major_rows(table)
    plot_counts = major_rows(count_table)
    fig, ax = plt.subplots(figsize=(8.8, 5.1))
    left = pd.Series(0.0, index=plot_data.index)

    dark_text_categories = {"Necessary for fitness", "Unannotated", "Missing"}
    for category in categories:
        values = plot_data[category].astype(float)
        ax.barh(
            plot_data["organism"],
            values,
            left=left,
            color=colors[category],
            edgecolor="white",
            linewidth=0.5,
        )
        for y_position, (start, value) in enumerate(zip(left, values, strict=True)):
            if value >= 8:
                ax.text(
                    start + value / 2,
                    y_position,
                    f"{value:.1f}%",
                    va="center",
                    ha="center",
                    fontsize=8.5,
                    color="#202020" if category in dark_text_categories else "white",
                )
        left = left + values

    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax.xaxis.grid(True, color="#D9D9D9", linewidth=0.7)
    ax.set_xlabel("Percentage of GenAge records within organism", labelpad=8)
    ax.set_title(title, loc="left", pad=12)
    style_species_axis(ax)

    for y_position, total in enumerate(plot_counts["record_total"]):
        ax.text(
            101.2,
            y_position,
            f"n={int(total):,}",
            va="center",
            ha="left",
            fontsize=9,
            clip_on=False,
        )

    add_legend(fig, categories, colors)
    fig.subplots_adjust(left=0.34, right=0.92, top=0.86, bottom=0.25)
    save_figure(fig, basename)


def write_outputs(
    longevity_counts: pd.DataFrame,
    longevity_percent: pd.DataFrame,
    lifespan_counts: pd.DataFrame,
    lifespan_percent: pd.DataFrame,
) -> None:
    """Write source tables and the four requested figure pairs."""

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    longevity_counts.to_csv(
        TABLE_DIR / "longevity_influence_by_organism_counts.csv", index=False
    )
    longevity_percent.to_csv(
        TABLE_DIR / "longevity_influence_by_organism_percent.csv", index=False
    )
    lifespan_counts.to_csv(
        TABLE_DIR / "lifespan_effect_by_organism_counts.csv", index=False
    )
    lifespan_percent.to_csv(
        TABLE_DIR / "lifespan_effect_by_organism_percent.csv", index=False
    )

    configure_plot_style()
    plot_stacked_counts(
        longevity_counts,
        LONGEVITY_CATEGORIES,
        LONGEVITY_COLORS,
        "Longevity influence classifications by model organism",
        "fig4a_longevity_influence_counts_by_species",
    )
    plot_stacked_percent(
        longevity_percent,
        longevity_counts,
        LONGEVITY_CATEGORIES,
        LONGEVITY_COLORS,
        "Composition of longevity influence classifications across model organisms",
        "fig4b_longevity_influence_percent_by_species",
    )
    plot_stacked_counts(
        lifespan_counts,
        LIFESPAN_CATEGORIES,
        LIFESPAN_COLORS,
        "Reported lifespan-effect directions by model organism",
        "fig5a_lifespan_effect_counts_by_species",
    )
    plot_stacked_percent(
        lifespan_percent,
        lifespan_counts,
        LIFESPAN_CATEGORIES,
        LIFESPAN_COLORS,
        "Composition of reported lifespan-effect directions across model organisms",
        "fig5b_lifespan_effect_percent_by_species",
    )


def main() -> None:
    """Run the full cross-species descriptive workflow."""

    data = load_and_validate_data()
    longevity_counts, longevity_percent = make_crosstab(
        data,
        source_column="longevity influence",
        categories=LONGEVITY_CATEGORIES,
    )
    lifespan_counts, lifespan_percent = make_crosstab(
        data,
        source_column="lifespan effect",
        categories=LIFESPAN_CATEGORIES,
        missing_label="Missing",
    )

    validate_global_category_totals(
        longevity_counts,
        LONGEVITY_CATEGORIES,
        EXPECTED_LONGEVITY_TOTALS,
        "Longevity influence",
    )
    validate_global_category_totals(
        lifespan_counts,
        LIFESPAN_CATEGORIES,
        EXPECTED_LIFESPAN_TOTALS,
        "Lifespan effect",
    )

    write_outputs(
        longevity_counts,
        longevity_percent,
        lifespan_counts,
        lifespan_percent,
    )

    print(f"Validated {len(data):,} GenAge model-organism records.")
    print("All organism count totals and within-organism percentages matched.")
    print(f"Wrote 4 cross-tabulation tables to: {TABLE_DIR}")
    print(f"Wrote 4 PNG and 4 SVG figures to: {FIGURE_DIR}")


if __name__ == "__main__":
    main()
