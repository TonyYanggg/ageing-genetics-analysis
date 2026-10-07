"""Create the Task 5 review figure from independently verified outputs."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# Editable text and publication-safe font defaults.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.7
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["legend.frameon"] = False

TASK = Path(__file__).resolve().parents[1]
RESULTS = TASK / "results"
FIGURES = TASK / "figures"
SOURCE = RESULTS / "figure_1_source_data"

BLUE = "#0F4D92"
BLUE_LIGHT = "#8FB4D9"
TEAL = "#42949E"
VIOLET = "#9A4D8E"
GOLD = "#C99324"
GREY = "#A8A8A8"
DARK = "#272727"
LIGHT = "#E8E8E8"

CLASS_ORDER = [
    "disease_or_clinical_condition",
    "phenotype_or_trait",
    "biological_process",
    "medical_procedure",
    "nonhuman_or_other",
]
CLASS_LABELS = {
    "disease_or_clinical_condition": "Disease/clinical condition",
    "phenotype_or_trait": "Phenotype/trait",
    "biological_process": "Biological process",
    "medical_procedure": "Medical procedure",
    "nonhuman_or_other": "Non-human/other",
}
CLASS_COLORS = {
    "disease_or_clinical_condition": BLUE,
    "phenotype_or_trait": TEAL,
    "biological_process": VIOLET,
    "medical_procedure": GREY,
    "nonhuman_or_other": GOLD,
}
DATATYPE_LABELS = {
    "literature": "Literature",
    "animal_model": "Animal model",
    "somatic_mutation": "Somatic mutation",
    "genetic_association": "Genetic association",
    "clinical": "Clinical",
    "rna_expression": "RNA expression",
    "genetic_literature": "Genetic literature",
    "affected_pathway": "Affected pathway",
}


def panel_label(ax, label):
    ax.text(-0.075, 1.04, label, transform=ax.transAxes, fontsize=8,
            fontweight="bold", ha="left", va="bottom")


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    verification = json.loads((RESULTS / "verification_report.json").read_text(encoding="utf-8"))
    assert verification["status"] == "PASS"
    collection = json.loads((RESULTS / "opentargets_collection_summary.json").read_text(encoding="utf-8"))
    assert collection["opentargets_release"]["dataVersion"] == "26.06"
    assert collection["direct_only"] is True

    genes = pd.read_csv(RESULTS / "gene_disease_summary.csv")
    assert len(genes) == 307 and genes["genage_id"].nunique() == 307
    associations = pd.read_csv(
        RESULTS / "opentargets_gene_disease.csv",
        usecols=["genage_id", "disease_id", "term_class", "overall_association_score"],
    )
    assert len(associations) == collection["unique_gene_disease_pairs"] == 422250
    datatypes = pd.read_csv(
        RESULTS / "opentargets_datatype_scores.csv",
        usecols=["genage_id", "disease_id", "datatype_id"],
    ).drop_duplicates()

    FIGURES.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)

    # Source Data b: complete scope and mutually exclusive term classes.
    class_counts = associations["term_class"].value_counts().reindex(CLASS_ORDER, fill_value=0)
    scope_rows = [
        {"metric": "GenAge genes", "category": "all", "count": 307, "percent": 100.0},
        {"metric": "Queryable Open Targets genes", "category": "all", "count": 306, "percent": 306 / 307 * 100},
        {"metric": "Unqueryable genes", "category": "MT-CO1", "count": 1, "percent": 1 / 307 * 100},
    ]
    scope_rows += [
        {"metric": "Unique direct relations", "category": key, "count": int(class_counts[key]),
         "percent": float(class_counts[key] / class_counts.sum() * 100)}
        for key in CLASS_ORDER
    ]
    write_csv(SOURCE / "figure_1b_scope_and_term_classes.csv", scope_rows,
              ["metric", "category", "count", "percent"])

    # Source Data a: all 307 genes, retaining the unqueryable state.
    gene_source_cols = [
        "genage_id", "genage_symbol", "current_symbol", "entrez_id", "ensembl_gene_id",
        "query_status", "direct_association_count", "disease_or_clinical_condition_count",
        "phenotype_or_trait_count", "biological_process_count", "medical_procedure_count",
        "nonhuman_or_other_count", "max_overall_association_score", "longevitymap_record_count",
        "longevitymap_association_labels",
    ]
    genes[gene_source_cols].to_csv(SOURCE / "figure_1a_gene_counts.csv", index=False)

    # Source Data c: compact exact-rank survival curve plus descriptive quantiles.
    scores = associations["overall_association_score"].to_numpy(float)
    assert np.all(np.isfinite(scores)) and np.all(scores >= 0)
    zero_score_count = int(np.sum(scores == 0))
    positive_scores = scores[scores > 0]
    sorted_scores = np.sort(positive_scores)
    sample_indices = np.unique(np.linspace(0, len(sorted_scores) - 1, 2000, dtype=int))
    survival = 100.0 * (len(sorted_scores) - sample_indices) / len(scores)
    score_source = pd.DataFrame({
        "overall_association_score": sorted_scores[sample_indices],
        "relations_at_or_above_percent": survival,
    })
    score_source.to_csv(SOURCE / "figure_1c_score_survival.csv", index=False)
    median_score = float(np.median(scores))
    p90_score = float(np.quantile(scores, 0.90))

    # Source Data d: evidence datatypes, counting unique pairs and genes.
    datatype_summary = (
        datatypes.groupby("datatype_id", as_index=False)
        .agg(unique_gene_term_pairs=("disease_id", "size"), genes=("genage_id", "nunique"))
        .sort_values("unique_gene_term_pairs", ascending=False)
    )
    datatype_summary["display_label"] = datatype_summary["datatype_id"].map(DATATYPE_LABELS).fillna(
        datatype_summary["datatype_id"]
    )
    datatype_summary.to_csv(SOURCE / "figure_1d_evidence_datatypes.csv", index=False)

    # Figure: double-column width, asymmetric quantitative grid.
    width_in = 183 / 25.4
    height_in = 155 / 25.4
    fig = plt.figure(figsize=(width_in, height_in), facecolor="white")
    gs = gridspec.GridSpec(
        3, 2, figure=fig, height_ratios=[1.45, 1.2, 1.35],
        width_ratios=[1.0, 1.18], hspace=0.62, wspace=0.42,
    )

    # a — hero panel: every queryable gene ranked by relation count.
    ax_b = fig.add_subplot(gs[0, :])
    ranked = genes.loc[genes["query_status"] == "queried_direct_associations"].copy()
    ranked = ranked.sort_values(["direct_association_count", "genage_symbol"], ascending=[False, True]).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    has_lm = ranked["longevitymap_record_count"] > 0
    ax_b.scatter(ranked.loc[~has_lm, "rank"], ranked.loc[~has_lm, "direct_association_count"],
                 s=9, color=GREY, alpha=0.72, linewidths=0, rasterized=False)
    ax_b.scatter(ranked.loc[has_lm, "rank"], ranked.loc[has_lm, "direct_association_count"],
                 s=11, color=BLUE, alpha=0.9, linewidths=0, rasterized=False)
    ax_b.set_yscale("log")
    ax_b.set_xlim(0, 307)
    ax_b.set_ylim(30, 9000)
    ax_b.set_xlabel("Rank among 306 queryable GenAge genes")
    ax_b.set_ylabel("Unique direct relations")
    ax_b.grid(axis="y", color=LIGHT, linewidth=0.6, which="both")
    ax_b.set_axisbelow(True)
    label_symbols = ["EGFR", "LMNA", "TNF", "VEGFA", "TP53", "GTF2H2"]
    offsets = {"EGFR": (3, 5), "LMNA": (10, -1), "TNF": (18, -9), "VEGFA": (28, 3), "TP53": (35, -8), "GTF2H2": (-4, 6)}
    for symbol in label_symbols:
        row = ranked.loc[ranked["genage_symbol"] == symbol].iloc[0]
        dx, dy = offsets[symbol]
        ax_b.annotate(symbol, (row["rank"], row["direct_association_count"]),
                      xytext=(dx, dy), textcoords="offset points", fontsize=6,
                      ha="left" if dx >= 0 else "right", va="bottom",
                      arrowprops={"arrowstyle": "-", "color": GREY, "lw": 0.5})
    ax_b.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgewidth=0,
               markersize=4, label=f"LongevityMap record(s), n={int(has_lm.sum())}"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY, markeredgewidth=0,
               markersize=4, label=f"No matched LongevityMap record, n={int((~has_lm).sum())}"),
    ], loc="lower left", ncol=2, fontsize=6, handletextpad=0.4, columnspacing=1.2)
    panel_label(ax_b, "a")

    # b — fixed universe and term-class composition.
    ax_a = fig.add_subplot(gs[1, 0])
    ax_a.set_axis_off()
    ax_a.text(0.00, 0.93, "Complete human catalogue", fontsize=7, fontweight="bold", transform=ax_a.transAxes)
    ax_a.text(0.00, 0.72, "307", fontsize=18, fontweight="bold", color=DARK, transform=ax_a.transAxes)
    ax_a.text(0.17, 0.76, "GenAge genes", fontsize=7, color=DARK, transform=ax_a.transAxes)
    ax_a.text(0.52, 0.72, "306", fontsize=18, fontweight="bold", color=BLUE, transform=ax_a.transAxes)
    ax_a.text(0.70, 0.76, "queried", fontsize=7, color=DARK, transform=ax_a.transAxes)
    ax_a.text(0.00, 0.58, "MT-CO1: no Ensembl target; association status unknown", fontsize=5.8,
              color="#B64342", transform=ax_a.transAxes)
    ax_a.text(0.00, 0.43, f"{len(associations):,} unique direct relations", fontsize=7,
              fontweight="bold", transform=ax_a.transAxes)
    left = 0.0
    for key in CLASS_ORDER:
        width = float(class_counts[key] / class_counts.sum())
        ax_a.add_patch(plt.Rectangle((left, 0.30), width, 0.08, transform=ax_a.transAxes,
                                     color=CLASS_COLORS[key], linewidth=0))
        left += width
    legend_y = [0.20, 0.11, 0.02]
    legend_items = list(CLASS_ORDER)
    for i, key in enumerate(legend_items):
        col = i % 2
        row_i = i // 2
        x = 0.00 if col == 0 else 0.52
        y = legend_y[row_i]
        ax_a.add_patch(plt.Rectangle((x, y), 0.025, 0.035, transform=ax_a.transAxes,
                                     color=CLASS_COLORS[key], linewidth=0))
        pct = class_counts[key] / class_counts.sum() * 100
        pct_text = f"{pct:.1f}%" if pct >= 1 else f"{pct:.2f}%"
        ax_a.text(x + 0.035, y - 0.001, f"{CLASS_LABELS[key]}  {pct_text}",
                  fontsize=5.2, va="bottom", transform=ax_a.transAxes)
    panel_label(ax_a, "b")

    # c — score survival distribution; score is explicitly not a probability.
    ax_c = fig.add_subplot(gs[1, 1])
    ax_c.plot(score_source["overall_association_score"], score_source["relations_at_or_above_percent"],
              color=TEAL, lw=1.5)
    ax_c.set_xscale("log")
    ax_c.set_xlim(max(positive_scores.min() * 0.72, 1e-5), 1.0)
    ax_c.set_ylim(0, 102)
    ax_c.set_xlabel("Open Targets overall score")
    ax_c.set_ylabel("Relations at or above score (%)")
    ax_c.grid(color=LIGHT, linewidth=0.6)
    for value, label in [(median_score, "median"), (p90_score, "90th percentile")]:
        y = 100 * np.mean(scores >= value)
        ax_c.axvline(value, color=GREY, lw=0.7, linestyle="--")
        ax_c.text(value, min(96, y + 5), f"{label}\n{value:.3g}", fontsize=5.5,
                  ha="left", va="bottom", color=DARK)
    ax_c.text(0.98, 0.03, "Score ≠ probability/effect size", transform=ax_c.transAxes,
              fontsize=5.3, ha="right", color="#B64342")
    ax_c.text(0.02, 0.96, f"{zero_score_count:,} relations score = 0\n(retained in denominator)",
              transform=ax_c.transAxes, fontsize=5.3, ha="left", va="top", color=DARK)
    panel_label(ax_c, "c")

    # d — heterogeneous evidence datatypes.
    ax_d = fig.add_subplot(gs[2, :])
    ds = datatype_summary.sort_values("unique_gene_term_pairs", ascending=True).reset_index(drop=True)
    colors = [BLUE_LIGHT if x not in {"literature", "animal_model"} else BLUE for x in ds["datatype_id"]]
    y = np.arange(len(ds))
    ax_d.barh(y, ds["unique_gene_term_pairs"], color=colors, height=0.66, edgecolor="none")
    ax_d.set_xscale("log")
    ax_d.set_yticks(y)
    ax_d.set_yticklabels(ds["display_label"])
    ax_d.set_xlabel("Unique gene–term pairs with datatype evidence")
    ax_d.grid(axis="x", color=LIGHT, linewidth=0.6, which="both")
    ax_d.set_axisbelow(True)
    for yi, row in ds.iterrows():
        ax_d.text(row["unique_gene_term_pairs"] * 1.10, yi,
                  f"{int(row['unique_gene_term_pairs']):,} pairs; {int(row['genes'])} genes",
                  fontsize=5.5, va="center", ha="left")
    ax_d.set_xlim(500, 450000)
    panel_label(ax_d, "d")

    fig.text(0.995, 0.006, "Open Targets Platform 26.06 · direct associations only",
             ha="right", va="bottom", fontsize=5.2, color="#606060")

    base = FIGURES / "task5_human_gene_disease_review"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    svg_text = base.with_suffix(".svg").read_text(encoding="utf-8")
    qa = {
        "figure": base.name,
        "backend": "Python/matplotlib",
        "final_size_mm": [183, 155],
        "exports": ["svg", "pdf", "tiff", "png"],
        "editable_svg_text_nodes": "<text" in svg_text,
        "source_data_files": sorted(p.name for p in SOURCE.glob("*.csv")),
        "gene_count": int(len(genes)),
        "queryable_gene_count": int((genes["query_status"] == "queried_direct_associations").sum()),
        "relation_count": int(len(associations)),
        "median_overall_score": median_score,
        "p90_overall_score": p90_score,
        "zero_score_relation_count": zero_score_count,
        "longevitymap_annotated_gene_count": int((genes["longevitymap_record_count"] > 0).sum()),
        "inferential_statistics": "none; complete resource-bounded descriptive catalogue",
    }
    (RESULTS / "figure_1_qa.json").write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qa), flush=True)


if __name__ == "__main__":
    main()
