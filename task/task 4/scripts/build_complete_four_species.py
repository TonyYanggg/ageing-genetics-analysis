"""Build the complete four-organism, direct-OMA-pair Task 4 review package."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


TASK = Path(__file__).resolve().parents[1]
SOURCE = TASK / "results/oma_pairwise/refined"
OUT = TASK / "results/complete_four_species"
FIG = TASK / "figures/complete_four_species"
SPECIES = (
    ("Saccharomyces cerevisiae", "Budding yeast"),
    ("Caenorhabditis elegans", "Nematode"),
    ("Drosophila melanogaster", "Fruit fly"),
    ("Mus musculus", "Mouse"),
)
NAMES = dict(SPECIES)
DIRECTIONS = ("Anti-Longevity", "Pro-Longevity", "discordant_or_mixed", "unresolved")
COLORS = ("#427F77", "#6687B1", "#B78363", "#BFC8D0")
EFFECTS = {"Increase", "Decrease", "Increase and Decrease"}


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    raw = read(TASK.parents[1] / "GenAge Project/models_genes/genage_models.csv")
    source = read(SOURCE / "longevity_pairs.csv")
    coverage = read(SOURCE / "coverage.csv")
    by_gene = defaultdict(list)
    for row in raw:
        by_gene[(row["organism"], row["entrez gene id"])].append(row)

    pairs = []
    genes = set()
    partners = defaultdict(set)
    for row in source:
        if row["species1"] not in NAMES or row["species2"] not in NAMES:
            continue
        assert row["both_lifespan"] == "True" and row["unit"] == "explicit_gene_pair"
        left = (row["species1"], row["entrez1"])
        right = (row["species2"], row["entrez2"])
        assert left < right
        assert any(r["lifespan effect"] in EFFECTS for r in by_gene[left])
        assert any(r["lifespan effect"] in EFFECTS for r in by_gene[right])
        genes.update((left, right))
        partners[left].add(right)
        partners[right].add(left)
        extended = dict(row)
        extended["lifespan_effects1"] = ";".join(sorted({r["lifespan effect"] for r in by_gene[left] if r["lifespan effect"] in EFFECTS}))
        extended["lifespan_effects2"] = ";".join(sorted({r["lifespan effect"] for r in by_gene[right] if r["lifespan effect"] in EFFECTS}))
        pairs.append(extended)
    assert len(pairs) == len({(r["species1"], r["entrez1"], r["species2"], r["entrez2"]) for r in pairs})
    pairs.sort(key=lambda r: (r["species1"], r["species2"], r["symbols1"], r["symbols2"], r["entrez1"], r["entrez2"]))
    write(OUT / "all_qualifying_pairs.csv", pairs, list(pairs[0]))

    order = [s for s, _ in SPECIES]
    counts = Counter(tuple(sorted((r["species1"], r["species2"]))) for r in pairs)
    directions = Counter((tuple(sorted((r["species1"], r["species2"]))), r["direction_agreement"]) for r in pairs)
    matrix = []
    for i, first in enumerate(order):
        for second in order[i + 1:]:
            key = tuple(sorted((first, second)))
            record = {"species1": first, "species2": second, "pairs": counts[key]}
            record.update({d: directions[(key, d)] for d in DIRECTIONS})
            assert sum(record[d] for d in DIRECTIONS) == record["pairs"]
            matrix.append(record)
    write(OUT / "species_pair_summary.csv", matrix, list(matrix[0]))

    coverage_rows = []
    for row in coverage:
        item = dict(row)
        item["analysis_role"] = "principal_four" if row["species"] in NAMES else "outside_main_comparison"
        coverage_rows.append(item)
    write(OUT / "all_nine_coverage.csv", coverage_rows, list(coverage_rows[0]))

    multiplicity = []
    for species, entrez in sorted(genes):
        other_species = {partner[0] for partner in partners[(species, entrez)]}
        multiplicity.append({"species": species, "entrez_id": entrez,
                             "symbols": ";".join(sorted({r["symbol"] for r in by_gene[(species, entrez)]})),
                             "genage_ids": ";".join(sorted({r["GenAge ID"] for r in by_gene[(species, entrez)]}, key=int)),
                             "partner_genes": len(partners[(species, entrez)]),
                             "partner_species": len(other_species)})
    write(OUT / "gene_pair_multiplicity.csv", multiplicity, list(multiplicity[0]))
    summary = {"source": "GenAge Build 21 and OMA All.May2026", "unit": "direct species+Entrez gene pair",
               "principal_species": order, "genage_records_principal": sum(int(r["records"]) for r in coverage if r["species"] in NAMES),
               "qualifying_pairs": len(pairs), "distinct_genes": len(genes),
               "genes_with_multiple_partners": sum(int(r["partner_genes"]) > 1 for r in multiplicity),
               "genes_with_partners_in_multiple_species": sum(int(r["partner_species"]) > 1 for r in multiplicity),
               "direction_counts": dict(Counter(r["direction_agreement"] for r in pairs)),
               "outside_main_comparison_records": sum(int(r["records"]) for r in coverage if r["species"] not in NAMES),
               "oma_unavailable_species": [r["species"] for r in coverage if r["source_status"] == "unavailable"]}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(183 / 25.4, 143 / 25.4), facecolor="white")
    fig.text(.06, .945, "Direct ortholog pairs with lifespan evidence", size=12, weight="bold", color="#203449")
    fig.text(.06, .905, f"Four major GenAge model organisms  |  {len(pairs)} pairs across {len(genes)} distinct genes",
             size=8, color="#4D6172")

    ax = fig.add_axes([.16, .55, .30, .27])
    rows = [next(r for r in coverage if r["species"] == s) for s in order]
    ypos = list(range(4))
    mapped = [int(r["mapped_genes"]) for r in rows]
    identified = [int(r["identified_genes"]) for r in rows]
    ax.barh(ypos, [100 * a / b for a, b in zip(mapped, identified)], color="#668CA0", height=.57)
    ax.set_yticks(ypos, [name for _, name in SPECIES])
    ax.invert_yaxis()
    ax.set_xlim(0, 123)
    ax.set_xticks([0, 50, 100], ["0", "50", "100%"])
    ax.set_title("a  Exact-ID mapping coverage", loc="left", size=9, weight="bold", pad=9)
    for y, a, b in zip(ypos, mapped, identified):
        ax.text(102, y, f"{a}/{b}", va="center", size=7, color="#203449")
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(length=0)

    bx = fig.add_axes([.55, .535, .39, .30])
    bx.set_xlim(0, 4)
    bx.set_ylim(4, 0)
    bx.set_xticks([i + .5 for i in range(4)], ["Yeast", "Worm", "Fly", "Mouse"])
    bx.set_yticks([i + .5 for i in range(4)], ["Yeast", "Worm", "Fly", "Mouse"])
    bx.xaxis.tick_top()
    bx.tick_params(length=0)
    bx.set_title("b  Direct pairs by organism", loc="left", size=9, weight="bold", pad=18)
    high = max(counts.values())
    for i, first in enumerate(order):
        for j, second in enumerate(order):
            if j <= i:
                continue
            count = counts[tuple(sorted((first, second)))]
            shade = plt.cm.Blues(.16 + .65 * count / high)
            bx.add_patch(Rectangle((j, i), 1, 1, facecolor=shade, edgecolor="white", linewidth=2))
            bx.text(j + .5, i + .5, str(count), ha="center", va="center", size=10, weight="bold",
                    color="white" if count > high * .55 else "#203449")
    for side in bx.spines.values():
        side.set_visible(False)

    cx = fig.add_axes([.31, .18, .61, .25])
    labels = [f"{NAMES[r['species1']]} – {NAMES[r['species2']]}" for r in matrix]
    for i, row in enumerate(matrix):
        start = 0
        for d, color in zip(DIRECTIONS, COLORS):
            value = row[d]
            cx.barh(i, value, left=start, height=.63, color=color)
            if value >= 6:
                cx.text(start + value / 2, i, str(value), ha="center", va="center", size=6.5,
                        color="white" if d != "unresolved" else "#203449")
            start += value
        cx.text(start + 1, i, str(row["pairs"]), va="center", size=7, color="#203449")
    cx.set_yticks(range(6), labels)
    cx.invert_yaxis()
    cx.set_xlim(0, max(r["pairs"] for r in matrix) + 13)
    cx.set_xlabel("Direct pairs (each pair counted once)")
    cx.set_title("c  GenAge direction annotation by organism pair", loc="left", size=9, weight="bold", pad=10)
    for side in ("top", "right", "left"):
        cx.spines[side].set_visible(False)
    cx.tick_params(axis="y", length=0)
    fig.legend([Rectangle((0, 0), 1, 1, color=c) for c in COLORS],
               ["Anti/Anti", "Pro/Pro", "Mixed/opposed", "Unresolved"],
               loc="lower center", bbox_to_anchor=(.55, .085), ncol=4, frameon=False, fontsize=7)
    fig.text(.06, .035, "Pairs overlap in genes; categories are catalogue labels. Missing pairs are not biological negatives.",
             size=6.7, color="#4D6172")
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf", "tiff"):
        opts = {"dpi": 600 if ext == "tiff" else 220}
        if ext == "tiff":
            opts["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(FIG / f"complete_ortholog_overview.{ext}", **opts)
    plt.close(fig)
    print(f"Built {len(pairs)} pairs, {len(genes)} genes, six organism comparisons")


if __name__ == "__main__":
    main()
