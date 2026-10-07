"""Expanded biological illustration with a reproducible direct-pair candidate audit.

Contract: six selected families illustrate recurring lifespan evidence and context
dependence. This evidence grid is not a set of source-defined OMA orthogroups.
Python only; 183 x 247 mm; editable SVG/PDF, PNG and 600-dpi TIFF.
"""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from plot_focused_examples import DISPLAY, SPECIES, NAMES, read_rows, write_rows

TASK = Path(__file__).resolve().parents[1]
OUT = TASK / "results/expanded"
FIG = TASK / "figures/expanded"
ORDER = ["TOR", "S6 kinase", "Insulin/IGF receptor", "AMPK", "ATG7", "SOD2"]
EDGE_NOTES = {
    "TOR": "All six species pairs; yeast relationships are one-to-many.",
    "S6 kinase": "Worm–fly and fly–mouse; no direct worm–mouse pair in the frozen selection.",
    "Insulin/IGF receptor": "Worm–fly and worm–mouse; Igf1r is one of two mouse receptor paralogs.",
    "AMPK": "Yeast–fly and worm–fly; no direct yeast–worm pair in the frozen selection.",
    "ATG7": "All three pairs among yeast, worm and fly; each is one-to-one.",
    "SOD2": "All six species pairs; worm relationships are one-to-many.",
}


def prepare():
    original = read_rows(TASK / "notes/focused_experiment_evidence.tsv", "\t")
    for row in original:
        _, intervention, outcome, _ = DISPLAY[row["genage_id"]]
        row.update(display_intervention=intervention, display_outcome=outcome)
    evidence = original + read_rows(TASK / "notes/expanded_additional_evidence.tsv", "\t")
    raw = {r["GenAge ID"]: r for r in read_rows(TASK.parents[1] / "GenAge Project/models_genes/genage_models.csv")}
    pairs = read_rows(TASK / "results/complete_four_species/all_qualifying_pairs.csv")
    groups = {name: {r["genage_id"] for r in evidence if r["example"] == name} for name in ORDER}
    for row in evidence:
        gene = raw[row["genage_id"]]
        assert row["species"] == gene["organism"] and row["gene"] == gene["symbol"]
        assert row["genage_category"] == gene["longevity influence"]
        row["genage_lifespan_effect"] = gene["lifespan effect"]
        row["primary_url"] = f"https://pubmed.ncbi.nlm.nih.gov/{row['primary_pmid']}/"
    selected = [{"example": name, **p} for name, ids in groups.items() for p in pairs
                if p["genage_ids1"] in ids and p["genage_ids2"] in ids]
    counts = Counter(p["example"] for p in selected)
    assert counts == {"TOR": 6, "S6 kinase": 2, "Insulin/IGF receptor": 2,
                      "AMPK": 2, "ATG7": 3, "SOD2": 6}
    write_rows(OUT / "displayed_evidence.csv", evidence, list(evidence[0]))
    write_rows(OUT / "selected_oma_pairs.csv", selected, list(selected[0]))

    # Connectivity only enumerates candidates; it never adds an orthology edge.
    graph = defaultdict(set)
    for p in pairs:
        a, b = p["genage_ids1"], p["genage_ids2"]
        graph[a].add(b)
        graph[b].add(a)
    seen, components = set(), []
    node_order = lambda node: tuple(int(part) for part in node.split(";"))
    for start in sorted(graph, key=node_order):
        if start in seen:
            continue
        stack, component = [start], set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(graph[node] - component)
        seen.update(component)
        components.append(component)
    audit = []
    for i, component in enumerate(components, 1):
        organisms = {raw[n.split(";")[0]]["organism"] for n in component}
        selected_names = [name for name, ids in groups.items() if ids <= component]
        audit.append(dict(candidate_id=f"C{i:03d}",
                          gene_ids="|".join(sorted(component, key=node_order)),
                          genes="|".join(raw[n.split(";")[0]]["symbol"] for n in sorted(component, key=node_order)),
                          species_count=len(organisms), gene_count=len(component),
                          direct_pair_count=sum(p["genage_ids1"] in component and p["genage_ids2"] in component for p in pairs),
                          selected_example=";".join(selected_names),
                          decision=("Selected; see selection rationale" if selected_names else
                                    "Eligible breadth; not selected or fully paper-reviewed" if len(organisms) >= 3 else
                                    "Two-species candidate; below preferred breadth"),
                          unit="connected candidate network; not an OMA orthogroup"))
    write_rows(OUT / "candidate_screen.csv", audit, list(audit[0]))
    return evidence, counts, audit


def draw(evidence, counts):
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(183 / 25.4, 247 / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 183), ylim=(0, 247))
    ax.axis("off")
    ink, muted, border = "#173047", "#4E6474", "#D5E0E7"
    checked = []

    def box(bounds, fill="white"):
        ax.add_patch(Rectangle(bounds[:2], bounds[2], bounds[3], facecolor=fill,
                               edgecolor=border, linewidth=.55))

    def text(label, x, y, bounds, size=7, color=ink, weight="normal"):
        artist = ax.text(x, y, label, fontsize=size, color=color, weight=weight, va="top")
        checked.append((artist, bounds))

    title = (5, 226, 173, 17)
    box(title, "#F7FAFB")
    text("Longevity-associated orthologs across major model organisms", 8, 240, title, 10.5, weight="bold")
    text("Six selected families • Primary lifespan experiments • Explicit OMA relationships", 8, 234, title, 7)
    text("Arrows show the reported effect of the stated intervention; they are not GenAge Pro/Anti labels.", 8, 230, title, 6.5, muted)
    label_width, cell_width, grid_x = 34, 34.75, 5
    for col, name in enumerate(NAMES):
        x = grid_x + label_width + col * cell_width
        bounds = (x, 216, cell_width - 1, 8)
        box(bounds, "#E8EFF3")
        text(name, x + 1.7, 221.7, bounds, 7.1, weight="bold")
    lookup = {(r["example"], r["species"]): r for r in evidence}
    fills = ["#E9F3F4", "#F1F2EC", "#EDF0F8", "#E9F3F4", "#F1F2EC", "#EDF0F8"]
    for index, name in enumerate(ORDER):
        top = 214 - index * 28
        bottom = top - 21
        bounds = (5, bottom, label_width - 1, 21)
        box(bounds, fills[index])
        label = "Insulin/IGF\nreceptor" if name == "Insulin/IGF receptor" else name
        text(label, 7, top - 3, bounds, 8.1, weight="bold")
        text(f"{counts[name]} direct OMA pairs", 7, bottom + 5, bounds, 6.1, muted)
        for col, species in enumerate(SPECIES):
            x = grid_x + label_width + col * cell_width
            cell = (x, bottom, cell_width - 1, 21)
            row = lookup.get((name, species))
            box(cell, "white" if row else "#F5F6F7")
            if not row:
                text("Not selected", x + 2, top - 9, cell, 7, "#82919A")
                continue
            outcome = row["display_outcome"]
            color = "#147268" if "↑" in outcome else "#995B28"
            text(row["gene"], x + 2, top - 2.3, cell, 8.2, weight="bold")
            text(row["display_intervention"], x + 2, top - 7.7, cell, 6.3, muted)
            text(outcome, x + 2, top - 12, cell, 6.5, color, "bold")
            text("PMID " + row["primary_pmid"], x + 2, top - 17, cell, 5.8, muted)
        edge = (5, bottom - 5.5, 173, 5)
        box(edge, "#F7FAFB")
        text("Direct edges: " + EDGE_NOTES[name], 7, bottom - 1.4, edge, 6.3, muted)
    footer = (5, 5, 173, 39)
    box(footer, "#FAFBFC")
    lines = [
        "* Mouse S6K1/Igf1r benefits shown are female-specific; corresponding male increases were not significant.",
        "† The original mouse Igf1r increase was not reproduced in a later cohort (PMID 22132081).",
        "‡ Worm ATG7: reduced longevity benefit in dietary-restricted eat-2 mutants; timing and background matter.",
        "§ Yeast SOD2: loss of sch9Δ chronological lifespan extension. Mouse SOD2: connective-tissue deficiency.",
        "AMPK and SOD2 illustrate context-dependent results. Yeast replicative and chronological life differ.",
        "TOR/S6K share a pathway. Rows are selected examples, not independent mechanisms or OMA orthogroups.",
        "Empty cells mean not selected, not absence of an ortholog or effect. See the legend for missing direct edges.",
        "Sources: GenAge Build 21, OMA All.May2026 and cited primary studies; full catalogue: 206 overlapping pairs.",
    ]
    for i, line in enumerate(lines):
        text(line, 7, 41 - i * 4.4, footer, 6.25, muted)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist, (x, y, w, h) in checked:
        bb = artist.get_window_extent(renderer)
        lo, hi = ax.transData.transform((x, y)), ax.transData.transform((x + w, y + h))
        assert bb.x0 >= lo[0] - .5 and bb.y0 >= lo[1] - .5 and bb.x1 <= hi[0] + .5 and bb.y1 <= hi[1] + .5, artist.get_text()
    for i, (a, _) in enumerate(checked):
        first = a.get_window_extent(renderer)
        for b, _ in checked[i + 1:]:
            second = b.get_window_extent(renderer)
            assert min(first.x1, second.x1) - max(first.x0, second.x0) <= .5 or min(first.y1, second.y1) - max(first.y0, second.y0) <= .5, (a.get_text(), b.get_text())
    for suffix in ("png", "svg", "pdf", "tiff"):
        opts = {"dpi": 600 if suffix == "tiff" else 220}
        if suffix == "tiff":
            opts["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(FIG / f"expanded_ortholog_examples.{suffix}", **opts)
    plt.close(fig)
    return len(checked)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    evidence, counts, audit = prepare()
    text_count = draw(evidence, counts)
    sources = [TASK / "results/complete_four_species/all_qualifying_pairs.csv",
               TASK / "notes/focused_experiment_evidence.tsv", TASK / "notes/expanded_additional_evidence.tsv"]
    summary = dict(families=len(ORDER), displayed_genes=len(evidence), direct_pairs=sum(counts.values()),
                   pairs_by_family=dict(counts), candidate_networks=len(audit),
                   candidates_with_at_least_three_species=sum(r["species_count"] >= 3 for r in audit),
                   text_elements_checked=text_count, final_mm=[183, 247],
                   source_sha256={str(p.relative_to(TASK)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
