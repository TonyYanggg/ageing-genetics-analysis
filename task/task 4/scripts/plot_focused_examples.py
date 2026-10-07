"""Pedro-facing, selected model-organism examples from explicit OMA pairs.

Contract: illustrative evidence, not an orthogroup census or proof of a
conserved mechanism. All displayed genes are checked against the frozen OMA
pair table and the manually reviewed primary-paper evidence ledger.
"""
import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


TASK = Path(__file__).resolve().parents[1]
EVIDENCE = TASK / "notes/focused_experiment_evidence.tsv"
PAIRS = TASK / "results/oma_pairwise/refined/longevity_pairs.csv"
OUT = TASK / "results/focused"
FIGURES = TASK / "figures/focused"
ORDER = ["TOR", "S6 kinase", "Insulin/IGF receptor"]
SPECIES = ["Saccharomyces cerevisiae", "Caenorhabditis elegans",
           "Drosophila melanogaster", "Mus musculus"]
NAMES = ["Budding yeast", "Nematode", "Fruit fly", "Mouse"]
SELECTED = {
    "TOR": {"679", "32", "678", "1809"},
    "S6 kinase": {"578", "583", "1737"},
    "Insulin/IGF receptor": {"163", "368", "365"},
}
DISPLAY = {
    "679": ("TOR1", "Deletion", "Replicative life ↑", "16293764"),
    "32": ("let-363", "RNAi", "Adult life ↑", "14668850"),
    "678": ("Tor", "Dominant-negative", "Adult life ↑", "15186745"),
    "1809": ("Mtor", "Hypomorph", "Median survival ↑", "23994476"),
    "578": ("rsks-1", "RNAi", "Mean life ↑", "17266680"),
    "583": ("S6k", "Dominant-negative", "Adult life ↑", "15186745"),
    "1737": ("Rps6kb1", "Global knockout", "Female life ↑*", "19797661"),
    "163": ("daf-2", "Hypomorph", "Adult life ↑", "8247153"),
    "368": ("InR", "Hypomorph", "Female life ↑", "11292875"),
    "365": ("Igf1r", "Heterozygous KO", "Female life ↑*†", "12483226"),
}


def read_rows(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def write_rows(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def check_and_select():
    evidence = read_rows(EVIDENCE, "\t")
    pairs = read_rows(PAIRS)
    assert len(evidence) == len(DISPLAY) == 10
    assert len({row["genage_id"] for row in evidence}) == len(evidence)
    for row in evidence:
        gid = row["genage_id"]
        assert gid in SELECTED[row["example"]]
        assert row["genage_category"] == "Anti-Longevity"
        assert (row["gene"], row["primary_pmid"]) == (DISPLAY[gid][0], DISPLAY[gid][3])
        assert row["primary_pmid"].isdigit() and row["review_depth"]
        assert row["species"] in SPECIES
    seen = {row["genage_id"] for row in evidence}
    assert seen == set(DISPLAY)
    chosen = []
    for example in ORDER:
        ids = SELECTED[example]
        for pair in pairs:
            if pair["genage_ids1"] in ids and pair["genage_ids2"] in ids:
                assert pair["both_lifespan"] == "True"
                assert pair["influences1"] == pair["influences2"] == "Anti-Longevity"
                chosen.append({"example": example, **pair})
    counts = Counter(row["example"] for row in chosen)
    assert counts == {"TOR": 6, "S6 kinase": 2, "Insulin/IGF receptor": 2}, counts
    for row in evidence:
        assert any(row["genage_id"] in (p["genage_ids1"], p["genage_ids2"])
                   for p in chosen)
    return evidence, chosen, counts


def draw(evidence, counts):
    lookup = {(row["example"], row["species"]): row for row in evidence}
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                         "pdf.fonttype": 42, "font.size": 7.2})
    fig = plt.figure(figsize=(183 / 25.4, 166 / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ink = "#173047"
    muted = "#4E6474"
    line = "#D5E0E7"
    checked_text = []

    def card(x, y, width, height, fill="white", edge=line):
        ax.add_patch(Rectangle((x, y), width, height, facecolor=fill,
                               edgecolor=edge, linewidth=.55))

    def put(x, y, label, bounds, **style):
        artist = ax.text(x, y, label, **style)
        checked_text.append((artist, bounds))
        return artist

    x0, width = .045, .91
    label_w = .185
    cell_w = (width - label_w) / 4

    # Header, evidence grid, edge note and caveats each have explicit bounds.
    title_box = (x0, .875, width, .095)
    card(*title_box, fill="#F7FAFB")
    put(.059, .936, "Selected longevity-associated orthologs in major model organisms",
        title_box, fontsize=11.1, fontweight="bold", color=ink)
    put(.059, .900, "Primary lifespan experiments; explicit OMA relationships; illustrative, not an exhaustive survey",
        title_box, fontsize=7.1, color=muted)

    for c, name in enumerate(NAMES):
        x = x0 + label_w + c * cell_w
        bounds = (x, .812, cell_w - .006, .047)
        card(*bounds, fill="#E8EFF3", edge="#E8EFF3")
        put(x + .009, .835, name, bounds, va="center",
            fontsize=7.4, weight="bold", color=ink)

    row_y = [.635, .454, .273]
    bg = ["#E9F3F4", "#F1F2EC", "#EDF0F8"]
    for i, example in enumerate(ORDER):
        y = row_y[i]
        label_box = (x0, y, label_w - .008, .16)
        card(*label_box, fill=bg[i], edge=bg[i])
        if example == "Insulin/IGF receptor":
            put(x0 + .012, y + .117, "Insulin/IGF", label_box,
                fontsize=8.2, weight="bold", color=ink)
            put(x0 + .012, y + .094, "receptor", label_box,
                fontsize=8.2, weight="bold", color=ink)
            pair_y = y + .048
        else:
            put(x0 + .012, y + .112, example, label_box,
                fontsize=8.2, weight="bold", color=ink)
            pair_y = y + .065
        put(x0 + .012, pair_y, f"{counts[example]} OMA pairs", label_box,
            fontsize=6.7, color=muted)
        for c, species in enumerate(SPECIES):
            x = x0 + label_w + c * cell_w
            row = lookup.get((example, species))
            cell_box = (x, y, cell_w - .006, .16)
            card(*cell_box, fill="white" if row else "#F5F6F7")
            if row is None:
                put(x + .011, y + .080, "Not selected", cell_box,
                    va="center", fontsize=7, color="#82919A")
                continue
            gene, intervention, outcome, pmid = DISPLAY[row["genage_id"]]
            put(x + .011, y + .123, gene, cell_box,
                fontsize=8.5, weight="bold", color=ink)
            put(x + .011, y + .083, intervention, cell_box,
                fontsize=6.9, color=muted)
            put(x + .011, y + .048, outcome, cell_box,
                fontsize=6.9, color="#147268", weight="bold")
            put(x + .011, y + .014, "PMID " + pmid, cell_box,
                fontsize=5.9, color=muted)

    edge_box = (x0, .166, width, .088)
    card(*edge_box, fill="#F7FAFB")
    put(.059, .222, "Direct OMA edges:", edge_box,
        fontsize=6.8, weight="bold", color=ink)
    put(.215, .222, "TOR: all six pairs (yeast edges 1:many); S6K: worm–fly, fly–mouse;",
        edge_box, fontsize=6.8, color=ink)
    put(.215, .190, "IIS: worm–fly, worm–mouse (1:many). No fly–mouse edge in this selection.",
        edge_box, fontsize=6.8, color=ink)

    caveat_box = (x0, .025, width, .128)
    card(*caveat_box, fill="#FAFBFC")
    put(.059, .123, "* Female-specific mouse benefit; corresponding male increase was not significant.",
        caveat_box, fontsize=6.6, color=muted)
    put(.059, .094, "† The original mouse Igf1r increase was not reproduced in a later cohort (PMID 22132081).",
        caveat_box, fontsize=6.6, color=muted)
    put(.059, .065, "TOR and S6K share a pathway; endpoints differ. Blank cells mean not selected, not no effect.",
        caveat_box, fontsize=6.6, color=muted)
    put(.059, .036, "Igf1r is one of two mouse daf-2 receptor paralogs. Sources: GenAge Build 21, OMA All.May2026, cited papers.",
        caveat_box, fontsize=6.5, color=muted)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for artist, (x, y, box_width, box_height) in checked_text:
        left, bottom = ax.transData.transform((x, y))
        right, top = ax.transData.transform((x + box_width, y + box_height))
        text_box = artist.get_window_extent(renderer=renderer)
        assert (text_box.x0 >= left - 1 and text_box.y0 >= bottom - 1
                and text_box.x1 <= right + 1 and text_box.y1 <= top + 1), (
                    f"Text outside its box: {artist.get_text()}"
                )
    for index, (first, _) in enumerate(checked_text):
        first_box = first.get_window_extent(renderer=renderer)
        for second, _ in checked_text[index + 1:]:
            second_box = second.get_window_extent(renderer=renderer)
            overlap_width = min(first_box.x1, second_box.x1) - max(first_box.x0, second_box.x0)
            overlap_height = min(first_box.y1, second_box.y1) - max(first_box.y0, second_box.y0)
            assert overlap_width <= .5 or overlap_height <= .5, (
                f"Overlapping text: {first.get_text()} / {second.get_text()}"
            )
    for extension in ("png", "svg", "pdf", "tiff"):
        options = {"dpi": 600 if extension == "tiff" else 220}
        if extension == "tiff":
            options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(FIGURES / f"focused_ortholog_examples.{extension}", **options)
    plt.close(fig)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    evidence, chosen, counts = check_and_select()
    write_rows(OUT / "selected_oma_pairs.csv", chosen, list(chosen[0]))
    draw(evidence, counts)
    (OUT / "verification.json").write_text(json.dumps({
        "selected_genage_records": len(evidence),
        "selected_explicit_pairs": len(chosen),
        "pairs_by_example": dict(counts),
        "all_gene_labels_and_pmids_reconciled": True,
        "figure_is_selected_examples_not_orthogroups": True,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Plotted {len(evidence)} records and {len(chosen)} explicit pairs")


if __name__ == "__main__":
    main()
