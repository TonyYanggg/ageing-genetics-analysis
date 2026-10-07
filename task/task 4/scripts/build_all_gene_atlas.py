"""Render every qualifying Task 4 gene in the approved four-column card style.

Contract: all 320 species/Entrez genes and 206 explicit pairs, with no selection
threshold. Connected networks are layout units, not inferred orthogroups. Use
reviewed interventions only where a ledger exists; otherwise label GenAge data.
Python evidence grid, 183 x 247 mm pages; PDF/SVG/PNG and 600-dpi TIFF.
"""
import csv
import hashlib
import html
import json
import textwrap
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from PIL import Image, ImageDraw

TASK = Path(__file__).resolve().parents[1]
PROJECT = TASK.parents[1]
OUT = TASK / "results/all_gene_atlas"
FIG = TASK / "figures/all_gene_atlas"
SPECIES = ["Saccharomyces cerevisiae", "Caenorhabditis elegans", "Drosophila melanogaster", "Mus musculus"]
NAMES = ["Budding yeast", "Nematode", "Fruit fly", "Mouse"]
CODES = dict(zip(SPECIES, "YWFM"))
ORDER = ["TOR", "S6 kinase", "Insulin/IGF receptor", "AMPK", "ATG7", "SOD2"]
NOTES = {
    "TOR": "Yeast replicative lifespan differs from animal survival. Yeast TOR1 pairs are non-1:1; TOR and S6K share a pathway.",
    "S6 kinase": "* Mouse benefit is female-specific; male increase was not significant. The 2011 correction concerns a worm strain, not mouse survival. No worm–mouse pair is implied.",
    "Insulin/IGF receptor": "* Original female Igf1r benefit; male increase not significant. † Not reproduced in a later cohort (PMID 22132081). Insr is a separate paralog; no fly–mouse pair is implied.",
    "AMPK": "Yeast result is replicative lifespan; worm/fly interventions differ. No direct yeast–worm pair is implied.",
    "ATG7": "‡ Worm cell: dietary-restricted eat-2 mutants, not a general wild-type result. GenAge labels differ from this selected experiment; timing and background matter.",
    "SOD2": "§ Yeast: loss of sch9Δ chronological lifespan extension. Mouse: connective-tissue deficiency. Worm deletion effects depend on genetic background.",
}


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def wrapped(value, width):
    return textwrap.wrap(value, width=width, break_long_words=False, break_on_hyphens=False) or [""]


def prepare():
    pairs = read(TASK / "results/complete_four_species/all_qualifying_pairs.csv")
    reviewed = {r["genage_id"]: r for r in read(TASK / "results/expanded/displayed_evidence.csv")}
    refs = defaultdict(set)
    reference_rows = read(PROJECT / "task/task 3/results/verified_timeline/entry_references.csv")
    for row in reference_rows:
        refs[row["genage_id"]].add(row["pubmed_url"].rstrip("/").split("/")[-1])
    nodes, graph = {}, defaultdict(set)
    for pair in pairs:
        keys = []
        for side in ("1", "2"):
            key = (pair["species" + side], pair["entrez" + side])
            data = {field: pair[field + side] for field in ("species", "entrez", "symbols", "genage_ids", "influences", "lifespan_effects")}
            if key in nodes:
                assert nodes[key] == data
            nodes[key] = data
            keys.append(key)
        graph[keys[0]].add(keys[1])
        graph[keys[1]].add(keys[0])
    screen = read(TASK / "results/expanded/candidate_screen.csv")
    by_ids = {frozenset(r["gene_ids"].split("|")): r for r in screen}
    unseen, groups = set(nodes), []
    while unseen:
        start = min(unseen)
        component, pending = set(), [start]
        while pending:
            node = pending.pop()
            if node in component:
                continue
            component.add(node)
            pending.extend(graph[node] - component)
        unseen -= component
        source = by_ids[frozenset(nodes[n]["genage_ids"] for n in component)]
        example = source["selected_example"]
        members = sorted(component, key=lambda n: (SPECIES.index(n[0]), nodes[n]["symbols"], int(n[1])))
        title = example or nodes[members[0]]["symbols"].split(";")[0]
        cards, code_counts = [], Counter()
        for key in members:
            n = nodes[key]
            code_counts[n["species"]] += 1
            code = CODES[n["species"]] + str(code_counts[n["species"]])
            entries = n["genage_ids"].split(";")
            pmids = sorted(set().union(*(refs[gid] for gid in entries)), key=int)
            assert pmids and all(p.isdigit() for p in pmids)
            matches = [reviewed[g] for g in entries if g in reviewed]
            assert len(matches) <= 1
            r = matches[0] if matches else None
            display_symbol = r["gene"] if r else n["symbols"].replace(";", " / ")
            card = dict(candidate_id=source["candidate_id"], candidate_label=title,
                        node_key="|".join(key), code=code, **n,
                        evidence_type="Reviewed example" if r else "GenAge annotation",
                        display_symbol=display_symbol,
                        intervention=r["display_intervention"] if r else "",
                        outcome=r["display_outcome"] if r else "",
                        primary_pmid=r["primary_pmid"] if r else "",
                        review_depth=r["review_depth"] if r else "Catalogue only; intervention not reviewed",
                        qualification=r["qualification"] if r else "Catalogue labels do not specify an intervention or comparable causal effect.",
                        catalogue_pmids=";".join(pmids))
            cards.append(card)
        code_lookup = {tuple(c["node_key"].split("|")): c["code"] for c in cards}
        edges = []
        for p in pairs:
            a, b = (p["species1"], p["entrez1"]), (p["species2"], p["entrez2"])
            if a in component:
                assert b in component
                edges.append(dict(candidate_id=source["candidate_id"],
                                  edge_label=f"{code_lookup[a]}–{code_lookup[b]} ({p['pair_rel_types']})", **p))
        groups.append(dict(id=source["candidate_id"], title=title, example=example,
                           species_count=len({n[0] for n in component}), cards=cards, edges=edges))
    # Approved examples first, then all remaining >=3-species and two-species sets.
    groups.sort(key=lambda g: (0, ORDER.index(g["example"])) if g["example"] else
                (1 if g["species_count"] >= 3 else 2, -g["species_count"], g["title"].casefold(), g["id"]))
    used_ids = {gid for n in nodes.values() for gid in n["genage_ids"].split(";")}
    write(OUT / "catalogue_references.csv", [r for r in reference_rows if r["genage_id"] in used_ids])
    return groups


def card_lines(card):
    # Text and vertical spacing form a data record as well as the rendered card.
    lines = [(card["code"] + "  " + line, 7.8, "bold", "ink") if i == 0 else
             (line, 7.8, "bold", "ink") for i, line in enumerate(card["display_symbol"].split(" / "))]
    lines += [("GenAge " + card["genage_ids"].replace(";", ", "), 5.8, "normal", "muted"),
              (card["evidence_type"], 6.0, "bold", "muted")]
    if card["evidence_type"] == "Reviewed example":
        lines += [(s, 6.4, "normal", "muted") for s in wrapped(card["intervention"], 25)]
        lines += [(s, 6.5, "bold", "up" if "↑" in s else "down") for s in wrapped(card["outcome"], 25)]
        lines += [("PMID " + card["primary_pmid"], 5.8, "normal", "muted")]
    else:
        lines += [(s, 6.4, "bold", "ink") for s in wrapped(card["influences"], 25)]
        effect = card["lifespan_effects"].replace("Increase and Decrease", "Increase + Decrease")
        lines += [(s, 6.2, "normal", "muted") for s in wrapped("Lifespan: " + effect, 25)]
        values = card["catalogue_pmids"].split(";")
        lines += [("Refs: " + ", ".join(values[:2]), 5.6, "normal", "muted")]
        lines += [("      " + ", ".join(values[i:i + 2]), 5.6, "normal", "muted") for i in range(2, len(values), 2)]
    return lines


def measure(groups):
    for group in groups:
        for card in group["cards"]:
            card["lines"] = card_lines(card)
            card["height"] = len(card["lines"]) * 3.65 + 4
        group["grid_height"] = max(sum(c["height"] + 1.5 for c in group["cards"] if c["species"] == s) - 1.5 for s in SPECIES)
        group["edge_lines"] = wrapped("Direct pairs: " + "; ".join(e["edge_label"] for e in group["edges"]), 114)
        group["note_lines"] = wrapped(NOTES[group["example"]], 115) if group["example"] else []
        group["note_height"] = 3 + 3.4 * (len(group["edge_lines"]) + len(group["note_lines"]))
        group["height"] = group["grid_height"] + 1 + group["note_height"] + 3
        assert group["height"] <= 186, (group["id"], group["height"])
    pages, page, remaining = [], [], 186
    for group in groups:
        if page and group["height"] > remaining:
            pages.append(page)
            page, remaining = [], 186
        page.append(group)
        remaining -= group["height"]
    if page:
        pages.append(page)
    return pages


def render_page(groups, page_no, page_count, combined, reuse_existing=False):
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(183 / 25.4, 247 / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 183), ylim=(0, 247))
    ax.axis("off")
    colours = dict(ink="#173047", muted="#4E6474", up="#147268", down="#995B28")
    checked = []

    def box(bounds, fill="white"):
        ax.add_patch(Rectangle(bounds[:2], bounds[2], bounds[3], facecolor=fill, edgecolor="#D5E0E7", linewidth=.55))

    def put(value, x, top, bounds, size=6.3, weight="normal", colour="ink"):
        artist = ax.text(x, top, value, va="top", fontsize=size, weight=weight, color=colours[colour])
        checked.append((artist, bounds))

    header = (5, 221, 173, 22)
    box(header, "#F7FAFB")
    put("Longevity orthologs across four major model organisms", 8, 240, header, 10.2, "bold")
    put(f"Complete visual catalogue • 320 genes • 206 direct pairs                       Page {page_no} / {page_count}", 8, 234, header, 7)
    put("Reviewed example = a checked experiment; GenAge annotation = catalogue evidence only.", 8, 229, header, 6.5, colour="muted")
    put("All genes in the qualifying pair set are included; rows organise direct links, not proven conserved mechanisms.", 8, 225, header, 6.15, colour="muted")
    label_w, cell_w = 30, 35.75
    for col, name in enumerate(NAMES):
        x = 5 + label_w + col * cell_w
        bounds = (x, 210, cell_w - 1, 8)
        box(bounds, "#E8EFF3")
        put(name, x + 1.5, 215.5, bounds, 7.0, "bold")
    top = 207
    all_cards, all_edges, contents = [], [], []
    for i, group in enumerate(groups):
        bottom = top - group["grid_height"]
        label = (5, bottom, label_w - 1, group["grid_height"])
        box(label, ["#E9F3F4", "#F1F2EC", "#EDF0F8"][i % 3])
        put(group["id"], 7, top - 2, label, 7.6, "bold")
        title_lines = wrapped(group["title"], 17)
        for j, line in enumerate(title_lines):
            put(line, 7, top - 7 - j * 3.5, label, 7.0, "bold")
        offset = 8 + len(title_lines) * 3.5
        put(f"{len(group['cards'])} genes / {group['species_count']} species", 7, top - offset, label, 5.7, colour="muted")
        pair_word = "pair" if len(group["edges"]) == 1 else "pairs"
        put(f"{len(group['edges'])} direct {pair_word}", 7, top - offset - 3.5, label, 5.9, colour="muted")
        for col, species in enumerate(SPECIES):
            x = 5 + label_w + col * cell_w
            members = [c for c in group["cards"] if c["species"] == species]
            if not members:
                bounds = (x, bottom, cell_w - 1, group["grid_height"])
                box(bounds, "#F5F6F7")
                put("No gene in this", x + 2, top - 9, bounds, 6.5, colour="muted")
                put("candidate network", x + 2, top - 13, bounds, 6.5, colour="muted")
            cursor = top
            for card in members:
                bounds = (x, cursor - card["height"], cell_w - 1, card["height"])
                box(bounds)
                for j, (value, size, weight, colour) in enumerate(card["lines"]):
                    put(value, x + 1.5, cursor - 2 - j * 3.65, bounds, size, weight, colour)
                payload = {k: v for k, v in card.items() if k not in {"lines", "height"}}
                payload.update(page=page_no, pdf_file="all_gene_ortholog_atlas.pdf", page_file=f"page_{page_no:02d}.png",
                               rendered_text="\n".join(t[0] for t in card["lines"]))
                all_cards.append(payload)
                cursor -= card["height"] + 1.5
        note = (5, bottom - 1 - group["note_height"], 173, group["note_height"])
        box(note, "#F7FAFB")
        for j, line in enumerate(group["edge_lines"] + group["note_lines"]):
            put(line, 7, bottom - 2.3 - j * 3.4, note, 6.1, colour="muted")
        all_edges.extend({"page": page_no, **e} for e in group["edges"])
        contents.append(dict(page=page_no, candidate_id=group["id"], label=group["title"],
                             species_count=group["species_count"], gene_count=len(group["cards"]),
                             pair_count=len(group["edges"]), symbols="; ".join(c["symbols"] for c in group["cards"])))
        top -= group["height"]
    footer = (5, 4, 173, 16)
    box(footer, "#FAFBFC")
    for j, line in enumerate([
        "Y/W/F/M identify genes in yeast/worm/fly/mouse columns. Pair codes retain the OMA source orientation.",
        "Refs are catalogue-linked PMIDs, not newly reviewed experiments. Empty cells do not establish biological absence.",
        "GenAge lifespan labels do not specify interventions. Only reviewed examples have experiment arrows and colours.",
        "Sources: GenAge Build 21; OMA All.May2026. See atlas legend, gene index and source tables for evidence details.",
    ]):
        put(line, 7, 18 - 3.3 * j, footer, 6.05, colour="muted")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    measurements = []
    for artist, bounds in checked:
        bb = artist.get_window_extent(renderer)
        x, y, w, h = bounds
        lo, hi = ax.transData.transform((x, y)), ax.transData.transform((x + w, y + h))
        assert bb.x0 >= lo[0] - .5 and bb.y0 >= lo[1] - .5 and bb.x1 <= hi[0] + .5 and bb.y1 <= hi[1] + .5, (page_no, artist.get_text())
        measurements.append(dict(page=page_no, text=artist.get_text(), box=list(bounds),
                                 text_bounds=ax.transData.inverted().transform([[bb.x0, bb.y0], [bb.x1, bb.y1]]).tolist()))
    measured = [(artist, artist.get_window_extent(renderer)) for artist, _ in checked]
    for i, (a, aa) in enumerate(measured):
        for b, bb in measured[i + 1:]:
            assert min(aa.x1, bb.x1) - max(aa.x0, bb.x0) <= .5 or min(aa.y1, bb.y1) - max(aa.y0, bb.y0) <= .5, (page_no, a.get_text(), b.get_text())
    reuse = reuse_existing and all((FIG / f"page_{page_no:02d}.{s}").exists() for s in ("png", "svg", "pdf", "tiff"))
    if reuse:
        root = ET.parse(FIG / f"page_{page_no:02d}.svg").getroot()
        old_labels = ["".join(n.itertext()) for n in root.iter() if n.tag.endswith("}text")]
        reuse = old_labels == [artist.get_text() for artist, _ in checked]
    if not reuse:
        for suffix in ("png", "svg", "pdf", "tiff"):
            kwargs = {"dpi": 600 if suffix == "tiff" else 220}
            if suffix == "tiff":
                kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
            fig.savefig(FIG / f"page_{page_no:02d}.{suffix}", **kwargs)
    combined.savefig(fig)
    plt.close(fig)
    print(f"Page {page_no}/{page_count}: {len(all_cards)} genes, {len(all_edges)} direct pairs", flush=True)
    return all_cards, all_edges, contents, measurements


def browse_and_preview(cards, contents, page_count):
    # A local index makes a long figure set searchable by every displayed gene.
    lines = ['<!doctype html><html lang="en"><meta charset="utf-8"><title>Task 4 all-gene atlas</title>',
             '<style>body{font:16px system-ui;max-width:1100px;margin:2em auto;color:#173047}input{padding:.7em;width:70%}table{border-collapse:collapse;width:100%}td,th{padding:.5em;border-bottom:1px solid #d5e0e7;text-align:left}a{color:#147268}iframe{width:100%;height:85vh;border:1px solid #d5e0e7}</style>',
             '<h1>All-gene ortholog atlas</h1><p>320 genes • 206 direct pairs • four major model organisms</p>',
             '<p>Same four-column card layout throughout. Reviewed experiments and catalogue annotations are explicitly labelled.</p>',
             '<p><a href="all_gene_ortholog_atlas.pdf">Open complete PDF</a> · <a href="../../notes/all_gene_atlas_legend.md">Legend</a></p>',
             '<input id="search" placeholder="Find a gene, species, GenAge ID or candidate">',
             '<table id="genes"><thead><tr><th>Gene</th><th>Organism</th><th>GenAge ID</th><th>Candidate</th><th>Page</th></tr></thead><tbody>']
    for card in sorted(cards, key=lambda c: (c["symbols"].casefold(), c["species"])):
        p = card["page"]
        cells = [card["symbols"], card["species"], card["genage_ids"], card["candidate_id"]]
        lines.append('<tr>' + ''.join('<td>' + html.escape(c) + '</td>' for c in cells) +
                     f'<td><a href="page_{p:02d}.png">{p} PNG</a> / <a href="all_gene_ortholog_atlas.pdf#page={p}">PDF</a></td></tr>')
    lines += ['</tbody></table><script>document.getElementById("search").addEventListener("input",e=>{const q=e.target.value.toLowerCase();document.querySelectorAll("#genes tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></html>']
    (FIG / "index.html").write_text("\n".join(lines), encoding="utf-8")
    # Contact sheets expose every rendered page for a visual layout sweep.
    for start in range(1, page_count + 1, 9):
        sheet = Image.new("RGB", (1050, 1500), "#e8eff3")
        pen = ImageDraw.Draw(sheet)
        for n, page in enumerate(range(start, min(start + 9, page_count + 1))):
            with Image.open(FIG / f"page_{page:02d}.png") as im:
                im.thumbnail((335, 460))
                x, y = 10 + (n % 3) * 350, 25 + (n // 3) * 500
                sheet.paste(im.convert("RGB"), (x, y))
                pen.text((x, y - 18), f"Page {page}", fill="#173047")
        sheet.save(OUT / f"contact_{start:02d}.png")


def main(reuse_existing=False):
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    groups = prepare()
    pages = measure(groups)
    cards, edges, contents, measurements = [], [], [], []
    pdf_path = FIG / "all_gene_ortholog_atlas.pdf"
    with PdfPages(pdf_path) as combined:
        for page_no, page in enumerate(pages, 1):
            a, b, c, d = render_page(page, page_no, len(pages), combined, reuse_existing=reuse_existing)
            cards.extend(a); edges.extend(b); contents.extend(c); measurements.extend(d)
    # Add navigable bookmarks without altering the matplotlib page graphics.
    with fitz.open(pdf_path) as pdf:
        pdf.set_toc([[1, f"{r['candidate_id']} | {r['label']}", r["page"]] for r in contents])
        pdf.saveIncr()
    write(OUT / "gene_cards.csv", cards)
    write(OUT / "all_displayed_pairs.csv", edges)
    write(OUT / "page_index.csv", contents)
    (OUT / "text_bounds.json").write_text(json.dumps(measurements, indent=2), encoding="utf-8")
    browse_and_preview(cards, contents, len(pages))
    inputs = [TASK / "results/complete_four_species/all_qualifying_pairs.csv",
              TASK / "results/expanded/displayed_evidence.csv", TASK / "results/expanded/candidate_screen.csv",
              PROJECT / "task/task 3/results/verified_timeline/entry_references.csv"]
    summary = dict(genes=len(cards), genage_records=len({g for c in cards for g in c["genage_ids"].split(";")}),
                   direct_pairs=len(edges), candidate_networks=len(groups), pages=len(pages),
                   reviewed_examples=sum(c["evidence_type"] == "Reviewed example" for c in cards),
                   catalogue_only=sum(c["evidence_type"] == "GenAge annotation" for c in cards),
                   text_elements=len(measurements), page_mm=[183, 247],
                   input_sha256={str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
