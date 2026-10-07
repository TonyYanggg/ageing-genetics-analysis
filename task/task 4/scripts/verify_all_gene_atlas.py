"""Prove full gene/pair coverage and inspect every atlas page and export."""
import csv
import hashlib
import html
import json
from collections import Counter, defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

import fitz
from PIL import Image

TASK = Path(__file__).resolve().parents[1]
PROJECT = TASK.parents[1]
OUT = TASK / "results/all_gene_atlas"
FIG = TASK / "figures/all_gene_atlas"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    pairs = read(TASK / "results/complete_four_species/all_qualifying_pairs.csv")
    cards = read(OUT / "gene_cards.csv")
    rendered_pairs = read(OUT / "all_displayed_pairs.csv")
    index = read(OUT / "page_index.csv")
    summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
    raw = {r["GenAge ID"]: r for r in read(PROJECT / "GenAge Project/models_genes/genage_models.csv")}
    reviewed = {r["genage_id"]: r for r in read(TASK / "results/expanded/displayed_evidence.csv")}
    refs = defaultdict(set)
    for row in read(PROJECT / "task/task 3/results/verified_timeline/entry_references.csv"):
        refs[row["genage_id"]].add(row["pubmed_url"].rstrip("/").split("/")[-1])
    key = lambda p: (p["species1"], p["entrez1"], p["species2"], p["entrez2"])
    expected_pairs = {key(p): p for p in pairs}
    assert len(rendered_pairs) == len(expected_pairs) == 206
    assert len({key(p) for p in rendered_pairs}) == 206
    for p in rendered_pairs:
        assert {k: p[k] for k in pairs[0]} == expected_pairs[key(p)]
    expected_nodes = {(p["species" + s], p["entrez" + s]) for p in pairs for s in ("1", "2")}
    actual = {(c["species"], c["entrez"]): c for c in cards}
    assert set(actual) == expected_nodes and len(cards) == len(actual) == 320
    graph = defaultdict(set)
    for p in pairs:
        a, b = (p["species1"], p["entrez1"]), (p["species2"], p["entrez2"])
        graph[a].add(b); graph[b].add(a)
    unseen, expected_groups = set(graph), set()
    while unseen:
        seen, pending = set(), {next(iter(unseen))}
        while pending:
            node = pending.pop()
            if node not in seen:
                seen.add(node)
                pending.update(graph[node] - seen)
        unseen -= seen
        expected_groups.add(frozenset(seen))
    by_group = defaultdict(set)
    by_page = defaultdict(list)
    all_genage = set()
    for card in cards:
        node = (card["species"], card["entrez"])
        by_group[card["candidate_id"]].add(node)
        by_page[int(card["page"])].append(card)
        gids = card["genage_ids"].split(";")
        all_genage.update(gids)
        assert all(raw[g]["organism"] == card["species"] for g in gids)
        assert all(raw[g]["entrez gene id"] == card["entrez"] for g in gids)
        assert set(card["symbols"].split(";")) == {raw[g]["symbol"] for g in gids}
        # Compare original catalogue fields against every source-pair occurrence.
        for p in pairs:
            for side in ("1", "2"):
                if (p["species" + side], p["entrez" + side]) == node:
                    for field in ("symbols", "genage_ids", "influences", "lifespan_effects"):
                        assert card[field] == p[field + side]
        assert set(card["catalogue_pmids"].split(";")) == set().union(*(refs[g] for g in gids))
        matches = [reviewed[g] for g in gids if g in reviewed]
        normalised_card = " ".join(card["rendered_text"].split())
        assert all(" ".join(symbol.split()) in normalised_card for symbol in card["display_symbol"].split(" / "))
        assert "GenAge " + ", ".join(gids) in normalised_card
        assert card["evidence_type"] in normalised_card
        if matches:
            r = matches[0]
            assert card["evidence_type"] == "Reviewed example"
            for here, there in [("intervention", "display_intervention"), ("outcome", "display_outcome"),
                                ("primary_pmid", "primary_pmid"), ("review_depth", "review_depth"),
                                ("qualification", "qualification"), ("display_symbol", "gene")]:
                assert card[here] == r[there]
            for value in (card["intervention"], card["outcome"], "PMID " + card["primary_pmid"]):
                assert value in card["rendered_text"].replace("\n", " ")
        else:
            assert card["evidence_type"] == "GenAge annotation"
            assert not card["intervention"] and not card["primary_pmid"] and not card["outcome"]
            assert "↑" not in card["rendered_text"] and "↓" not in card["rendered_text"]
            assert card["influences"] in card["rendered_text"].replace("\n", " ")
            assert card["lifespan_effects"].replace("Increase and Decrease", "Increase + Decrease") in card["rendered_text"].replace("\n", " ")
            assert all(pmid in card["rendered_text"] for pmid in card["catalogue_pmids"].split(";"))
    assert {frozenset(g) for g in by_group.values()} == expected_groups
    assert len(by_group) == len(index) == 133
    assert len(all_genage) == 321
    expected_refs = [r for r in read(PROJECT / "task/task 3/results/verified_timeline/entry_references.csv") if r["genage_id"] in all_genage]
    assert read(OUT / "catalogue_references.csv") == expected_refs
    assert Counter(c["evidence_type"] for c in cards) == {"Reviewed example": 20, "GenAge annotation": 300}
    for row in index:
        members = [c for c in cards if c["candidate_id"] == row["candidate_id"]]
        edges = [p for p in rendered_pairs if p["candidate_id"] == row["candidate_id"]]
        assert len(members) == int(row["gene_count"]) and len(edges) == int(row["pair_count"])
        assert {c["page"] for c in members} == {row["page"]}
        assert all(p["page"] == row["page"] for p in edges)
        for p in edges:
            a, b = actual[(p["species1"], p["entrez1"])], actual[(p["species2"], p["entrez2"])]
            assert a["candidate_id"] == b["candidate_id"] == row["candidate_id"]
            assert p["edge_label"] == f"{a['code']}–{b['code']} ({p['pair_rel_types']})"
    page_count = summary["pages"]
    assert set(by_page) == set(range(1, page_count + 1))
    html_index = (FIG / "index.html").read_text(encoding="utf-8")
    assert html_index.count("<tr>") == 321
    for card in cards:
        assert html.escape(card["symbols"]) in html_index
        assert f"page_{int(card['page']):02d}.png" in html_index
    with fitz.open(FIG / "all_gene_ortholog_atlas.pdf") as pdf:
        assert len(pdf) == page_count
        assert len(pdf.get_toc()) == 133
        for page_no, page in enumerate(pdf, 1):
            assert abs(page.rect.width * 25.4 / 72 - 183) < .05
            assert abs(page.rect.height * 25.4 / 72 - 247) < .05
            for word in page.get_text("words"):
                assert page.rect.contains(fitz.Rect(word[:4])), (page_no, word)
            root = ET.parse(FIG / f"page_{page_no:02d}.svg").getroot()
            labels = ["".join(n.itertext()) for n in root.iter() if n.tag.endswith("}text")]
            pdf_text = " ".join(page.get_text().split())
            with fitz.open(FIG / f"page_{page_no:02d}.pdf") as single:
                assert len(single) == 1 and single[0].get_text() == page.get_text()
            for card in by_page[page_no]:
                for line in card["rendered_text"].splitlines():
                    assert line in labels, (page_no, line)
                    assert " ".join(line.split()) in pdf_text, (page_no, line)
            svg_text = " ".join(labels)
            for p in rendered_pairs:
                if int(p["page"]) == page_no:
                    assert p["edge_label"] in svg_text
            for suffix, dpi in (("png", 220), ("tiff", 600)):
                path = FIG / f"page_{page_no:02d}.{suffix}"
                with Image.open(path) as im:
                    im.verify()
                with Image.open(path) as im:
                    assert abs(im.width - 183 / 25.4 * dpi) <= 1
                    assert abs(im.height - 247 / 25.4 * dpi) <= 1
                    assert abs(im.info["dpi"][0] - dpi) < 1
                    im.load()
                    assert any(lo < hi for lo, hi in im.convert("RGB").getextrema())
            if page_no in {1, 3, page_count} or any(len(g) == 6 and next(iter(g)) in {(c["species"], c["entrez"]) for c in by_page[page_no]} for g in expected_groups):
                page.get_pixmap(matrix=fitz.Matrix(1.8, 1.8)).save(OUT / f"pdf_preview_{page_no:02d}.png")
    bounds = json.loads((OUT / "text_bounds.json").read_text(encoding="utf-8"))
    for item in bounds:
        x, y, w, h = item["box"]
        (left, bottom), (right, top) = item["text_bounds"]
        assert left >= x - .15 and bottom >= y - .15 and right <= x + w + .15 and top <= y + h + .15
    assert len(bounds) == summary["text_elements"]
    for path, digest in summary["input_sha256"].items():
        assert hashlib.sha256((PROJECT / path).read_bytes()).hexdigest() == digest
    report = dict(status="PASS", genes=320, genage_records=321, direct_pairs=206,
                  networks=133, pages=page_count, reviewed_examples=20, catalogue_annotations=300,
                  all_pages_all_card_text_checked=True, all_pairs_rendered=True,
                  all_svg_pdf_png_tiff_exports_checked=True, text_elements=len(bounds))
    (OUT / "verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
