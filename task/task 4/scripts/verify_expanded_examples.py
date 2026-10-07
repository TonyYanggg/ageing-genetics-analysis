"""Independently reconcile the expanded illustration, screen and export bundle."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

import fitz
from PIL import Image

TASK = Path(__file__).resolve().parents[1]


def read(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def main():
    pairs = read(TASK / "results/complete_four_species/all_qualifying_pairs.csv")
    selected = read(TASK / "results/expanded/selected_oma_pairs.csv")
    evidence = read(TASK / "results/expanded/displayed_evidence.csv")
    expected = {"TOR": 6, "S6 kinase": 2, "Insulin/IGF receptor": 2, "AMPK": 2, "ATG7": 3, "SOD2": 6}
    assert len(pairs) == 206 and len(evidence) == 20
    assert Counter(r["example"] for r in selected) == expected
    key = lambda p: (p["species1"], p["entrez1"], p["species2"], p["entrez2"])
    original = {key(p): p for p in pairs}
    assert len({key(p) for p in selected}) == 21
    for p in selected:
        assert {k: v for k, v in p.items() if k != "example"} == original[key(p)]
    raw = {r["GenAge ID"]: r for r in read(TASK.parents[1] / "GenAge Project/models_genes/genage_models.csv")}
    assert len({r["genage_id"] for r in evidence}) == 20
    for r in evidence:
        gene = raw[r["genage_id"]]
        assert (r["gene"], r["species"], r["genage_category"], r["genage_lifespan_effect"]) == (
            gene["symbol"], gene["organism"], gene["longevity influence"], gene["lifespan effect"])
        assert r["review_depth"] and r["qualification"] and r["primary_pmid"].isdigit()
        assert any(r["genage_id"] in (p["genage_ids1"], p["genage_ids2"]) for p in selected if p["example"] == r["example"])
    # Every original evidence field is preserved; display fields are added only.
    previous = read(TASK / "notes/focused_experiment_evidence.tsv", "\t")
    lookup = {r["genage_id"]: r for r in evidence}
    for old in previous:
        assert all(lookup[old["genage_id"]][k] == v for k, v in old.items())
    for new in read(TASK / "notes/expanded_additional_evidence.tsv", "\t"):
        assert all(lookup[new["genage_id"]][k] == v for k, v in new.items())
        cache = TASK / "results/expanded/papers" / (new["primary_pmid"] + ".json")
        assert json.loads(cache.read_text(encoding="utf-8"))["paper"]["id"] == new["primary_pmid"]
    # Independently enumerate components using taxon/Entrez keys (not builder IDs).
    graph, ids = defaultdict(set), {}
    for p in pairs:
        a, b = (p["species1"], p["entrez1"]), (p["species2"], p["entrez2"])
        graph[a].add(b)
        graph[b].add(a)
        ids[a], ids[b] = p["genage_ids1"], p["genage_ids2"]
    unseen, components = set(graph), []
    while unseen:
        component, pending = set(), {next(iter(unseen))}
        while pending:
            node = pending.pop()
            if node not in component:
                component.add(node)
                pending.update(graph[node] - component)
        unseen -= component
        components.append(component)
    screen = read(TASK / "results/expanded/candidate_screen.csv")
    assert len(components) == len(screen) == 133
    assert sum(len({n[0] for n in c}) >= 3 for c in components) == 22
    assert {frozenset(ids[n] for n in c) for c in components} == {frozenset(r["gene_ids"].split("|")) for r in screen}
    assert sum(int(r["direct_pair_count"]) for r in screen) == 206
    stem = TASK / "figures/expanded/expanded_ortholog_examples"
    svg = ET.parse(stem.with_suffix(".svg"))
    labels = ["".join(n.itertext()) for n in svg.getroot().iter() if n.tag.endswith("}text")]
    for row in evidence:
        for label in (row["gene"], row["display_intervention"], row["display_outcome"], "PMID " + row["primary_pmid"]):
            assert label in labels, label
    with fitz.open(stem.with_suffix(".pdf")) as pdf:
        assert len(pdf) == 1
        page = pdf[0]
        assert abs(page.rect.width * 25.4 / 72 - 183) < .05
        assert abs(page.rect.height * 25.4 / 72 - 247) < .05
        for word in page.get_text("words"):
            assert page.rect.contains(fitz.Rect(word[:4])), word
        for row in evidence:
            assert row["primary_pmid"] in page.get_text()
        # Python-rendered PDF preview for visual inspection, same physical aspect.
        page.get_pixmap(matrix=fitz.Matrix(1.6, 1.6)).save(TASK / "results/expanded/pdf_preview.png")
    for suffix in ("png", "tiff"):
        with Image.open(stem.with_suffix("." + suffix)) as image:
            image.verify()
        with Image.open(stem.with_suffix("." + suffix)) as image:
            dpi = 600 if suffix == "tiff" else 220
            assert abs(image.width - 183 / 25.4 * dpi) <= 1
            assert abs(image.height - 247 / 25.4 * dpi) <= 1
            assert abs(image.info["dpi"][0] - dpi) < 1
    print("PASS: 6 families, 20 genes, 21 exact pairs; 133 candidate networks/22 with >=3 species; source ledger and four exports reconciled.")


if __name__ == "__main__":
    main()
