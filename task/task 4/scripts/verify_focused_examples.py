"""Independent reconciliation of the selected review figure and source tables."""
import csv
from collections import Counter
from pathlib import Path

from PIL import Image


TASK = Path(__file__).resolve().parents[1]


def rows(path, delimiter=","):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def main():
    original = rows(TASK / "results/oma_pairwise/refined/longevity_pairs.csv")
    selected = rows(TASK / "results/focused/selected_oma_pairs.csv")
    evidence = rows(TASK / "notes/focused_experiment_evidence.tsv", "\t")
    assert len(selected) == len(evidence) == 10
    original_keys = {(r["species1"], r["entrez1"], r["species2"], r["entrez2"])
                     for r in original}
    selected_keys = [(r["species1"], r["entrez1"], r["species2"], r["entrez2"])
                     for r in selected]
    assert len(set(selected_keys)) == len(selected)
    assert set(selected_keys) <= original_keys
    expected = {"TOR": 6, "S6 kinase": 2, "Insulin/IGF receptor": 2}
    assert Counter(r["example"] for r in selected) == expected
    for r in evidence:
        assert r["genage_category"] == "Anti-Longevity"
        assert r["genage_id"] and r["primary_pmid"].isdigit()
        assert r["review_depth"] and r["qualification"]
        assert any(r["genage_id"] in (p["genage_ids1"], p["genage_ids2"])
                   for p in selected if p["example"] == r["example"])
    for r in selected:
        assert r["source"] == "OMA All.May2026"
        assert r["unit"] == "explicit_gene_pair"
        assert r["influences1"] == r["influences2"] == "Anti-Longevity"
    figure = TASK / "figures/focused/focused_ortholog_examples"
    for suffix in ("png", "pdf", "svg", "tiff"):
        output = figure.with_suffix("." + suffix)
        assert output.is_file() and output.stat().st_size > 1000
    with Image.open(figure.with_suffix(".png")) as im:
        assert im.width >= 1500 and im.height >= 1200
    with Image.open(figure.with_suffix(".tiff")) as im:
        assert im.width >= 4000 and im.height >= 3000
    svg = figure.with_suffix(".svg").read_text(encoding="utf-8")
    for label in ("TOR1", "let-363", "Rps6kb1", "Igf1r", "PMID 12483226"):
        assert label in svg
    print("PASS: 10 distinct records, 10 exact source pairs, selection counts, labels and four figure exports")


if __name__ == "__main__":
    main()
