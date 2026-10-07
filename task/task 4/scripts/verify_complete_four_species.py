"""Independently reconcile the complete four-species tables and figure exports."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import fitz
from PIL import Image


TASK = Path(__file__).resolve().parents[1]
ROOT = TASK.parents[1]
SOURCE = TASK / "results/oma_pairwise/refined"
OUT = TASK / "results/complete_four_species"
FIG = TASK / "figures/complete_four_species/complete_ortholog_overview"
FOUR = {"Saccharomyces cerevisiae", "Caenorhabditis elegans",
        "Drosophila melanogaster", "Mus musculus"}
EFFECTS = {"Increase", "Decrease", "Increase and Decrease"}
DIRECTIONS = ("Anti-Longevity", "Pro-Longevity", "discordant_or_mixed", "unresolved")


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def key(row):
    return row["species1"], row["entrez1"], row["species2"], row["entrez2"]


def main():
    raw = read(ROOT / "GenAge Project/models_genes/genage_models.csv")
    original = read(SOURCE / "longevity_pairs.csv")
    complete = read(OUT / "all_qualifying_pairs.csv")
    coverage = read(OUT / "all_nine_coverage.csv")
    matrix = read(OUT / "species_pair_summary.csv")
    multiplicity = read(OUT / "gene_pair_multiplicity.csv")
    summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))

    expected = {key(r): r for r in original if r["species1"] in FOUR and r["species2"] in FOUR}
    actual = {key(r): r for r in complete}
    assert len(expected) == len(actual) == 206 and expected.keys() == actual.keys()
    assert len({(r["species1"], r["species2"]) for r in matrix}) == 6
    gene_records = defaultdict(list)
    for r in raw:
        gene_records[(r["organism"], r["entrez gene id"])].append(r)
    genes = set()
    partners = defaultdict(set)
    by_pair = Counter()
    by_direction = Counter()
    for k, row in actual.items():
        source = expected[k]
        for field in source:
            assert row[field] == source[field], (k, field)
        left, right = (k[0], k[1]), (k[2], k[3])
        assert left < right and left[0] != right[0]
        for side, gene in (("1", left), ("2", right)):
            records = gene_records[gene]
            assert records
            effects = sorted({r["lifespan effect"] for r in records if r["lifespan effect"] in EFFECTS})
            assert effects and row[f"lifespan_effects{side}"] == ";".join(effects)
            assert set(row[f"genage_ids{side}"].split(";")) == {r["GenAge ID"] for r in records}
        genes.update((left, right))
        partners[left].add(right)
        partners[right].add(left)
        species_pair = tuple(sorted((left[0], right[0])))
        by_pair[species_pair] += 1
        by_direction[(species_pair, row["direction_agreement"])] += 1
    assert len(genes) == len(multiplicity) == summary["distinct_genes"] == 320
    assert { (r["species"], r["entrez_id"]) for r in multiplicity } == genes
    for row in multiplicity:
        gene = row["species"], row["entrez_id"]
        assert int(row["partner_genes"]) == len(partners[gene])
        assert int(row["partner_species"]) == len({p[0] for p in partners[gene]})
    assert summary["genes_with_multiple_partners"] == sum(len(p) > 1 for p in partners.values())
    assert summary["genes_with_partners_in_multiple_species"] == sum(len({x[0] for x in p}) > 1 for p in partners.values())
    assert summary["direction_counts"] == dict(Counter(r["direction_agreement"] for r in actual.values()))
    for row in matrix:
        species_pair = tuple(sorted((row["species1"], row["species2"])))
        assert int(row["pairs"]) == by_pair[species_pair]
        for direction in DIRECTIONS:
            assert int(row[direction]) == by_direction[(species_pair, direction)]
    assert sum(int(r["pairs"]) for r in matrix) == len(actual)
    assert len(coverage) == 9
    source_coverage = {r["species"]: r for r in read(SOURCE / "coverage.csv")}
    for row in coverage:
        assert all(row[field] == source_coverage[row["species"]][field] for field in source_coverage[row["species"]])
        assert row["analysis_role"] == ("principal_four" if row["species"] in FOUR else "outside_main_comparison")
    assert summary["genage_records_principal"] == sum(int(r["records"]) for r in coverage if r["species"] in FOUR) == 2138
    assert summary["outside_main_comparison_records"] == 67
    assert set(summary["oma_unavailable_species"]) == {"Podospora anserina", "Mesocricetus auratus"}
    for suffix in ("png", "svg", "pdf", "tiff"):
        path = FIG.with_suffix("." + suffix)
        assert path.is_file() and path.stat().st_size > 1000
    with Image.open(FIG.with_suffix(".tiff")) as image:
        assert image.width >= 4000 and image.height >= 3000
    assert "206 pairs" in FIG.with_suffix(".svg").read_text(encoding="utf-8")
    pdf = fitz.open(FIG.with_suffix(".pdf"))
    assert len(pdf) == 1
    page = pdf[0]
    assert len(page.get_text("words")) > 70
    assert all(0 <= w[0] < w[2] <= page.rect.width and 0 <= w[1] < w[3] <= page.rect.height
               for w in page.get_text("words"))
    print(f"PASS: {len(actual)} exact pairs, {len(genes)} genes, six species comparisons, nine-species coverage, four exports")


if __name__ == "__main__":
    main()
