"""Independently verify Task 5 Open Targets outputs against cached responses."""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path


TASK = Path(__file__).resolve().parents[1]
OUT = TASK / "results"
CACHE = OUT / "source_cache" / "opentargets_26_06"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def checked_json(path):
    meta_path = path.with_suffix(".meta.json")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    assert sha256(path) == meta["sha256"], f"Checksum mismatch: {path}"
    return json.loads(path.read_text(encoding="utf-8"))


def source_union(ensembl_id):
    pages = []
    index = 0
    while True:
        path = CACHE / f"{ensembl_id}_page_{index:04d}.json"
        if not path.exists():
            break
        pages.append(checked_json(path))
        index += 1
    assert pages, f"Missing primary cache: {ensembl_id}"
    target = pages[0]["data"]["target"]
    assert target is not None and target["id"] == ensembl_id
    declared = int(target["associatedDiseases"]["count"])
    rows = [r for p in pages for r in p["data"]["target"]["associatedDiseases"]["rows"]]
    assert len(rows) == declared, f"Primary page total differs from declared count: {ensembl_id}"
    by_id = {}
    for row in rows:
        old = by_id.get(row["disease"]["id"])
        assert old is None or old == row, f"Conflicting primary rows: {ensembl_id}"
        by_id[row["disease"]["id"]] = row

    repair_dir = CACHE / "pagination_reconciliation"
    repair_paths = (
        [p for p in repair_dir.glob(f"{ensembl_id}_size_*_page_*.json") if not p.name.endswith(".meta.json")]
        if repair_dir.exists() else []
    )
    for path in sorted(repair_paths):
        data = checked_json(path)
        repair_target = data["data"]["target"]
        assert repair_target["id"] == ensembl_id
        assert int(repair_target["associatedDiseases"]["count"]) == declared
        for row in repair_target["associatedDiseases"]["rows"]:
            old = by_id.get(row["disease"]["id"])
            assert old is None or old == row, f"Conflicting repair rows: {ensembl_id}"
            by_id[row["disease"]["id"]] = row
    assert len(by_id) == declared, f"Unique source union is {len(by_id)}/{declared}: {ensembl_id}"
    return set(by_id), declared


def main():
    audit = read_csv(OUT / "human_gene_audit.csv")
    summaries = read_csv(OUT / "gene_disease_summary.csv")
    assert len(audit) == len(summaries) == 307
    assert len({r["genage_id"] for r in summaries}) == 307
    summary_by_gene = {r["genage_id"]: r for r in summaries}

    expected_by_gene = {}
    source_pairs = set()
    reconciliation_needed = 0
    for gene in audit:
        gene_id = gene["genage_id"]
        ensembl_id = gene["ensembl_gene_ids"]
        summary = summary_by_gene[gene_id]
        if not ensembl_id:
            assert summary["query_status"] == "missing_ensembl_id"
            assert int(summary["direct_association_count"]) == 0
            continue
        disease_ids, declared = source_union(ensembl_id)
        expected_by_gene[gene_id] = disease_ids
        source_pairs.update((gene_id, disease_id) for disease_id in disease_ids)
        assert int(summary["direct_association_count"]) == declared
        primary_paths = [
            p for p in CACHE.glob(f"{ensembl_id}_page_*.json") if not p.name.endswith(".meta.json")
        ]
        if len(primary_paths) > 1:
            primary_ids = []
            for path in sorted(primary_paths):
                data = checked_json(path)
                primary_ids.extend(r["disease"]["id"] for r in data["data"]["target"]["associatedDiseases"]["rows"])
            reconciliation_needed += len(set(primary_ids)) < declared

    output_pairs = set()
    class_by_gene = defaultdict(Counter)
    with (OUT / "opentargets_gene_disease.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            pair = (row["genage_id"], row["disease_id"])
            assert pair not in output_pairs, f"Duplicate output pair: {pair}"
            output_pairs.add(pair)
            assert row["direct_association"] == "true"
            assert row["opentargets_data_version"] == "26.06"
            class_by_gene[row["genage_id"]][row["term_class"]] += 1
    assert output_pairs == source_pairs, "Output gene-disease pairs differ from cached source union"

    class_columns = {
        "disease_or_clinical_condition": "disease_or_clinical_condition_count",
        "phenotype_or_trait": "phenotype_or_trait_count",
        "measurement": "measurement_count",
        "biological_process": "biological_process_count",
        "medical_procedure": "medical_procedure_count",
        "nonhuman_or_other": "nonhuman_or_other_count",
    }
    for gene_id, summary in summary_by_gene.items():
        for class_name, column in class_columns.items():
            assert class_by_gene[gene_id][class_name] == int(summary[column])
        assert sum(class_by_gene[gene_id].values()) == int(summary["direct_association_count"])

    datatype_rows = 0
    with (OUT / "opentargets_datatype_scores.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            assert (row["genage_id"], row["disease_id"]) in output_pairs
            datatype_rows += 1
    datasource_rows = 0
    with (OUT / "opentargets_datasource_scores.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            assert (row["genage_id"], row["disease_id"]) in output_pairs
            datasource_rows += 1

    collection = json.loads((OUT / "opentargets_collection_summary.json").read_text(encoding="utf-8"))
    for name, expected_hash in collection["output_sha256"].items():
        assert sha256(OUT / name) == expected_hash, f"Output checksum mismatch: {name}"
    assert collection["unique_gene_disease_pairs"] == len(output_pairs)
    assert collection["datatype_score_rows"] == datatype_rows
    assert collection["datasource_score_rows"] == datasource_rows

    report = {
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "checks": {
            "all_307_genage_genes_present": True,
            "queryable_genes_reconstructed_from_checksummed_cache": len(expected_by_gene),
            "explicit_unqueryable_genes": 307 - len(expected_by_gene),
            "unique_output_pairs_equal_source_union": len(output_pairs),
            "per_gene_summary_counts_match": True,
            "evidence_score_rows_reference_valid_pairs": True,
            "output_hashes_match_collection_manifest": True,
            "genes_needing_pagination_reconciliation": reconciliation_needed,
        },
        "datatype_score_rows": datatype_rows,
        "datasource_score_rows": datasource_rows,
    }
    (OUT / "verification_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
