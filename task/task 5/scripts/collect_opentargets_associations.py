"""Collect direct Open Targets associations for every queryable GenAge human gene.

The raw GraphQL responses are cached and checksummed before flat result tables are
materialised.  This makes the relatively slow collection resumable and keeps the
source evidence available for an independent verification pass.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import time

import requests


TASK = Path(__file__).resolve().parents[1]
OUT = TASK / "results"
CACHE = OUT / "source_cache" / "opentargets_26_06"
AUDIT = OUT / "human_gene_audit.csv"
LONGEVITY = OUT / "longevitymap_matches.csv"
API = "https://api.platform.opentargets.org/api/v4/graphql"
META_QUERY = """
query ReleaseMetadata {
  meta {
    product
    name
    dataPrefix
    enableDataReleasePrefix
    downloads
    apiVersion { x y z suffix }
    dataVersion { year month iteration }
  }
}
"""
EXPECTED_DATA_VERSION = "26.06"
PAGE_SIZE = 2000
RECONCILIATION_PAGE_SIZES = (3000, 2500, 1500, 1000)

QUERY = """
query DirectAssociations($ensemblId: String!, $pageIndex: Int!, $pageSize: Int!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    associatedDiseases(
      enableIndirect: false
      page: {index: $pageIndex, size: $pageSize}
    ) {
      count
      rows {
        score
        disease {
          id
          name
          ancestors
          therapeuticAreas { id name }
        }
        datasourceScores { id score }
        datatypeScores { id score }
      }
    }
  }
}
"""


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def request_json(session, method, url, **kwargs):
    last_error = None
    for attempt in range(6):
        try:
            response = session.request(method, url, timeout=120, **kwargs)
            response.raise_for_status()
            data = response.json()
            if isinstance(data, dict) and data.get("errors"):
                raise RuntimeError(json.dumps(data["errors"], ensure_ascii=False))
            return data, dict(response.headers)
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            last_error = exc
            if attempt == 5:
                break
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Request failed after retries: {last_error}")


def freeze_release(session):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / "source_release.json"
    if path.exists():
        record = json.loads(path.read_text(encoding="utf-8"))
    else:
        data, headers = request_json(session, "POST", API, json={"query": META_QUERY})
        raw = data["data"]["meta"]
        dv = raw["dataVersion"]
        av = raw["apiVersion"]
        data_version = f"{dv['year']}.{dv['month']}" + (f".{dv['iteration']}" if dv.get("iteration") else "")
        api_version = f"{av['x']}.{av['y']}.{av['z']}" + (f"-{av['suffix']}" if av.get("suffix") else "")
        normalised = dict(raw)
        normalised["dataVersionComponents"] = normalised.pop("dataVersion")
        normalised["apiVersionComponents"] = normalised.pop("apiVersion")
        normalised["dataVersion"] = data_version
        normalised["apiVersion"] = api_version
        record = {
            "endpoint": API,
            "query": META_QUERY,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "response_headers": headers,
            "response": normalised,
        }
        path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    meta = record["response"]
    if meta.get("dataVersion") != EXPECTED_DATA_VERSION:
        raise RuntimeError(
            f"Expected Open Targets {EXPECTED_DATA_VERSION}, received {meta.get('dataVersion')!r}. "
            "Review the source before mixing releases."
        )
    return record


def cache_page(session, ensembl_id, page_index):
    CACHE.mkdir(parents=True, exist_ok=True)
    stem = f"{ensembl_id}_page_{page_index:04d}"
    path = CACHE / f"{stem}.json"
    meta_path = CACHE / f"{stem}.meta.json"
    variables = {
        "ensemblId": ensembl_id,
        "pageIndex": page_index,
        "pageSize": PAGE_SIZE,
    }
    if path.exists():
        record = json.loads(meta_path.read_text(encoding="utf-8"))
        if record["endpoint"] != API or record["variables"] != variables:
            raise RuntimeError(f"Cached request metadata mismatch: {path}")
        if sha256(path) != record["sha256"]:
            raise RuntimeError(f"Cached response checksum mismatch: {path}")
        return json.loads(path.read_text(encoding="utf-8"))

    data, headers = request_json(
        session,
        "POST",
        API,
        json={"query": QUERY, "variables": variables},
    )
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    path.write_bytes(payload)
    meta = {
        "endpoint": API,
        "variables": variables,
        "direct_only": True,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "response_headers": headers,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
    }
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return data


def cache_reconciliation_page(session, ensembl_id, page_index, page_size):
    folder = CACHE / "pagination_reconciliation"
    folder.mkdir(parents=True, exist_ok=True)
    stem = f"{ensembl_id}_size_{page_size:04d}_page_{page_index:04d}"
    path = folder / f"{stem}.json"
    meta_path = folder / f"{stem}.meta.json"
    variables = {"ensemblId": ensembl_id, "pageIndex": page_index, "pageSize": page_size}
    if path.exists():
        record = json.loads(meta_path.read_text(encoding="utf-8"))
        if record["endpoint"] != API or record["variables"] != variables:
            raise RuntimeError(f"Reconciliation metadata mismatch: {path}")
        if sha256(path) != record["sha256"]:
            raise RuntimeError(f"Reconciliation response checksum mismatch: {path}")
        return json.loads(path.read_text(encoding="utf-8"))
    data, headers = request_json(
        session, "POST", API, json={"query": QUERY, "variables": variables}
    )
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    path.write_bytes(payload)
    meta = {
        "endpoint": API,
        "variables": variables,
        "direct_only": True,
        "purpose": "recover rows obscured by unstable ordering among tied association scores",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "response_headers": headers,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
    }
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return data


def collect_target(session, ensembl_id):
    pages = []
    page_index = 0
    expected_count = None
    while True:
        data = cache_page(session, ensembl_id, page_index)
        target = data.get("data", {}).get("target")
        if target is None:
            return None, pages
        block = target["associatedDiseases"]
        if expected_count is None:
            expected_count = int(block["count"])
        elif expected_count != int(block["count"]):
            raise RuntimeError(f"Association count changed during paging for {ensembl_id}")
        pages.append(data)
        obtained = sum(len(p["data"]["target"]["associatedDiseases"]["rows"]) for p in pages)
        if obtained >= expected_count:
            if obtained != expected_count:
                raise RuntimeError(f"Retrieved {obtained} rows for declared count {expected_count}: {ensembl_id}")
            return target, pages
        if not block["rows"]:
            raise RuntimeError(f"Empty page before count was reached: {ensembl_id}")
        page_index += 1


def classify_term(disease):
    disease_id = disease.get("id", "")
    areas = disease.get("therapeuticAreas") or []
    area_names = {a.get("name", "").lower() for a in areas}
    if disease_id.startswith(("MONDO_", "Orphanet_")):
        return "disease_or_clinical_condition"
    if disease_id.startswith("HP_"):
        return "phenotype_or_trait"
    if disease_id.startswith(("GO_", "OBA_", "GSSO_", "UBERON_")):
        return "biological_process"
    if disease_id.startswith("MP_"):
        return "nonhuman_or_other"
    if any(a.get("id") == "EFO_0001444" or a.get("name", "").lower() == "measurement" for a in areas):
        return "measurement"
    nonclinical = {"phenotype", "biological_process", "medical procedure", "non-human animal disease"}
    if area_names - nonclinical:
        return "disease_or_clinical_condition"
    if "phenotype" in area_names:
        return "phenotype_or_trait"
    if "biological_process" in area_names:
        return "biological_process"
    if "medical procedure" in area_names:
        return "medical_procedure"
    return "nonhuman_or_other"


def load_cached_pages(ensembl_id):
    pages = []
    index = 0
    while True:
        path = CACHE / f"{ensembl_id}_page_{index:04d}.json"
        if not path.exists():
            break
        meta_path = CACHE / f"{ensembl_id}_page_{index:04d}.meta.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if sha256(path) != meta["sha256"]:
            raise RuntimeError(f"Checksum failure while materialising {path}")
        pages.append(json.loads(path.read_text(encoding="utf-8")))
        index += 1
    return pages


def reconcile_rows(session, ensembl_id, pages):
    """Recover a complete unique set when offset pages overlap on tied scores."""
    target = pages[0]["data"]["target"]
    declared = int(target["associatedDiseases"]["count"])
    by_disease = {}
    raw_rows = [
        row for page in pages
        for row in page["data"]["target"]["associatedDiseases"]["rows"]
    ]

    def add_rows(rows):
        for row in rows:
            disease_id = row["disease"]["id"]
            previous = by_disease.get(disease_id)
            if previous is not None and previous != row:
                raise RuntimeError(f"Conflicting rows for {ensembl_id} and {disease_id}")
            by_disease[disease_id] = row

    add_rows(raw_rows)
    initial_unique = len(by_disease)
    attempts = []
    if initial_unique == declared:
        return list(by_disease.values()), {
            "required": False, "declared": declared, "initial_unique": initial_unique,
            "final_unique": initial_unique, "attempts": attempts,
        }
    for page_size in RECONCILIATION_PAGE_SIZES:
        page_index = 0
        variant_rows = []
        while len(variant_rows) < declared:
            data = cache_reconciliation_page(session, ensembl_id, page_index, page_size)
            variant_target = data.get("data", {}).get("target")
            if variant_target is None:
                raise RuntimeError(f"Target disappeared during reconciliation: {ensembl_id}")
            block = variant_target["associatedDiseases"]
            if int(block["count"]) != declared:
                raise RuntimeError(f"Association count changed during reconciliation: {ensembl_id}")
            variant_rows.extend(block["rows"])
            if not block["rows"] and len(variant_rows) < declared:
                raise RuntimeError(f"Empty reconciliation page: {ensembl_id}")
            page_index += 1
        if len(variant_rows) != declared:
            raise RuntimeError(f"Reconciliation returned {len(variant_rows)} of {declared}: {ensembl_id}")
        add_rows(variant_rows)
        attempts.append({"page_size": page_size, "pages": page_index, "unique_after_union": len(by_disease)})
        if len(by_disease) == declared:
            return list(by_disease.values()), {
                "required": True, "declared": declared, "initial_unique": initial_unique,
                "final_unique": len(by_disease), "attempts": attempts,
            }
        if len(by_disease) > declared:
            raise RuntimeError(f"Union exceeded declared association count: {ensembl_id}")
    raise RuntimeError(
        f"Could not recover all unique associations for {ensembl_id}: {len(by_disease)}/{declared}"
    )


def main():
    audit = read_csv(AUDIT)
    if len(audit) != 307 or len({r["genage_id"] for r in audit}) != 307:
        raise RuntimeError("The audited GenAge universe is not exactly 307 unique genes")
    longevity = read_csv(LONGEVITY)
    longevity_by_gene = defaultdict(list)
    for row in longevity:
        longevity_by_gene[row["genage_id"]].append(row)

    session = requests.Session()
    session.headers.update({"User-Agent": "GenAge-Task5-evidence-audit/1.0"})
    release = freeze_release(session)

    queryable = [r for r in audit if r["ensembl_gene_ids"]]
    for position, gene in enumerate(queryable, 1):
        ensembl_ids = [x for x in gene["ensembl_gene_ids"].split(";") if x]
        if len(ensembl_ids) != 1:
            raise RuntimeError(f"Manual resolution required for multiple Ensembl IDs: {gene['genage_symbol']}")
        target, pages = collect_target(session, ensembl_ids[0])
        state = "target_not_found" if target is None else f"{target['associatedDiseases']['count']} direct associations"
        print(f"Open Targets {position:03d}/{len(queryable)} {gene['genage_symbol']} {state}", flush=True)

    association_rows = []
    datatype_rows = []
    datasource_rows = []
    summary_rows = []
    seen_pairs = set()
    retrieval_times = []
    reconciliation_audit = []

    for gene in audit:
        ensembl_id = gene["ensembl_gene_ids"]
        route_status = "missing_ensembl_id"
        target_id = ""
        target_symbol = ""
        rows = []
        if ensembl_id:
            pages = load_cached_pages(ensembl_id)
            if not pages:
                raise RuntimeError(f"No cached Open Targets response: {ensembl_id}")
            target = pages[0].get("data", {}).get("target")
            if target is None:
                route_status = "target_not_found"
            else:
                route_status = "queried_direct_associations"
                target_id = target["id"]
                target_symbol = target["approvedSymbol"]
                if target_id != ensembl_id:
                    raise RuntimeError(f"Target ID mismatch: requested {ensembl_id}, received {target_id}")
                declared = int(target["associatedDiseases"]["count"])
                rows, reconciliation = reconcile_rows(session, ensembl_id, pages)
                reconciliation["ensembl_gene_id"] = ensembl_id
                reconciliation["genage_symbol"] = gene["genage_symbol"]
                reconciliation_audit.append(reconciliation)
                if len(rows) != declared:
                    raise RuntimeError(f"Reconciled {len(rows)} of {declared}: {ensembl_id}")
                for page_index in range(len(pages)):
                    meta = json.loads((CACHE / f"{ensembl_id}_page_{page_index:04d}.meta.json").read_text())
                    retrieval_times.append(meta["retrieved_utc"])

        class_counts = Counter()
        scores = []
        for row in rows:
            disease = row["disease"]
            pair = (gene["genage_id"], disease["id"])
            if pair in seen_pairs:
                raise RuntimeError(f"Duplicate gene-disease pair: {pair}")
            seen_pairs.add(pair)
            term_class = classify_term(disease)
            class_counts[term_class] += 1
            score = float(row["score"])
            scores.append(score)
            base = {
                "genage_id": gene["genage_id"],
                "genage_symbol": gene["genage_symbol"],
                "current_symbol": gene["current_symbol"],
                "entrez_id": gene["entrez_id"],
                "ensembl_gene_id": ensembl_id,
                "opentargets_target_id": target_id,
                "opentargets_target_symbol": target_symbol,
                "disease_id": disease["id"],
                "disease_name": disease["name"],
                "term_class": term_class,
                "therapeutic_area_ids": ";".join(a["id"] for a in disease.get("therapeuticAreas") or []),
                "therapeutic_area_names": ";".join(a["name"] for a in disease.get("therapeuticAreas") or []),
                "ancestor_ids": ";".join(disease.get("ancestors") or []),
                "overall_association_score": format(score, ".12g"),
                "direct_association": "true",
                "opentargets_data_version": release["response"]["dataVersion"],
            }
            association_rows.append(base)
            for item in row.get("datatypeScores") or []:
                datatype_rows.append({
                    "genage_id": gene["genage_id"],
                    "genage_symbol": gene["genage_symbol"],
                    "ensembl_gene_id": ensembl_id,
                    "disease_id": disease["id"],
                    "datatype_id": item["id"],
                    "datatype_score": format(float(item["score"]), ".12g"),
                })
            for item in row.get("datasourceScores") or []:
                datasource_rows.append({
                    "genage_id": gene["genage_id"],
                    "genage_symbol": gene["genage_symbol"],
                    "ensembl_gene_id": ensembl_id,
                    "disease_id": disease["id"],
                    "datasource_id": item["id"],
                    "datasource_score": format(float(item["score"]), ".12g"),
                })

        lm = longevity_by_gene.get(gene["genage_id"], [])
        summary_rows.append({
            "genage_id": gene["genage_id"],
            "genage_symbol": gene["genage_symbol"],
            "current_symbol": gene["current_symbol"],
            "entrez_id": gene["entrez_id"],
            "ensembl_gene_id": ensembl_id,
            "query_status": route_status,
            "direct_association_count": len(rows),
            "disease_or_clinical_condition_count": class_counts["disease_or_clinical_condition"],
            "phenotype_or_trait_count": class_counts["phenotype_or_trait"],
            "measurement_count": class_counts["measurement"],
            "biological_process_count": class_counts["biological_process"],
            "medical_procedure_count": class_counts["medical_procedure"],
            "nonhuman_or_other_count": class_counts["nonhuman_or_other"],
            "max_overall_association_score": format(max(scores), ".12g") if scores else "",
            "longevitymap_record_count": len({x["longevitymap_id"] for x in lm}),
            "longevitymap_association_labels": ";".join(sorted({x["association"] for x in lm})),
        })

    association_fields = [
        "genage_id", "genage_symbol", "current_symbol", "entrez_id", "ensembl_gene_id",
        "opentargets_target_id", "opentargets_target_symbol", "disease_id", "disease_name",
        "term_class", "therapeutic_area_ids", "therapeutic_area_names", "ancestor_ids",
        "overall_association_score", "direct_association", "opentargets_data_version",
    ]
    datatype_fields = ["genage_id", "genage_symbol", "ensembl_gene_id", "disease_id", "datatype_id", "datatype_score"]
    datasource_fields = ["genage_id", "genage_symbol", "ensembl_gene_id", "disease_id", "datasource_id", "datasource_score"]
    summary_fields = list(summary_rows[0])
    write_csv(OUT / "opentargets_gene_disease.csv", association_rows, association_fields)
    write_csv(OUT / "opentargets_datatype_scores.csv", datatype_rows, datatype_fields)
    write_csv(OUT / "opentargets_datasource_scores.csv", datasource_rows, datasource_fields)
    write_csv(OUT / "gene_disease_summary.csv", summary_rows, summary_fields)

    status_counts = Counter(r["query_status"] for r in summary_rows)
    term_counts = Counter(r["term_class"] for r in association_rows)
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "opentargets_release": release["response"],
        "direct_only": True,
        "page_size": PAGE_SIZE,
        "genage_gene_count": len(audit),
        "query_status_counts": dict(status_counts),
        "unique_gene_disease_pairs": len(association_rows),
        "term_class_counts": dict(term_counts),
        "datatype_score_rows": len(datatype_rows),
        "datasource_score_rows": len(datasource_rows),
        "genes_with_at_least_one_direct_association": sum(int(r["direct_association_count"]) > 0 for r in summary_rows),
        "pagination_reconciliation": {
            "genes_requiring_reconciliation": sum(bool(r["required"]) for r in reconciliation_audit),
            "genes_reconciled_to_declared_count": sum(
                bool(r["required"]) and r["final_unique"] == r["declared"] for r in reconciliation_audit
            ),
            "audit": reconciliation_audit,
        },
        "retrieval_utc_min": min(retrieval_times) if retrieval_times else None,
        "retrieval_utc_max": max(retrieval_times) if retrieval_times else None,
        "output_sha256": {
            name: sha256(OUT / name)
            for name in [
                "opentargets_gene_disease.csv", "opentargets_datatype_scores.csv",
                "opentargets_datasource_scores.csv", "gene_disease_summary.csv",
            ]
        },
        "interpretation_limit": (
            "Open Targets associations summarise heterogeneous evidence. Scores are not causal effect sizes "
            "or probabilities, and association with disease does not establish an effect on human lifespan."
        ),
    }
    (OUT / "opentargets_collection_summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
