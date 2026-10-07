"""Cache primary publication metadata for the expanded review illustration."""
import json
import time
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from pathlib import Path

TASK = Path(__file__).resolve().parents[1]
OUT = TASK / "results/expanded/papers"
PMIDS = ["10921902", "15574588", "23312286", "25199830", "19302372",
         "18056421", "17912023", "19469880", "17204841", "20421943",
         "23939249", "12586694", "12456885", "19197346", "21108731"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for pmid in PMIDS:
        path = OUT / f"{pmid}.json"
        if not path.exists():
            query = urllib.parse.urlencode({"query": f"EXT_ID:{pmid} AND SRC:MED",
                                           "format": "json", "resultType": "core"})
            url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + query
            with urllib.request.urlopen(url, timeout=45) as response:
                data = json.load(response)
            assert data["hitCount"] == 1
            paper = data["resultList"]["result"][0]
            assert paper["id"] == pmid
            path.write_text(json.dumps({"retrieved": "2026-09-25", "source_url": url,
                                        "paper": paper}, indent=2), encoding="utf-8")
            time.sleep(.4)
        paper = json.loads(path.read_text(encoding="utf-8"))["paper"]
        print(pmid, paper["title"], paper.get("pmcid", ""), flush=True)
        print(paper.get("abstractText", "NO ABSTRACT"), flush=True)
        if pmid in {"10921902", "15574588"}:
            path = OUT / f"{pmid}.xml"
            if not path.exists():
                url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/"
                       + paper["pmcid"] + "/fullTextXML")
                try:
                    with urllib.request.urlopen(url, timeout=45) as response:
                        payload = response.read()
                except urllib.error.URLError as error:
                    print(f"Full text unavailable for {pmid}: {error}; use recorded web review.", flush=True)
                    continue
                root = ET.fromstring(payload)
                assert any(x.text == pmid for x in root.findall('.//article-id[@pub-id-type="pmid"]'))
                path.write_bytes(payload)
            root = ET.fromstring(path.read_bytes())
            for paragraph in root.findall(".//p"):
                content = "".join(paragraph.itertext())
                if any(term in content.lower() for term in ("overexpress", "extra copies", "forced expression")):
                    print(content, flush=True)


if __name__ == "__main__":
    main()
