#!/usr/bin/env python3
"""Reproducible log of the literature searches behind the novelty statement (read-only network queries).
PubMed via NCBI E-utilities; preprints (bioRxiv, medRxiv and others) via Europe PMC (SRC:PPR). Output: data/mycotb_extended_out/literature_search_log.csv"""
import csv, datetime, json, os, time, urllib.parse, urllib.request
from pathlib import Path
OUT = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data"))) / "mycotb_extended_out" / "literature_search_log.csv"
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PUBMED = {
 "PubMed 1": '(tuberculosis[ti] OR "Mycobacterium tuberculosis"[ti]) AND (genom*[ti] OR sequenc*[ti]) AND (represent*[ti] OR equit*[ti] OR bias*[ti] OR disparit*[ti] OR inequal*[ti] OR gaps[ti] OR coverage[ti])',
 "PubMed 2": '("Mycobacterium tuberculosis"[tiab]) AND ("public database*"[tiab] OR "public genome*"[tiab] OR "sequence read archive"[tiab] OR "European Nucleotide Archive"[tiab]) AND (incidence[tiab] OR burden[tiab]) AND (country[tiab] OR geographic*[tiab])',
 "PubMed 3": '("Mycobacterium tuberculosis"[tiab]) AND (genom*[tiab] OR "whole genome"[tiab]) AND ("sampling bias"[tiab] OR representativeness[tiab] OR "geographic bias"[tiab] OR underrepresent*[tiab] OR "under-represent*"[tiab])',
 "PubMed 4": '(tuberculosis[tiab]) AND ("genomic surveillance"[tiab]) AND (equity[tiab] OR inequit*[tiab] OR disparit*[tiab] OR gap[tiab] OR gaps[tiab])',
}
EPMC = {
 "Europe PMC preprints 1": '("Mycobacterium tuberculosis" AND (genome OR genomes OR sequencing) AND (representativeness OR "sampling bias" OR inequity OR disparity OR "public database" OR "public databases")) AND SRC:PPR',
 "Europe PMC preprints 2": 'TITLE:("tuberculosis" AND (genomes OR sequencing OR genomic) AND (representation OR representativeness OR equity OR bias OR disparities OR gaps)) AND SRC:PPR',
}
rows = []; today = datetime.date.today().isoformat()
for name, q in PUBMED.items():
    r = json.load(urllib.request.urlopen(EU + "esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=" + urllib.parse.quote(q)))["esearchresult"]
    rows.append({"date": today, "database": name, "query": q, "records": int(r["count"])}); time.sleep(.5)
for name, q in EPMC.items():
    try:
        r = json.load(urllib.request.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=1&query=" + urllib.parse.quote(q)))
        rows.append({"date": today, "database": name, "query": q, "records": int(r["hitCount"])})
    except Exception as e:  # noqa: BLE001
        rows.append({"date": today, "database": name, "query": q, "records": f"error: {e}"})
    time.sleep(.5)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["date", "database", "query", "records"]); w.writeheader(); w.writerows(rows)
for r in rows: print(r["database"], r["records"])
