#!/usr/bin/env python3
"""Provenance of Turkey-labelled public M. tuberculosis records and search for Turkish submitters.
Outputs data/mycotb_extended_out/turkey_provenance.csv and turkey_provenance.txt (offline)."""
import os, sys, re
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mycotb_analysis as A
BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))
OUT = BASE / "mycotb_extended_out"; OUT.mkdir(exist_ok=True)
_attr, samples, _ = A.load_and_clean()
strict = [s for s in samples if s["country"] == "Turkey"]
raw = pd.read_csv(A.SRC, sep="\t", dtype=str, low_memory=False).fillna("")
lines = [f"Strict samples labelled Turkey: {len(strict)}"]
rows = []
for s in strict:
    r = raw[raw.sample_accession == s["sample"]].iloc[0]
    rows.append({"sample": s["sample"], "study": s["study"], "center": r.center_name, "collection_year": s["cyear"] or "", "release_year": s["pyear"] or "", "study_title": r.study_title})
pd.DataFrame(rows).to_csv(OUT / "turkey_provenance.csv", index=False)
for r in rows: lines.append(f"  {r['sample']} | {r['study']} | {r['center'][:60]} | {r['study_title'][:90]}")
t = raw[raw.country.str.contains("turk", case=False)]
lines.append(f"All runs with Turkey country label (any taxon label): {len(t)} runs, {t.sample_accession.nunique()} samples, studies: {sorted(t.study_accession.unique())}")
pat = (r"kirikkale|kırıkkale|coban|çoban|akvuam|ankara|izmir|konya|selcuk|pamukkale|gaziantep|sakarya|bursa|uludag|trakya|samsun|kocaeli|baskent|yeditepe|bezmialem|acibadem|koc univ|sabanci|bilkent|middle east tech|tubitak|turkish|turkiye|türkiye|hacettepe|istanbul|marmara|dokuz eylul|ege univ|ataturk|erciyes|cukurova|akdeniz|ondokuz|gazi univ|aydin|adnan menderes|gebze|yildiz|medipol|biruni|umraniye|suat|sureyyapasa|yedikule")
hit = raw[raw.center_name.str.contains(pat, case=False, regex=True) | raw.study_title.str.contains(pat, case=False, regex=True) | raw.strain.str.contains(pat, case=False, regex=True)]
lines.append(f"Runs (of {len(raw)}) whose submitting centre, study title or strain contains a Turkish institution/place token: {len(hit)}")
(OUT / "turkey_provenance.txt").write_text("\n".join(lines), encoding="utf-8"); print("\n".join(lines))
