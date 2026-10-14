#!/usr/bin/env python3
"""Export the cleaned analysis unit (one row per BioSample) so results can be audited without the 119 MB raw ENA table.
Columns: sample, country (normalised), collection_year, release_year, study. Public INSDC metadata only."""
import gzip, csv, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mycotb_analysis as A
MBTC_BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))
OUT = MBTC_BASE / "clean_samples" / "mtb_clean_samples.csv.gz"
OUT.parent.mkdir(parents=True, exist_ok=True)
_attr, samples, _ = A.load_and_clean()
with gzip.open(OUT, "wt", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample", "country", "collection_year", "release_year", "study"])
    for s in samples:
        w.writerow([s["sample"], s["country"], s["cyear"] or "", s["pyear"] or "", s["study"]])
print(len(samples), "samples ->", OUT, f"{OUT.stat().st_size/1e6:.1f} MB")
