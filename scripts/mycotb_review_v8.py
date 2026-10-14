#!/usr/bin/env python3
"""MycoTB v8 - responses to external review.

1. Flow table of the analysis set (every drop with one reason).
2. Missing-country scenarios: Gini, representation index (RI) by income group and zero-genome countries for six scenarios
   (S0 excluded, S1 proportional to observed, S2 proportional to burden, S3 low-coverage countries, S4 study-informed, S5 adverse).
3. Probe of missingness: share of samples without country by the (modal) income group of the study they belong to.
4. Publication-lag sensitivity: thresholds of 3, 5 and 7 years, equal 3-year windows, annual series; income composition of fast and slow depositors.
5. Year-specific World Bank income classification (Our World in Data copy of the World Bank series) versus the current classification.
Offline except for the classification file in data/mycotb_extended_out/external/. Outputs in data/mycotb_extended_out/.
"""
import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mycotb_analysis as A  # noqa: E402
import mycotb_robustness as R  # noqa: E402

BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))
OUT = BASE / "mycotb_extended_out"
EXT = OUT / "external" / "wb_income_groups_by_year.csv"
DS = BASE / "mycotb_regression_out" / "analysis_dataset.csv"
P0, P1 = 2010, 2024
LINES = []


def out(s=""):
    print(s, flush=True)
    LINES.append(s)


def ri_by_income(gen, inc, income_of):
    """RI per income group from dicts gen (genomes, may be fractional) and inc (burden); income_of: country -> group or None."""
    gs, bs = defaultdict(float), defaultdict(float)
    for c, v in inc.items():
        g = income_of.get(c)
        if g:
            bs[g] += v; gs[g] += gen.get(c, 0)
    tg, tb = sum(gs.values()), sum(bs.values())
    return {g: (gs[g] / tg) / (bs[g] / tb) if bs[g] > 0 and tg > 0 else float("nan") for g in ["LIC", "LMC", "UMC", "HIC"]}


def main():
    attr, samples, _ = A.load_and_clean()
    inc, _lo, _hi, region = R.load_who_full()
    ds = pd.read_csv(DS)
    iso_of = dict(zip(ds.country, ds.iso3))
    income_now = {c: (g if isinstance(g, str) else None) for c, g in zip(ds.country, ds.income)}
    main_inc = R.psum(inc, P0, P1)
    out(f"MycoTB v8 review responses; samples {len(samples)}")

    # ------------------------------------------------ 1. flow table
    known = [s for s in samples if not A.is_missing_country(s["country"])]
    known_dated = [s for s in known if s["cyear"] and P0 <= s["cyear"] <= P1]
    matched = [s for s in known_dated if s["country"] in main_inc]
    unmatched = Counter(s["country"] for s in known_dated if s["country"] not in main_inc)
    n_reg = len(ds[ds.gdp.notna() & ds.artpc.notna() & ds.rate.notna() & ds.income.notna() & ds.region.notna() & (ds.gdp > 0) & (ds.artpc > 0) & (ds.rate > 0) & (ds.cases > 0)])
    flow = [
        ("ENA runs, taxon tree 1773 (BioSample-level: %d)" % attr[0][2], attr[0][1], "retrieved 14 October 2026; NCBI Entrez count of 378990 records was higher by 19385 and could not be reconciled"),
        ("after strict M. tuberculosis label", attr[1][1], "excluded variants and complex members (bovis, africanum, caprae, BCG and others)"),
        ("after excluding laboratory reference strains", attr[2][1], "H37Rv, H37Ra, Erdman, HN878, CDC1551, ATCC"),
        ("after whole-genome, genomic-source runs", attr[3][1], "excluded RNA-seq, Tn-seq and amplicon assays"),
        ("BioSamples (analysis unit)", len(samples), "runs collapsed to BioSample"),
        ("with a single known country", len(known), "country missing or multiple for %d (%.1f%%)" % (len(samples) - len(known), 100 * (len(samples) - len(known)) / len(samples))),
        ("with known country and collection year 2010-2024", len(known_dated), "collection year missing or outside 2010-2024 for %d" % (len(known) - len(known_dated))),
        ("matched to a WHO reporting country (denominator of all representation statistics)", len(matched), "unmatched: " + ", ".join(f"{c} {n}" for c, n in unmatched.most_common(4)) + ", other %d" % (sum(unmatched.values()) - sum(n for _, n in unmatched.most_common(4)))),
        ("countries with WHO incidence, 2010-2024", len(main_inc), "of which %d with a World Bank income group (5 small territories unclassified, all without genomes)" % ds.income.notna().sum()),
        ("countries in the regression models", n_reg, "complete GDP, articles per million, incidence rate, region and income"),
    ]
    pd.DataFrame(flow, columns=["step", "count", "reason_or_note"]).to_csv(OUT / "flow_table.csv", index=False)
    out("\n=== 1. Flow table ===")
    for r in flow: out(f"  {r[1]:>8}  {r[0]}  [{r[2]}]")

    # ------------------------------------------------ 2. scenarios
    out("\n=== 2. Missing-country scenarios (2010-2024) ===")
    miss = [s for s in samples if A.is_missing_country(s["country"]) and s["cyear"] and P0 <= s["cyear"] <= P1]
    nm = len(miss)
    base = Counter(s["country"] for s in matched)
    ti = sum(main_inc.values()); tw = sum(base.get(c, 0) for c in main_inc)
    scen = {"S0 missing excluded (main analysis)": dict(base)}
    k = (tw + nm) / tw
    scen["S1 missing proportional to observed genomes"] = {c: base.get(c, 0) * k for c in main_inc}
    scen["S2 missing proportional to burden (most favourable to equality)"] = {c: base.get(c, 0) + nm * main_inc[c] / ti for c in main_inc}
    low = [c for c in main_inc if 1000.0 * base.get(c, 0) / main_inc[c] < A.THR]
    li = sum(main_inc[c] for c in low)
    g3 = dict(base)
    for c in low: g3[c] = g3.get(c, 0) + nm * main_inc[c] / li
    scen["S3 missing to low-coverage countries by burden"] = g3
    # S4 study-informed: allocate each missing sample to the known-country distribution of its own study
    by_study = defaultdict(Counter)
    for s in known:
        if s["country"] in main_inc: by_study[s["study"]][s["country"]] += 1
    g4 = dict(base); allocated = 0
    for s in miss:
        cnt = by_study.get(s["study"])
        if cnt:
            tot = sum(cnt.values()); allocated += 1
            for c, v in cnt.items(): g4[c] = g4.get(c, 0) + v / tot
    resid = nm - allocated
    kk = (sum(g4.get(c, 0) for c in main_inc) + resid) / sum(g4.get(c, 0) for c in main_inc)
    scen["S4 missing follow their own study's countries (rest proportional)"] = {c: g4.get(c, 0) * kk for c in main_inc}
    high = [c for c in main_inc if 1000.0 * base.get(c, 0) / main_inc[c] >= 10]
    hw = sum(base.get(c, 0) for c in high)
    g5 = dict(base)
    for c in high: g5[c] = g5.get(c, 0) + nm * base.get(c, 0) / hw
    scen["S5 missing all to high-coverage countries (least favourable)"] = g5
    rows = []
    for name, g in scen.items():
        m = A.metrics(g, main_inc); ri = ri_by_income(g, main_inc, income_now)
        rows.append({"scenario": name, "genomes": round(sum(g.get(c, 0) for c in main_inc)), "genomes_per_1000": round(m["per1000"], 2), "gini": round(m["gini"], 3),
                     "zero_genome_countries": m["zero_countries"], **{f"RI_{k_}": round(v, 2) for k_, v in ri.items()}})
        out(f"  {name}: {rows[-1]}")
    out(f"  (samples with valid 2010-2024 year but missing country: {nm}; allocated to a study with known countries: {allocated}, residual {resid})")
    pd.DataFrame(rows).to_csv(OUT / "missing_country_scenarios.csv", index=False)

    # ------------------------------------------------ 3. probe of missingness by study income
    out("\n=== 3. Share of samples without country, by income group of the study ===")
    study_mod = {}
    for st, cnt in by_study.items():
        c = cnt.most_common(1)[0][0]; study_mod[st] = income_now.get(c)
    tot_by, miss_by = Counter(), Counter()
    for s in samples:
        g = study_mod.get(s["study"])
        if g is None: continue
        tot_by[g] += 1
        if A.is_missing_country(s["country"]): miss_by[g] += 1
    probe = [{"study_income_group": g, "samples": tot_by[g], "without_country": miss_by[g], "share_without_country_pct": round(100 * miss_by[g] / tot_by[g], 1)} for g in ["LIC", "LMC", "UMC", "HIC"] if tot_by[g]]
    no_known = [s for s in samples if s["study"] not in by_study]
    probe.append({"study_income_group": "studies with no sample of known country", "samples": len(no_known), "without_country": sum(1 for s in no_known if A.is_missing_country(s["country"])), "share_without_country_pct": round(100 * len(no_known) / len(samples), 1)})
    pd.DataFrame(probe).to_csv(OUT / "missing_country_probe.csv", index=False)
    for r in probe: out(f"  {r}")

    # ------------------------------------------------ 4. lag sensitivity
    out("\n=== 4. Publication-lag sensitivity ===")
    WINDOWS = [(2010, 2012), (2013, 2015), (2016, 2018), (2019, 2021)]
    rows = []
    for kth in (3, 5, 7):
        for (a, b) in WINDOWS:
            if b + kth > 2025: continue
            flt = lambda s, kth=kth: s["pyear"] is not None and s["cyear"] is not None and 0 <= s["pyear"] - s["cyear"] <= kth  # noqa: E731
            gen = R.gen_counts(samples, "c", a, b, filt=flt); ic = R.psum(inc, a, b); m = A.metrics(gen, ic)
            rows.append({"lag_threshold_years": kth, "window": f"{a}-{b}", "genomes": m["genomes"], "genomes_per_1000": round(m["per1000"], 3), "gini": round(m["gini"], 3), "zero_genome_countries": m["zero_countries"]})
    lagdf = pd.DataFrame(rows); lagdf.to_csv(OUT / "lag_threshold_windows.csv", index=False)
    out(lagdf.to_string(index=False))
    ann = []
    for y in range(2010, 2023):
        flt = lambda s: s["pyear"] is not None and s["cyear"] is not None and 0 <= s["pyear"] - s["cyear"] <= 3  # noqa: E731
        gen = R.gen_counts(samples, "c", y, y, filt=flt); m = A.metrics(gen, R.psum(inc, y, y))
        ann.append({"collection_year": y, "genomes": m["genomes"], "gini": round(m["gini"], 3), "genomes_per_1000": round(m["per1000"], 3)})
    pd.DataFrame(ann).to_csv(OUT / "lag_annual_k3.csv", index=False)
    out("annual (k=3): " + ", ".join(f"{r['collection_year']}:{r['gini']}" for r in ann))
    # composition of fast vs slow depositors (collections 2010-2018, fully observable up to 7 years)
    comp = []
    for label, flt in [("fast (lag <= 3 y)", lambda s: s["pyear"] and s["cyear"] and 0 <= s["pyear"] - s["cyear"] <= 3),
                       ("slow (lag > 3 y)", lambda s: s["pyear"] and s["cyear"] and s["pyear"] - s["cyear"] > 3),
                       ("all", lambda s: s["pyear"] and s["cyear"] and s["pyear"] - s["cyear"] >= 0)]:
        gen = R.gen_counts(samples, "c", 2010, 2018, filt=flt); ic = R.psum(inc, 2010, 2018); m = A.metrics(gen, ic)
        tot = sum(gen.get(c, 0) for c in ic if income_now.get(c))
        sh = {g: 100 * sum(v for c, v in gen.items() if income_now.get(c) == g) / tot for g in ["LIC", "LMC", "UMC", "HIC"]}
        ri = ri_by_income(gen, ic, income_now)
        comp.append({"depositor_class": label, "genomes": m["genomes"], "gini": round(m["gini"], 3), **{f"share_{g}_pct": round(v, 1) for g, v in sh.items()}, **{f"RI_{g}": round(v, 2) for g, v in ri.items()}})
    pd.DataFrame(comp).to_csv(OUT / "lag_fast_slow_composition.csv", index=False)
    for r in comp: out("  " + str(r))

    # ------------------------------------------------ 5. year-specific income classification
    out("\n=== 5. Year-specific World Bank income classification ===")
    if EXT.exists():
        cl = pd.read_csv(EXT)
        mapg = {"Low-income countries": "LIC", "Lower-middle-income countries": "LMC", "Upper-middle-income countries": "UMC", "High-income countries": "HIC"}
        cl["g"] = cl.classification.map(mapg)
        by_iso_year = {(r.code, int(r.year)): r.g for r in cl.itertuples() if isinstance(r.g, str)}
        rows = []
        for label, (a, b) in [("2010-2014", (2010, 2014)), ("2015-2019", (2015, 2019)), ("2020-2024", (2020, 2024)), ("2010-2024", (2010, 2024))]:
            gs, bs = defaultdict(float), defaultdict(float)
            gs_c, bs_c = defaultdict(float), defaultdict(float)
            for s in samples:
                if A.is_missing_country(s["country"]) or not s["cyear"] or not (a <= s["cyear"] <= b) or s["country"] not in inc: continue
                g = by_iso_year.get((iso_of.get(s["country"]), s["cyear"]))
                if g: gs[g] += 1
                gc = income_now.get(s["country"])
                if gc: gs_c[gc] += 1
            for c, yy in inc.items():
                for y, v in yy.items():
                    if a <= y <= b:
                        g = by_iso_year.get((iso_of.get(c), y))
                        if g: bs[g] += v
                        gc = income_now.get(c)
                        if gc: bs_c[gc] += v
            tg, tb, tgc, tbc = sum(gs.values()), sum(bs.values()), sum(gs_c.values()), sum(bs_c.values())
            for g in ["LIC", "LMC", "UMC", "HIC"]:
                rows.append({"period": label, "income_group": g, "RI_year_specific": round((gs[g] / tg) / (bs[g] / tb), 3), "RI_current_classification": round((gs_c[g] / tgc) / (bs_c[g] / tbc), 3),
                             "burden_share_year_specific_pct": round(100 * bs[g] / tb, 2), "genome_share_year_specific_pct": round(100 * gs[g] / tg, 2)})
        yd = pd.DataFrame(rows); yd.to_csv(OUT / "income_year_specific.csv", index=False)
        out(yd.to_string(index=False))
    else:
        out("  classification file missing")
    (OUT / "review_v8_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
