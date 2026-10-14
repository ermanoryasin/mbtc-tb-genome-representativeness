#!/usr/bin/env python3
"""MycoTB v6 - extended robustness: lag-free regional analysis and study-dependence checks.

Offline (uses the cleaned samples from mycotb_analysis, WHO table, regression dataset).
Outputs go to data/mycotb_extended_out/.

A. WHO region x release-year cohort (removes publication-lag exposure of Fig 1c / regions.csv).
B. Study dependence (samples are not independent isolates):
   B1 cap on genomes per study within a country (k = 25, 100);
   B2 study-level cluster bootstrap CIs for Gini, genomes/1,000 cases and income-group RI;
   B3 Poisson models with (i) HC1 and region-clustered SE, (ii) capped counts, (iii) number of distinct studies.
Requires numpy and pandas.
"""
import csv
import datetime
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mycotb_analysis as A  # noqa: E402
import mycotb_robustness as R  # noqa: E402

MBTC_BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))
OUT = Path(os.environ.get("MYCO_OUT4", str(MBTC_BASE / "mycotb_extended_out")))
DS = MBTC_BASE / "mycotb_regression_out" / "analysis_dataset.csv"
P0, P1 = 2010, 2024
BOOT_N = 1000
SEED = 20261014
LINES = []


def out(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


def gini_vec(gen, inc):
    """Gini of the Lorenz curve of genome share vs burden share (countries ordered by representation).
    Same definition as mycotb_analysis.metrics; gen and inc are 1-D arrays over countries with inc > 0."""
    tw, ti = gen.sum(), inc.sum()
    if tw == 0:
        return np.nan
    order = np.argsort(gen / inc, kind="stable")
    pi, pw = inc[order] / ti, gen[order] / tw
    cw = np.cumsum(pw) - pw
    return 1 - 2 * np.sum(pi * (cw + pw / 2))


def income_ri(gen, inc, income_idx, n_groups=4):
    """Representation index by income group from arrays; income_idx -1 = unmatched (ignored)."""
    ok = income_idx >= 0
    tw, ti = gen[ok].sum(), inc[ok].sum()
    ri = np.full(n_groups, np.nan)
    for g in range(n_groups):
        m = income_idx == g
        if m.any() and tw > 0:
            ri[g] = (gen[m].sum() / tw) / (inc[m].sum() / ti)
    return ri


# ----------------------------------------------------------------- Poisson (IRLS) with robust / clustered SE
def poisson_irls(X, y, offset, iters=100):
    b = np.zeros(X.shape[1])
    b[0] = np.log(max(y.sum(), 1) / np.exp(offset).sum())
    for _ in range(iters):
        eta = np.clip(X @ b + offset, -30, 30)
        mu = np.exp(eta)
        z = eta - offset + (y - mu) / mu
        XtW = X.T * mu
        nb = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(nb - b)) < 1e-10:
            b = nb
            break
        b = nb
    mu = np.exp(np.clip(X @ b + offset, -30, 30))
    return b, mu


def sandwich(X, y, mu, cluster=None):
    n, k = X.shape
    bread = np.linalg.inv((X.T * mu) @ X)
    u = X * (y - mu)[:, None]
    if cluster is None:
        meat = u.T @ u
        corr = n / (n - k)  # HC1
    else:
        meat = np.zeros((k, k))
        groups = pd.Series(cluster).astype("category").cat.codes.to_numpy()
        G = groups.max() + 1
        for g in range(G):
            s = u[groups == g].sum(axis=0)
            meat += np.outer(s, s)
        corr = (G / (G - 1)) * ((n - 1) / (n - k))  # CR1
    return np.sqrt(np.diag(bread @ meat @ bread * corr))


def fit_row(label, X, y, off, names, cluster=None, dist="norm"):
    from scipy import stats
    b, mu = poisson_irls(X, y, off)
    se = sandwich(X, y, mu, cluster)
    G = None if cluster is None else len(set(cluster))
    out(f"  [{label}]  n={len(y)}, dispersion={((y - mu) ** 2 / mu).sum() / (len(y) - X.shape[1]):.1f}"
        + (f", clusters={G}" if G else ""))
    rows = []
    for j, nm in enumerate(names):
        if nm == "const":
            continue
        crit = stats.t.ppf(0.975, G - 1) if (G and dist == "t") else 1.959964
        lo, hi = b[j] - crit * se[j], b[j] + crit * se[j]
        p = 2 * (stats.t.sf(abs(b[j] / se[j]), G - 1) if (G and dist == "t") else stats.norm.sf(abs(b[j] / se[j])))
        out(f"     {nm:9s} RR={np.exp(b[j]):.2f} ({np.exp(lo):.2f}-{np.exp(hi):.2f}) p={p:.3g}")
        rows.append({"model": label, "term": nm, "rr": round(float(np.exp(b[j])), 3),
                     "ci_lo": round(float(np.exp(lo)), 3), "ci_hi": round(float(np.exp(hi)), 3),
                     "p": float(f"{p:.3g}"), "n": len(y)})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out(f"MycoTB v6 extended robustness - {datetime.datetime.now():%Y-%m-%d %H:%M}")
    _attr, samples, _ = A.load_and_clean()
    inc, _lo, _hi, region = R.load_who_full()
    out(f"Analysis unit: {len(samples)} clean WGS samples")

    # ============================================ A. region x release cohort
    out("")
    out("=== A. WHO region x release-year cohort (burden of the same calendar window) ===")
    rows = []
    for (y0, y1) in R.REL_WINDOWS:
        gen = R.gen_counts(samples, "p", y0, y1)
        inc_w = R.psum(inc, y0, y1)
        ti = sum(inc_w.values())
        tw = sum(gen.get(c, 0) for c in inc_w)
        for reg in ["AFR", "AMR", "EMR", "EUR", "SEA", "WPR"]:
            cs = [c for c in inc_w if region.get(c) == reg]
            b = sum(inc_w[c] for c in cs)
            g = sum(gen.get(c, 0) for c in cs)
            bs, gs = 100 * b / ti, 100 * g / tw
            rows.append({"release_window": f"{y0}-{y1}", "who_region": reg, "burden_share_pct": round(bs, 2),
                         "genome_share_pct": round(gs, 2), "representation_index": round(gs / bs, 3),
                         "genomes_per_1000": round(1000 * g / b, 3), "genomes": g})
            out(f"  {y0}-{y1} {reg}: burden {bs:5.1f}% genomes {gs:5.1f}% RI {gs / bs:6.2f}  {1000 * g / b:7.3f}/1000")
    pd.DataFrame(rows).to_csv(OUT / "regions_release.csv", index=False)

    # ============================================ B1. cap per study within country
    out("")
    out("=== B1. Cap on genomes per study within country (2010-2024, collection year) ===")
    sc = defaultdict(Counter)  # country -> study -> n
    for s in samples:
        if A.is_missing_country(s["country"]) or not s["cyear"] or not (P0 <= s["cyear"] <= P1):
            continue
        sc[s["country"]][s["study"]] += 1
    inc_main = R.psum(inc, P0, P1)
    countries = sorted(inc_main)
    inc_arr = np.array([inc_main[c] for c in countries])
    ds = pd.read_csv(DS)
    inc_map = {"LIC": 0, "LMC": 1, "UMC": 2, "HIC": 3}
    income_of = dict(zip(ds.country, ds.income))
    income_idx = np.array([inc_map.get(income_of.get(c), -1) for c in countries])
    rows = []
    for label, k in [("no cap", None), ("cap 100 per study", 100), ("cap 25 per study", 25)]:
        gen = np.array([sum(min(n, k) if k else n for n in sc[c].values()) if c in sc else 0 for c in countries], float)
        ri = income_ri(gen, inc_arr, income_idx)
        top = np.argsort(-inc_arr)[:30]
        rows.append({"scenario": label, "genomes": int(gen.sum()), "genomes_per_1000": round(1000 * gen.sum() / inc_arr.sum(), 3),
                     "gini": round(gini_vec(gen, inc_arr), 3), "zero_countries": int((gen == 0).sum()),
                     "top30_genome_pct": round(100 * gen[top].sum() / gen.sum(), 1),
                     "RI_LIC": round(ri[0], 2), "RI_LMC": round(ri[1], 2), "RI_UMC": round(ri[2], 2), "RI_HIC": round(ri[3], 2)})
        out(f"  {label}: {rows[-1]}")
    pd.DataFrame(rows).to_csv(OUT / "study_cap_sensitivity.csv", index=False)

    # ============================================ B2. study-level cluster bootstrap
    out("")
    out(f"=== B2. Study-level cluster bootstrap ({BOOT_N} resamples of studies) ===")
    studies = sorted({st for c in sc for st in sc[c]})
    sidx = {s: i for i, s in enumerate(studies)}
    cidx = {c: i for i, c in enumerate(countries)}
    M = np.zeros((len(studies), len(countries)))
    for c, cnt in sc.items():
        if c in cidx:
            for st, n in cnt.items():
                M[sidx[st], cidx[c]] = n
    rng = np.random.default_rng(SEED)
    ginis, per1000, ris, top30 = [], [], [], []
    top = np.argsort(-inc_arr)[:30]
    for _ in range(BOOT_N):
        w = np.bincount(rng.integers(0, len(studies), len(studies)), minlength=len(studies))
        gen = w @ M
        if gen.sum() == 0:
            continue
        ginis.append(gini_vec(gen, inc_arr))
        per1000.append(1000 * gen.sum() / inc_arr.sum())
        ris.append(income_ri(gen, inc_arr, income_idx))
        top30.append(100 * gen[top].sum() / gen.sum())
    ris = np.array(ris)
    q = lambda v: (np.nanpercentile(v, 2.5), np.nanpercentile(v, 97.5))  # noqa: E731
    base_gen = M.sum(axis=0)
    res = [("gini", gini_vec(base_gen, inc_arr), q(ginis)),
           ("genomes_per_1000", 1000 * base_gen.sum() / inc_arr.sum(), q(per1000)),
           ("top30_genome_pct", 100 * base_gen[top].sum() / base_gen.sum(), q(top30))]
    base_ri = income_ri(base_gen, inc_arr, income_idx)
    for g, nm in enumerate(["LIC", "LMC", "UMC", "HIC"]):
        res.append((f"RI_{nm}", base_ri[g], q(ris[:, g])))
    tab = []
    for nm, pt, (lo_, hi_) in res:
        out(f"  {nm:18s} {pt:8.3f}  95% CI {lo_:.3f}-{hi_:.3f}")
        tab.append({"metric": nm, "point": round(float(pt), 3), "ci_lo": round(float(lo_), 3), "ci_hi": round(float(hi_), 3)})
    pd.DataFrame(tab).to_csv(OUT / "study_bootstrap.csv", index=False)

    # ============================================ B3. regressions
    out("")
    out("=== B3. Poisson models: SE type and outcome definition (P3 specification) ===")
    d = ds.copy()
    d = d[d.gdp.notna() & d.artpc.notna() & d.rate.notna() & d.income.notna() & d.region.notna()]
    d = d[(d.gdp > 0) & (d.artpc > 0) & (d.rate > 0) & (d.cases > 0)].reset_index(drop=True)
    out(f"  analysis set n = {len(d)} (published P3: 187)")
    X = np.column_stack([np.ones(len(d)), np.log2(d.gdp), np.log2(d.rate), np.log2(d.artpc)])
    names = ["const", "log2 GDP", "log2 rate", "log2 art/M"]
    off = np.log(d.cases.to_numpy())
    obs = d.obs.to_numpy(float)
    cap25 = np.array([sum(min(n, 25) for n in sc[c].values()) if c in sc else 0 for c in d.country], float)
    nstud = np.array([len(sc[c]) if c in sc else 0 for c in d.country], float)
    clu = d.region.to_numpy()
    allrows = []
    allrows += fit_row("P3 HC1 (published)", X, obs, off, names)
    allrows += fit_row("P3 region-clustered CR1 (normal)", X, obs, off, names, cluster=clu)
    allrows += fit_row("P3 region-clustered CR1 (t, G-1 df)", X, obs, off, names, cluster=clu, dist="t")
    allrows += fit_row("Outcome: genomes capped at 25/study, HC1", X, cap25, off, names)
    allrows += fit_row("Outcome: number of distinct studies, HC1", X, nstud, off, names)
    pd.DataFrame(allrows).to_csv(OUT / "study_dependence_models.csv", index=False)

    (OUT / "summary.txt").write_text("\n".join(LINES), encoding="utf-8")
    out(f"\nOutputs: {OUT}")


if __name__ == "__main__":
    main()
