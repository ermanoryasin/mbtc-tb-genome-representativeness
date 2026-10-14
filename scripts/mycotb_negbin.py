#!/usr/bin/env python3
"""MycoTB v7 - negative binomial (NB2) and quasi-Poisson sensitivity analysis of the country-level genome-count models.

Why: the Poisson model P3 has a Pearson dispersion of about 11487, so interval estimates from Poisson are approximate.
Models fitted (same sample as P3, n = 187; offset = log cumulative incident cases):
  P3  : log2 GDP, log2 incidence rate, log2 articles per million
  P4  : P3 + log2 cumulative cases (case volume)
Estimators: Poisson (HC1), quasi-Poisson (Poisson coefficients, SE scaled by sqrt(Pearson dispersion)), NB2 by maximum likelihood
with alpha estimated jointly (ML standard errors and HC1 sandwich with alpha fixed at its estimate).
Turkey: expected genomes with a 95% interval (coefficient draws from the sandwich covariance) and the probability
of observing at most the observed count under the NB2 model.
Offline; outputs in data/mycotb_extended_out/.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))
DS = BASE / "mycotb_regression_out" / "analysis_dataset.csv"
OUT = BASE / "mycotb_extended_out"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20261014)
LINES = []


def out(s=""):
    print(s, flush=True)
    LINES.append(s)


def poisson_fit(X, y, off):
    b = np.zeros(X.shape[1]); b[0] = np.log(y.sum() / np.exp(off).sum())
    for _ in range(200):
        mu = np.exp(np.clip(X @ b + off, -30, 30))
        z = X @ b + (y - mu) / mu
        XtW = X.T * mu
        nb = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(nb - b)) < 1e-10:
            b = nb; break
        b = nb
    return b, np.exp(np.clip(X @ b + off, -30, 30))


def nb_ll(params, X, y, off):
    b, la = params[:-1], params[-1]
    a = np.exp(la)
    mu = np.exp(np.clip(X @ b + off, -30, 30))
    r = 1.0 / a
    return -np.sum(special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1)
                   + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))


def nb_fit(X, y, off):
    b0, mu0 = poisson_fit(X, y, off)
    a0 = max(np.sum(((y - mu0) ** 2 - mu0) / mu0 ** 2) / len(y), 0.1)
    best = None
    for start in (np.log(a0), 0.0, 2.0):
        res = optimize.minimize(nb_ll, np.r_[b0, start], args=(X, y, off), method="BFGS", options={"maxiter": 2000, "gtol": 1e-7})
        if best is None or res.fun < best.fun: best = res
    b, a = best.x[:-1], float(np.exp(best.x[-1]))
    mu = np.exp(np.clip(X @ b + off, -30, 30))
    k = X.shape[1]
    # ML covariance for beta (alpha fixed): (X' W X)^-1, W = mu/(1+a mu)
    W = mu / (1 + a * mu)
    bread = np.linalg.inv((X.T * W) @ X)
    # numerical full Hessian for ML SE including alpha
    H = optimize.approx_fprime(best.x, lambda p: optimize.approx_fprime(p, nb_ll, 1e-5, X, y, off)[0], 1e-4) if False else None
    s = X * ((y - mu) / (1 + a * mu))[:, None]
    meat = s.T @ s * (len(y) / (len(y) - k))
    cov_rob = bread @ meat @ bread
    return b, a, mu, bread, cov_rob, -best.fun


def row(label, term, b, se, n, extra=""):
    lo, hi = b - 1.959964 * se, b + 1.959964 * se
    p = 2 * stats.norm.sf(abs(b / se))
    return {"model": label, "term": term, "rr": round(float(np.exp(b)), 3), "ci_lo": round(float(np.exp(lo)), 3),
            "ci_hi": round(float(np.exp(hi)), 3), "p": float(f"{p:.3g}"), "n": n}


def main():
    ds = pd.read_csv(DS)
    d = ds[ds.gdp.notna() & ds.artpc.notna() & ds.rate.notna() & ds.income.notna() & ds.region.notna()]
    d = d[(d.gdp > 0) & (d.artpc > 0) & (d.rate > 0) & (d.cases > 0)].reset_index(drop=True)
    y = d.obs.to_numpy(float); off = np.log(d.cases.to_numpy())
    n = len(d)
    out(f"MycoTB v7 negative binomial sensitivity; n = {n}")
    base_cols = [np.log2(d.gdp), np.log2(d.rate), np.log2(d.artpc)]
    specs = {"P3": (base_cols, ["log2 GDP", "log2 rate", "log2 art/M"]),
             "P4": (base_cols + [np.log2(d.cases)], ["log2 GDP", "log2 rate", "log2 art/M", "log2 cases"])}
    rows, turk = [], []
    ti = d.index[d.country == "Turkey"][0]
    for name, (cols, terms) in specs.items():
        X = np.column_stack([np.ones(n)] + cols)
        b, mu = poisson_fit(X, y, off)
        k = X.shape[1]
        disp = float(np.sum((y - mu) ** 2 / mu) / (n - k))
        bread = np.linalg.inv((X.T * mu) @ X)
        u = X * (y - mu)[:, None]
        cov_hc1 = bread @ (u.T @ u) @ bread * n / (n - k)
        out(f"\n[{name}] Pearson dispersion (Poisson) = {disp:.1f}")
        for j, t in enumerate(terms, start=1):
            rows.append(row(f"{name} Poisson HC1", t, b[j], np.sqrt(cov_hc1[j, j]), n))
            rows.append(row(f"{name} quasi-Poisson", t, b[j], np.sqrt(cov_hc1[j, j] * 0 + bread[j, j] * disp), n))
        bn, a, mun, bread_n, cov_rob, ll = nb_fit(X, y, off)
        out(f"[{name}] NB2 alpha = {a:.2f}  (log-likelihood {ll:.1f})")
        for j, t in enumerate(terms, start=1):
            rows.append(row(f"{name} NB2 ML", t, bn[j], np.sqrt(bread_n[j, j]), n))
            rows.append(row(f"{name} NB2 robust", t, bn[j], np.sqrt(cov_rob[j, j]), n))
        for r_ in rows:
            if r_["model"].startswith(name) and r_["term"] == "log2 GDP":
                out(f"   {r_['model']:22s} GDP RR {r_['rr']:.2f} ({r_['ci_lo']:.2f}-{r_['ci_hi']:.2f})")
        for r_ in rows:
            if r_["model"].startswith(name) and r_["term"] == "log2 rate":
                out(f"   {r_['model']:22s} rate RR {r_['rr']:.2f} ({r_['ci_lo']:.2f}-{r_['ci_hi']:.2f})")
        # Turkey: expected under Poisson and NB2, and probability of <= observed
        draws = RNG.multivariate_normal(bn, cov_rob, size=20000)
        mu_draw = np.exp(draws @ X[ti] + off[ti])
        r_ = 1.0 / a
        p_le = stats.nbinom.cdf(int(y[ti]), r_, r_ / (r_ + mun[ti]))
        # probability averaged over coefficient uncertainty
        p_draw = stats.nbinom.cdf(int(y[ti]), r_, r_ / (r_ + mu_draw)).mean()
        turk.append({"model": name, "observed": int(y[ti]), "expected_poisson": round(float(mu[ti]), 0),
                     "expected_nb2": round(float(mun[ti]), 0), "expected_nb2_lo": round(float(np.percentile(mu_draw, 2.5)), 0),
                     "expected_nb2_hi": round(float(np.percentile(mu_draw, 97.5)), 0), "alpha": round(a, 2),
                     "p_obs_le_under_nb2": float(f"{p_le:.3g}"), "p_obs_le_averaged": float(f"{p_draw:.3g}")})
        # shortfall rank under NB2 (expected - observed) among 187
        sf = pd.Series(mun - y, index=d.country).sort_values(ascending=False)
        turk[-1]["shortfall_rank_nb2"] = int(list(sf.index).index("Turkey") + 1)
        sfp = pd.Series(mu - y, index=d.country).sort_values(ascending=False)
        turk[-1]["shortfall_rank_poisson"] = int(list(sfp.index).index("Turkey") + 1)
        out(f"   Turkey: observed {int(y[ti])}, expected Poisson {mu[ti]:.0f}, NB2 {mun[ti]:.0f} "
            f"(95% {np.percentile(mu_draw, 2.5):.0f}-{np.percentile(mu_draw, 97.5):.0f}); P(Y<={int(y[ti])}) under NB2 = {p_le:.3g}; "
            f"rank by shortfall NB2 {turk[-1]['shortfall_rank_nb2']} vs Poisson {turk[-1]['shortfall_rank_poisson']}")
    # ---- additional specifications: Poisson (HC1) vs NB2 (robust) ----
    L = np.log2
    reg = pd.get_dummies(d.region).drop(columns=["AFR"]).astype(float)
    inc = pd.get_dummies(d.income).drop(columns=["LIC"]).astype(float)
    extra = [
        ("GDP only", [L(d.gdp)], ["log2 GDP"]),
        ("Incidence rate only", [L(d.rate)], ["log2 rate"]),
        ("Articles per million only", [L(d.artpc)], ["log2 art/M"]),
        ("GDP + rate", [L(d.gdp), L(d.rate)], ["log2 GDP", "log2 rate"]),
        ("GDP + articles", [L(d.gdp), L(d.artpc)], ["log2 GDP", "log2 art/M"]),
        ("P5 (P3 + WHO region)", base_cols + [reg[c].to_numpy() for c in reg.columns], ["log2 GDP", "log2 rate", "log2 art/M"] + ["region " + c for c in reg.columns]),
        ("Income groups (ref low)", [inc[c].to_numpy() for c in inc.columns], ["income " + c for c in inc.columns]),
    ]
    out("\n=== Additional specifications (rate ratios; Poisson HC1 vs NB2 robust) ===")
    for name, cols, terms in extra:
        X = np.column_stack([np.ones(n)] + cols)
        b, mu = poisson_fit(X, y, off); k = X.shape[1]
        br_ = np.linalg.inv((X.T * mu) @ X); u = X * (y - mu)[:, None]
        cov = br_ @ (u.T @ u) @ br_ * n / (n - k)
        bn, a, mun, bread_n, cov_rob, ll = nb_fit(X, y, off)
        line = f"{name:28s} alpha={a:5.2f}"
        for j, t in enumerate(terms, start=1):
            rows.append(row(f"{name} Poisson HC1", t, b[j], np.sqrt(cov[j, j]), n))
            rows.append(row(f"{name} NB2 robust", t, bn[j], np.sqrt(cov_rob[j, j]), n))
            if not t.startswith("region"):
                r1, r2 = rows[-2], rows[-1]
                line += f" | {t}: Poisson {r1['rr']:.2f} ({r1['ci_lo']:.2f}-{r1['ci_hi']:.2f}); NB2 {r2['rr']:.2f} ({r2['ci_lo']:.2f}-{r2['ci_hi']:.2f})"
        out(line)
    out(f"\nSpearman correlations among predictors: GDP-rate {stats.spearmanr(d.gdp, d.rate)[0]:.2f}; GDP-articles/M {stats.spearmanr(d.gdp, d.artpc)[0]:.2f}; rate-articles/M {stats.spearmanr(d.rate, d.artpc)[0]:.2f}")
    pd.DataFrame(rows).to_csv(OUT / "negbin_models.csv", index=False)
    pd.DataFrame(turk).to_csv(OUT / "negbin_turkey.csv", index=False)
    (OUT / "negbin_summary.txt").write_text("\n".join(LINES), encoding="utf-8")


if __name__ == "__main__":
    main()
