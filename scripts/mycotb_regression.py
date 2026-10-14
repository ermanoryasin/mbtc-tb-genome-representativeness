#!/usr/bin/env python3
"""
MycoTB v5: "Kim halka acik depolara genom yukluyor?" -- Dunya Bankasi degiskenleri ile regresyon.

Gereksinim: mycotb_analysis.py ve mycotb_robustness.py ile AYNI klasorde olmali.
Internet gerekir (yalniz Dunya Bankasi API'si, birkac KB). Python standart kutuphanesi + curl.
Girdi: ~/Downloads/mycotb_pilot_out/ (onceki taramalar)
Cikti: ~/Downloads/mycotb_regression_out/ (+ mycotb_regression_results.zip)

Yapilanlar:
  1. Dunya Bankasi: ulke meta (gelir grubu), kisi basi GSYH (PPP), kisi basi saglik harcamasi,
     bilimsel makale sayisi (2010-2024 ortalamasi, mevcut yillar).
  2. Gelir grubuna gore: yuk payi vs genom payi.
  3. Poisson regresyon (offset = log kumulatif vaka, saglam/robust SE): genom sayisi ~ GSYH + bilimsel
     uretim + TB insidans hizi (+ bolge). Sonuc: "bir katina cikinca" oran orani (IRR).
  4. Lojistik regresyon: ulkede en az 1 halka acik genom olma olasiligi.
  5. Duyarlilik: ABD haric, yalniz dusuk/orta gelirli ulkeler, saglik harcamasi ile alternatif model.
  6. Beklenen vs gozlenen genom: Turkiye dahil en buyuk "acik" ulkeler.

Calistirma:  cd ~/Downloads && python3 mycotb_regression.py
"""
import csv
import datetime
import json
import math
import os
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
MBTC_BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import mycotb_analysis as A
    import mycotb_robustness as R
except ImportError:
    sys.exit("HATA: mycotb_analysis.py ve mycotb_robustness.py bu scriptle ayni klasorde olmali.")

OUT = Path(os.environ.get("MYCO_OUT3", str(MBTC_BASE / "mycotb_regression_out")))
WB = OUT / "worldbank"
P0, P1 = 2010, 2024
INDICATORS = {"gdp": "NY.GDP.PCAP.PP.KD", "health": "SH.XPD.CHEX.PC.CD", "articles": "IP.JRN.ARTC.SC"}
LINES = []


def out(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


# ------------------------------------------------------------------ indirme
def fetch(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = Path(str(dest) + ".part")
    cmd = ["curl", "-sS", "--fail", "--retry", "3", "--retry-delay", "5", "--max-time", "180",
           "-o", str(tmp), url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(f"curl cikis kodu {r.returncode}: {r.stderr.strip()[:300]}")
    tmp.replace(dest)


def load_wb():
    WB.mkdir(parents=True, exist_ok=True)
    fetch("https://api.worldbank.org/v2/country?format=json&per_page=400", WB / "countries.json")
    meta = {}
    with open(WB / "countries.json", encoding="utf-8") as f:
        data = json.load(f)
    for it in data[1]:
        if (it.get("region") or {}).get("id") in (None, "NA"):
            continue  # toplulastirmalar (World, Low income vb.) atla
        meta[it["id"]] = {"name": it.get("name", ""), "income": (it.get("incomeLevel") or {}).get("id", "")}
    ind = {}
    for key, code in INDICATORS.items():
        fetch(f"https://api.worldbank.org/v2/country/all/indicator/{code}?format=json&per_page=20000"
              f"&date={P0}:{P1}", WB / f"{code}.json")
        with open(WB / f"{code}.json", encoding="utf-8") as f:
            data = json.load(f)
        sums, cnt = defaultdict(float), defaultdict(int)
        for r in (data[1] or []):
            v = r.get("value")
            iso = r.get("countryiso3code")
            if v is not None and iso:
                sums[iso] += float(v)
                cnt[iso] += 1
        ind[key] = {k: sums[k] / cnt[k] for k in sums if cnt[k]}
    return meta, ind


def who_extra():
    iso, pop = {}, defaultdict(dict)
    with open(A.WHO, "r", encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            c = A.WHO_MAP.get(r["country"], r["country"])
            iso[c] = r.get("iso3", "")
            try:
                pop[c][int(r["year"])] = float(r["e_pop_num"])
            except (ValueError, KeyError, TypeError):
                pass
    return iso, pop


# ------------------------------------------------------------------ lineer cebir / GLM
def xtwx(X, w):
    p = len(X[0])
    M = [[0.0] * p for _ in range(p)]
    for row, wi in zip(X, w):
        for a in range(p):
            ra = row[a] * wi
            for b in range(a, p):
                M[a][b] += ra * row[b]
    for a in range(p):
        for b in range(a):
            M[a][b] = M[b][a]
    return M


def solve(M, b):
    n = len(M)
    A_ = [row[:] + [b[i]] for i, row in enumerate(M)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(A_[r][i]))
        if abs(A_[piv][i]) < 1e-12:
            raise ValueError("tekil matris (cok ilintili/degiskenlik yok)")
        A_[i], A_[piv] = A_[piv], A_[i]
        for r in range(i + 1, n):
            f = A_[r][i] / A_[i][i]
            for c in range(i, n + 1):
                A_[r][c] -= f * A_[i][c]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (A_[i][n] - sum(A_[i][j] * x[j] for j in range(i + 1, n))) / A_[i][i]
    return x


def invert(M):
    n = len(M)
    cols = []
    for j in range(n):
        e = [0.0] * n
        e[j] = 1.0
        cols.append(solve(M, e))
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def matmul(A_, B_):
    return [[sum(A_[i][k] * B_[k][j] for k in range(len(B_))) for j in range(len(B_[0]))] for i in range(len(A_))]


def glm(X, y, family, offset=None, iters=200):
    n, p = len(X), len(X[0])
    off = offset or [0.0] * n
    beta = [0.0] * p
    if family == "poisson":
        beta[0] = math.log(max(sum(y), 0.5) / sum(math.exp(o) for o in off))
    for _ in range(iters):
        eta = [max(-30.0, min(30.0, sum(b * x for b, x in zip(beta, row)) + o)) for row, o in zip(X, off)]
        if family == "poisson":
            mu = [math.exp(e) for e in eta]
            w = [max(m, 1e-9) for m in mu]
            z = [(e - o) + (yy - m) / max(m, 1e-9) for e, o, yy, m in zip(eta, off, y, mu)]
        else:
            pr = [1 / (1 + math.exp(-e)) for e in eta]
            w = [max(q * (1 - q), 1e-9) for q in pr]
            z = [e + (yy - q) / wi for e, yy, q, wi in zip(eta, y, pr, w)]
        M = xtwx(X, w)
        rhs = [sum(X[i][a] * w[i] * z[i] for i in range(n)) for a in range(p)]
        new = solve(M, rhs)
        delta = max(abs(a - b) for a, b in zip(new, beta))
        beta = new
        if delta < 1e-8:
            break
    eta = [max(-30.0, min(30.0, sum(b * x for b, x in zip(beta, row)) + o)) for row, o in zip(X, off)]
    if family == "poisson":
        mu = [math.exp(e) for e in eta]
        w = [max(m, 1e-9) for m in mu]
        res = [yy - m for yy, m in zip(y, mu)]
    else:
        mu = [1 / (1 + math.exp(-e)) for e in eta]
        w = [max(q * (1 - q), 1e-9) for q in mu]
        res = [yy - q for yy, q in zip(y, mu)]
    bread = invert(xtwx(X, w))
    meat = [[0.0] * p for _ in range(p)]
    for row, r_ in zip(X, res):
        for a in range(p):
            for b in range(p):
                meat[a][b] += r_ * r_ * row[a] * row[b]
    V = matmul(matmul(bread, meat), bread)
    c = n / max(n - p, 1)
    se = [math.sqrt(max(V[i][i] * c, 0.0)) for i in range(p)]
    return beta, se, mu


def pnorm2(z):
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


# ------------------------------------------------------------------ model kurma
def build_X(rows, terms):
    names = ["(intercept)"]
    cols = []
    for t in terms:
        if t == "log2gdp":
            names.append(t)
            cols.append([math.log2(r["gdp"]) for r in rows])
        elif t == "log2health":
            names.append(t)
            cols.append([math.log2(r["health"]) for r in rows])
        elif t == "log2art":
            names.append(t)
            cols.append([math.log2(r["articles"]) for r in rows])
        elif t == "log2rate":
            names.append(t)
            cols.append([math.log2(r["rate"]) for r in rows])
        elif t == "log2artpc":
            names.append(t)
            cols.append([math.log2(r["artpc"]) for r in rows])
        elif t == "log2cases":
            names.append(t)
            cols.append([math.log2(r["cases"]) for r in rows])
        elif t.startswith("inc_"):
            names.append(t)
            cols.append([1.0 if r["income"] == t[4:] else 0.0 for r in rows])
        elif t.startswith("reg_"):
            names.append(t)
            cols.append([1.0 if r["region"] == t[4:] else 0.0 for r in rows])
    X = [[1.0] + [cols[j][i] for j in range(len(cols))] for i in range(len(rows))]
    return X, names


def fit(rows, terms, family, label, writer):
    X, names = build_X(rows, terms)
    try:
        if family == "poisson":
            y = [float(r["obs"]) for r in rows]
            off = [math.log(r["cases"]) for r in rows]
            beta, se, mu = glm(X, y, "poisson", off)
        else:
            y = [1.0 if r["obs"] > 0 else 0.0 for r in rows]
            beta, se, mu = glm(X, y, "logit")
    except ValueError as e:
        out(f"[{label}] hesaplanamadi: {e}")
        return rows, None, None
    out(f"[{label}] n={len(rows)}")
    for nm, b, s in zip(names, beta, se):
        if nm == "(intercept)":
            continue
        z = b / s if s else float("nan")
        est = math.exp(b)
        lo_, hi_ = math.exp(b - 1.96 * s), math.exp(b + 1.96 * s)
        tag = "IRR" if family == "poisson" else "OR"
        out(f"   {nm}: {tag}={est:.2f} [{lo_:.2f}-{hi_:.2f}], p={pnorm2(z):.3g}")
        writer.writerow([label, family, nm, f"{b:.4f}", f"{est:.4f}", f"{lo_:.4f}", f"{hi_:.4f}",
                         f"{s:.4f}", f"{pnorm2(z):.4g}", len(rows)])
    return rows, beta, mu


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out(f"MycoTB regresyon analizi v5 - {datetime.datetime.now():%Y-%m-%d %H:%M}")
    try:
        meta, ind = load_wb()
    except Exception as e:  # noqa: BLE001
        sys.exit(f"HATA: Dunya Bankasi verisi indirilemedi: {e}\n"
                 f"Internet baglantinizi kontrol edip ayni komutu tekrar calistirin.")
    out(f"Dunya Bankasi: {len(meta)} ulke; GSYH {len(ind['gdp'])}, saglik {len(ind['health'])}, "
        f"makale {len(ind['articles'])} ulke icin veri")

    attr, samples, _ = A.load_and_clean()
    inc, _lo, _hi, region = R.load_who_full()
    iso, pop = who_extra()
    gen = R.gen_counts(samples, "c", P0, P1)
    cases = R.psum(inc, P0, P1)

    rows = []
    for c, cs in cases.items():
        i3 = iso.get(c, "")
        m = meta.get(i3)
        pop_sum = sum(v for y, v in pop.get(c, {}).items() if P0 <= y <= P1)
        rows.append({
            "country": c, "iso3": i3, "cases": cs, "obs": gen.get(c, 0), "region": region.get(c, ""),
            "income": (m or {}).get("income", ""), "gdp": ind["gdp"].get(i3), "health": ind["health"].get(i3),
            "articles": ind["articles"].get(i3), "rate": (1e5 * cs / pop_sum) if pop_sum else None,
            "artpc": (ind["articles"].get(i3) / (pop_sum / (P1 - P0 + 1) / 1e6))
            if (pop_sum and ind["articles"].get(i3)) else None,
        })
    with open(OUT / "analysis_dataset.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    out(f"WHO ulke sayisi: {len(rows)}; ulkesi Dunya Bankasi ile eslesen: {sum(1 for r in rows if r['income'])}")

    # ---------------- gelir grubu tablosu
    out("")
    out(f"=== 1. GELIR GRUBUNA GORE ({P0}-{P1}; Dunya Bankasi guncel siniflamasi) ===")
    grp = {}
    tot_c = sum(r["cases"] for r in rows if r["income"])
    tot_g = sum(r["obs"] for r in rows if r["income"])
    names = {"LIC": "Dusuk", "LMC": "Orta-alt", "UMC": "Orta-ust", "HIC": "Yuksek"}
    with open(OUT / "income_groups.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["income_group", "countries", "zero_genome_countries", "burden_share_pct",
                    "genome_share_pct", "genomes_per_1000_cases", "representation_index"])
        for g in ["LIC", "LMC", "UMC", "HIC"]:
            rs = [r for r in rows if r["income"] == g]
            if not rs:
                continue
            cg = sum(r["cases"] for r in rs)
            gg = sum(r["obs"] for r in rs)
            bs, gs = 100.0 * cg / tot_c, (100.0 * gg / tot_g if tot_g else 0)
            ri = gs / bs if bs else 0
            z = sum(1 for r in rs if r["obs"] == 0)
            out(f"  {names[g]} ({g}): {len(rs)} ulke, sifir genomlu {z} | yuk %{bs:.1f} | genom %{gs:.1f} | "
                f"genom/1000 vaka={1000.0 * gg / cg:.2f} | RI={ri:.2f}")
            w.writerow([g, len(rs), z, f"{bs:.2f}", f"{gs:.2f}", f"{1000.0 * gg / cg:.3f}", f"{ri:.3f}"])

    # ---------------- regresyonlar
    base = [r for r in rows if r["gdp"] and r["artpc"] and r["rate"] and r["income"] and r["region"]
            and r["gdp"] > 0 and r["artpc"] > 0 and r["rate"] > 0 and r["cases"] > 0]
    out("")
    out(f"=== 2. REGRESYON ANALIZLERI (analiz kumesi: {len(base)} ulke; eksik GSYH/makale/bolge dislandi) ===")
    out("Poisson: genom sayisi ~ ..., offset=log(kumulatif vaka), robust SE. IRR = prediktor bir KATINA "
        "cikinca genom/vaka oraninin carpani.")
    out("Bilimsel uretim = milyon nufus basina yillik makale (toplam makale degil; ulke buyuklugu ile karismasin).")
    regs = sorted({r["region"] for r in base})
    reg_terms = [f"reg_{x}" for x in regs[1:]]
    main_terms = ["log2gdp", "log2rate", "log2artpc"]
    with open(OUT / "models.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model", "family", "term", "coef", "ratio_exp_coef", "ci_lo", "ci_hi", "robust_se", "p", "n"])
        fit(base, ["inc_LMC", "inc_UMC", "inc_HIC"], "poisson", "P1 gelir grubu (ref=LIC)", w)
        fit(base, ["log2gdp", "log2rate"], "poisson", "P2 GSYH + insidans hizi", w)
        rows3, beta3, mu3 = fit(base, main_terms, "poisson", "P3 ANA MODEL: GSYH + insidans hizi + makale/milyon", w)
        rows4, beta4, mu4 = fit(base, main_terms + ["log2cases"], "poisson",
                                "P4 P3 + log2 kumulatif vaka (olcek)", w)
        fit(base, main_terms + reg_terms, "poisson", f"P5 P3 + WHO bolgesi (ref={regs[0]})", w)
        hb = [r for r in base if r["health"] and r["health"] > 0]
        fit(hb, ["log2health", "log2rate", "log2artpc"], "poisson", "P6 saglik harcamasi ile (GSYH yerine)", w)
        fit([r for r in base if r["country"] != "USA"], main_terms, "poisson", "S1 ABD haric", w)
        fit([r for r in base if r["income"] != "HIC"], main_terms, "poisson", "S2 yalniz dusuk/orta gelirli", w)
        out("")
        out("Lojistik: ulkede >=1 halka acik genom olma olasiligi (OR = bir katina cikinca odds carpani)")
        try:
            fit(base, main_terms + ["log2cases"], "logit", "L1 >=1 genom (vaka olcegi dahil)", w)
        except ValueError as e:
            out(f"L1 hesaplanamadi: {e}")

    # ---------------- beklenen vs gozlenen
    out("")
    out("=== 3. BEKLENEN vs GOZLENEN GENOM (P3 ana model; P4 olcek duyarliligi) ===")
    if mu4 is None:
        out("(P4 hesaplanamadi; P4 beklentisi yerine P3 kullanildi)")
        mu4 = mu3
    res = []
    for r, m, m4 in zip(rows3, mu3, mu4):
        res.append({"country": r["country"], "region": r["region"], "income": r["income"],
                    "cases": int(r["cases"]), "observed": r["obs"], "expected": m, "expected_p4": m4,
                    "obs_over_exp": (r["obs"] / m) if m else None, "gap": m - r["obs"], "gap4": m4 - r["obs"]})
    res.sort(key=lambda x: -x["gap"])
    with open(OUT / "expected_vs_observed.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "who_region", "income_group", "cumulative_cases", "observed_genomes",
                    "expected_genomes_P3", "expected_genomes_P4", "observed_over_expected_P3",
                    "expected_minus_observed_P3"])
        for x in res:
            w.writerow([x["country"], x["region"], x["income"], x["cases"], x["observed"],
                        f"{x['expected']:.1f}", f"{x['expected_p4']:.1f}",
                        f"{x['obs_over_exp']:.3f}" if x["obs_over_exp"] is not None else "",
                        f"{x['gap']:.1f}"])
    out("Beklentinin en cok altinda kalan 15 ulke (P3: beklenen - gozlenen):")
    for x in res[:15]:
        out(f"  {x['country']} ({x['income']}, {x['region']}): gozlenen={x['observed']}, "
            f"beklenen={x['expected']:.0f}, vaka={x['cases']}")
    res4 = sorted(res, key=lambda x: -x["gap4"])
    tr = [i for i, x in enumerate(res) if x["country"] == "Turkey"]
    tr4 = [i for i, x in enumerate(res4) if x["country"] == "Turkey"]
    if tr:
        x = res[tr[0]]
        out(f"Turkey: gozlenen={x['observed']}, beklenen P3={x['expected']:.0f} (sira {tr[0] + 1}/{len(res)}), "
            f"P4={x['expected_p4']:.0f} (sira {tr4[0] + 1}/{len(res)})")
    else:
        out("Turkey analiz kumesinde yok (GSYH/makale/bolge eksik olabilir).")
    over = sorted(res, key=lambda x: -(x["observed"] - x["expected"]))[:8]
    out("Beklentinin en cok ustunde 8 ulke: " + ", ".join(
        f"{x['country']} ({x['observed']} vs {x['expected']:.0f})" for x in over))

    out("")
    out("Yorum notu: modeller aciklayici/ilintisel, nedensel degil; genomlar bagimsiz izolat degil; "
        "bagimli degisken 'halka acik depoya yukleme', sekanslama kapasitesi degil.")
    text = "\n".join(LINES)
    (OUT / "summary.txt").write_text(text, encoding="utf-8")
    z = OUT / "mycotb_regression_results.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in ["summary.txt", "income_groups.csv", "models.csv", "expected_vs_observed.csv",
                  "analysis_dataset.csv"]:
            p = OUT / n
            if p.exists():
                zf.write(p, n)
    print(f"\nBITTI. Bu dosyayi Claude'a yukleyin: {z}")


if __name__ == "__main__":
    main()
