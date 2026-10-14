#!/usr/bin/env python3
"""
MycoTB analiz scripti (v3, OFFLINE): temizleme + temsil/kapsama metrikleri.
(representativeness / coverage of publicly available M. tuberculosis genomes vs TB burden)

Girdi (onceki taramadan, indirme YOK):
  ~/Downloads/mycotb_pilot_out/mtb_read_runs_all.tsv
  ~/Downloads/mycotb_pilot_out/who_tb_estimates.csv
Cikti: ~/Downloads/mycotb_analysis_out/  (+ mycotb_analysis_results.zip)

Neler yapar:
  1. Temizleme hatti ve her adimda kayip tablosu (attrition):
     siki etiket -> laboratuvar suslari disi -> WGS/genomik -> sample duzeyine indirgeme
  2. Ulke bazli: koleksiyon yili ile eslestirilmis genom / 1.000 tahmini vaka (2010-2024)
  3. Donemsel egilim (2010-14, 2015-19, 2020-24): kapsama, Gini, yuk payi
  4. Eksik ulke bilgisi icin duyarlilik senaryolari
  5. Proje (study) yogunlasmasi, metadata doluluğu, yukleme gecikmesi
  6. Turkiye ve komsulari, "genomik kor nokta" listesi

Calistirma:  cd ~/Downloads && python3 mycotb_analysis.py
"""
import csv
import datetime
import os
import re
import statistics
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
MBTC_BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))

csv.field_size_limit(10 ** 8)

HOME = Path.home()
IN_DIR = Path(os.environ.get("MYCO_IN", str(MBTC_BASE / "mycotb_pilot_out")))
OUT = Path(os.environ.get("MYCO_OUT", str(MBTC_BASE / "mycotb_analysis_out")))
SRC = IN_DIR / "mtb_read_runs_all.tsv"
WHO = IN_DIR / "who_tb_estimates.csv"

NOW_YEAR = datetime.date.today().year
P0, P1 = 2010, 2024           # ana donem (WHO verisi 2024'e kadar)
PERIODS = [(2010, 2014), (2015, 2019), (2020, 2024)]
THR = 10.0                    # esik: 1.000 yillik vakaya <10 genom = dusuk kapsama

MISSING = {
    "", "missing", "not applicable", "not collected", "not provided", "na", "n/a",
    "unknown", "none", "restricted access", "not available", "null", "-", "not determined",
    "not recorded", "unspecified", "not known",
}
COUNTRY_FIX = {
    "turkey": "Turkey", "turkiye": "Turkey", "türkiye": "Turkey", "türkei": "Turkey",
    "republic of turkey": "Turkey",
    "usa": "USA", "united states": "USA", "united states of america": "USA", "us": "USA",
    "uk": "United Kingdom", "great britain": "United Kingdom", "england": "United Kingdom",
    "viet nam": "Vietnam", "russian federation": "Russia",
    "republic of korea": "South Korea", "korea": "South Korea",
    "people's republic of china": "China", "czech republic": "Czechia", "burma": "Myanmar",
    "iran, islamic republic of": "Iran", "islamic republic of iran": "Iran",
    "syrian arab republic": "Syria", "republic of georgia": "Georgia", "swaziland": "Eswatini",
}
WHO_MAP = {
    "Türkiye": "Turkey", "Russian Federation": "Russia", "Iran (Islamic Republic of)": "Iran",
    "Viet Nam": "Vietnam", "United States of America": "USA",
    "United Kingdom of Great Britain and Northern Ireland": "United Kingdom",
    "Republic of Korea": "South Korea", "Republic of Moldova": "Moldova",
    "United Republic of Tanzania": "Tanzania", "Côte d'Ivoire": "Cote d'Ivoire",
    "Bolivia (Plurinational State of)": "Bolivia", "Venezuela (Bolivarian Republic of)": "Venezuela",
    "Lao People's Democratic Republic": "Laos", "Syrian Arab Republic": "Syria",
    "Democratic People's Republic of Korea": "North Korea",
    "Netherlands (Kingdom of the)": "Netherlands",
    "China, Hong Kong SAR": "Hong Kong", "China, Macao SAR": "Macao",
}
NON_STRICT_TOKENS = ("variant", "complex", "bovis", "caprae", "africanum", "microti", "canettii",
                     "orygis", "pinnipedii", "mungi", "suricattae", "bcg")
LAB_RE = re.compile(r"(h37rv|h37ra|erdman|hn878|cdc1551|atcc)", re.IGNORECASE)
NEIGHBORS = ["Turkey", "Georgia", "Armenia", "Azerbaijan", "Iran", "Iraq", "Syria", "Greece",
             "Bulgaria", "Romania", "Russia", "Ukraine", "Moldova", "Lebanon", "Egypt"]

LINES = []


def out(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


def present(v):
    return (v or "").strip().lower() not in MISSING


def norm_country(raw):
    c = (raw or "").strip()
    if c.lower() in MISSING:
        return "(missing)"
    c = c.split(":")[0].strip()
    if c.lower() in MISSING:
        return "(missing)"
    if ";" in c or "|" in c:
        return "(multiple)"
    if c.lower() in COUNTRY_FIX:
        return COUNTRY_FIX[c.lower()]
    if c.isupper() and len(c) > 3:
        c = c.title()
    return COUNTRY_FIX.get(c.lower(), c)


def is_strict_mtb(name):
    n = (name or "").lower().strip()
    return n.startswith("mycobacterium tuberculosis") and not any(t in n for t in NON_STRICT_TOKENS)


def year_int(v, lo=1900, hi=None):
    s = (v or "").strip()
    if len(s) >= 4 and s[:4].isdigit():
        y = int(s[:4])
        if lo <= y <= (hi or NOW_YEAR):
            return y
    return None


def pct(a, b):
    return f"{(100.0 * a / b):.1f}" if b else "0.0"


def is_missing_country(c):
    return c in ("(missing)", "(multiple)")


# ------------------------------------------------------------ 1. temizleme
def load_and_clean():
    if not SRC.exists():
        sys.exit(f"HATA: {SRC} bulunamadi. Once mycotb_pilot.py'yi calistirin.")
    attr = []
    runs_all = 0
    all_samples = set()
    cur = []
    with open(SRC, "r", encoding="utf-8", errors="replace", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            runs_all += 1
            all_samples.add(r.get("sample_accession") or r.get("run_accession"))
            cur.append({
                "run": r.get("run_accession", ""),
                "sample": r.get("sample_accession") or r.get("run_accession", ""),
                "study": r.get("study_accession", ""),
                "name": r.get("scientific_name", ""),
                "strain": r.get("strain", ""),
                "country": r.get("country", ""),
                "cdate": r.get("collection_date", ""),
                "pub": r.get("first_public", ""),
                "host": r.get("host", ""),
                "iso": r.get("isolation_source", ""),
                "strat": (r.get("library_strategy") or "").upper(),
                "lsrc": (r.get("library_source") or "").upper(),
                "center": r.get("center_name", ""),
                "title": r.get("study_title", ""),
            })
    attr.append(("0. Tum run (taxid 1773 agaci)", len(cur), len({x["sample"] for x in cur})))

    def step(label, rows, pred):
        rows = [x for x in rows if pred(x)]
        attr.append((label, len(rows), len({x["sample"] for x in rows})))
        return rows

    cur = step("1. Siki etiket (bovis/africanum/BCG/complex vb. haric)", cur,
               lambda x: is_strict_mtb(x["name"]))
    cur = step("2. Laboratuvar referans suslari haric (H37Rv/Erdman/HN878/CDC1551/ATCC)", cur,
               lambda x: not (LAB_RE.search(x["name"] or "") or LAB_RE.search(x["strain"] or "")))
    cur = step("3. WGS + genomik kaynak (RNA-Seq/Tn-Seq/amplicon haric)", cur,
               lambda x: x["strat"] == "WGS" and x["lsrc"] in ("GENOMIC", ""))

    # sample duzeyine indirgeme
    by = defaultdict(list)
    for x in cur:
        by[x["sample"]].append(x)
    samples = []
    conflicts = 0
    for sa, xs in by.items():
        cs = [norm_country(x["country"]) for x in xs]
        good = [c for c in cs if not is_missing_country(c)]
        if len(set(good)) > 1:
            conflicts += 1
        if good:
            c = good[0]
        else:
            c = "(multiple)" if "(multiple)" in cs else "(missing)"
        cy = next((y for y in (year_int(x["cdate"]) for x in xs) if y), None)
        pubs = [year_int(x["pub"]) for x in xs if year_int(x["pub"])]
        samples.append({
            "sample": sa, "country": c, "cyear": cy,
            "pyear": min(pubs) if pubs else None,
            "study": xs[0]["study"], "center": xs[0]["center"], "title": xs[0]["title"],
            "has_date": cy is not None,
            "has_host": any(present(x["host"]) for x in xs),
            "has_iso": any(present(x["iso"]) for x in xs),
        })
    attr.append(("4. Sample (BioSample) duzeyine indirgeme = analiz birimi", len(cur), len(samples)))
    return attr, samples, conflicts


# ------------------------------------------------------------ WHO
def load_who():
    if not WHO.exists():
        return None
    inc = defaultdict(dict)
    with open(WHO, "r", encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            try:
                y = int(r["year"])
                v = float(r["e_inc_num"])
            except (ValueError, KeyError, TypeError):
                continue
            inc[WHO_MAP.get(r["country"], r["country"])][y] = v
    return inc


def period_inc(who, y0, y1):
    d = {}
    for c, yy in who.items():
        s = sum(v for y, v in yy.items() if y0 <= y <= y1)
        if s > 0:
            d[c] = s
    return d


def metrics(gen, inc):
    cs = [c for c, i in inc.items() if i > 0]
    ti = sum(inc[c] for c in cs)
    tw = sum(gen.get(c, 0) for c in cs)
    if ti == 0 or tw == 0:
        return None
    lo = sum(inc[c] for c in cs if 1000.0 * gen.get(c, 0) / inc[c] < THR)
    order = sorted(cs, key=lambda c: gen.get(c, 0) / inc[c])
    cw = area = 0.0
    for c in order:
        pi, pw = inc[c] / ti, gen.get(c, 0) / tw
        area += pi * (cw + pw / 2)
        cw += pw
    top = sorted(cs, key=lambda c: -inc[c])[:30]
    return {
        "countries": len(cs), "per1000": 1000.0 * tw / ti, "burden_low_pct": 100.0 * lo / ti,
        "gini": 1 - 2 * area, "top30_burden_pct": 100.0 * sum(inc[c] for c in top) / ti,
        "top30_genome_pct": 100.0 * sum(gen.get(c, 0) for c in top) / tw,
        "zero_countries": sum(1 for c in cs if gen.get(c, 0) == 0),
        "genomes": tw,
    }


def fmt_metrics(m):
    return (f"ulke={m['countries']}, genom/1000 vaka={m['per1000']:.1f}, "
            f"yukun %{m['burden_low_pct']:.1f}'i <{THR:.0f}/1000 ulkelerde, Gini={m['gini']:.3f}, "
            f"ilk30 yuk %{m['top30_burden_pct']:.1f} / genom %{m['top30_genome_pct']:.1f}, "
            f"sifir genomlu ulke={m['zero_countries']}")


# ------------------------------------------------------------ ana analiz
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out(f"MycoTB offline analiz v3 - {datetime.datetime.now():%Y-%m-%d %H:%M}")
    attr, samples, conflicts = load_and_clean()
    out("")
    out("=== 1. TEMIZLEME HATTI (kayip tablosu) ===")
    with open(OUT / "attrition.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["step", "runs", "samples"])
        for a in attr:
            w.writerow(a)
            out(f"{a[0]}: run={a[1]}, sample={a[2]}")
    out(f"Ayni sample icinde celisen ulke etiketi: {conflicts}")
    n = len(samples)
    miss = [s for s in samples if is_missing_country(s["country"])]
    out(f"Analiz birimi: {n} WGS sample | ulke eksik/cok-ulkeli: {len(miss)} (%{pct(len(miss), n)})")
    valid_date = sum(1 for s in samples if s["has_date"])
    out(f"Gecerli koleksiyon yili olan sample: {valid_date} (%{pct(valid_date, n)})")

    # ulke bazli yapilar
    by_c = defaultdict(list)
    for s in samples:
        by_c[s["country"]].append(s)
    gen_all = {c: len(v) for c, v in by_c.items() if not is_missing_country(c)}
    gen_per = {}
    for p in PERIODS + [(P0, P1)]:
        g = Counter()
        for s in samples:
            if not is_missing_country(s["country"]) and s["cyear"] and p[0] <= s["cyear"] <= p[1]:
                g[s["country"]] += 1
        gen_per[p] = g

    who = load_who()
    if who is None:
        out("UYARI: WHO tablosu bulunamadi; adalet/kapsama metrikleri atlandi.")
    main_inc = period_inc(who, P0, P1) if who else {}
    inc_2024 = {c: yy.get(2024, 0) for c, yy in who.items()} if who else {}
    in_period = sum(gen_per[(P0, P1)].values())
    out(f"Ulkesi bilinen ve koleksiyon yili {P0}-{P1} olan sample: {in_period}")

    # --- 2. ulke tablosu
    studies_c = {}
    for c, v in by_c.items():
        cnt = Counter(s["study"] for s in v)
        top1 = cnt.most_common(1)[0][1] if cnt else 0
        studies_c[c] = (len(cnt), top1 / len(v) if v else 0.0, cnt)
    rows = []
    for c, v in sorted(by_c.items(), key=lambda kv: -len(kv[1])):
        k = len(v)
        d = sum(1 for s in v if s["has_date"])
        h = sum(1 for s in v if s["has_host"])
        i = sum(1 for s in v if s["has_iso"])
        lags = [s["pyear"] - s["cyear"] for s in v if s["cyear"] and s["pyear"] and -1 <= s["pyear"] - s["cyear"] <= 60]
        inc_p = main_inc.get(c, 0)
        gp = gen_per[(P0, P1)].get(c, 0)
        rows.append({
            "country": c, "wgs_samples_clean": k, f"wgs_samples_{P0}_{P1}": gp,
            f"who_incidence_{P0}_{P1}": int(inc_p), "who_incidence_2024": int(inc_2024.get(c, 0)),
            "genomes_per_1000_cases_period": f"{1000.0 * gp / inc_p:.2f}" if inc_p else "",
            "n_studies": studies_c[c][0], "top_study_share_pct": f"{100 * studies_c[c][1]:.1f}",
            "pct_collection_date": pct(d, k), "pct_host": pct(h, k), "pct_isolation_source": pct(i, k),
            "median_publication_lag_years": f"{statistics.median(lags):.1f}" if len(lags) >= 10 else "",
        })
    with open(OUT / "country_table.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    unmatched = [(len(v), c) for c, v in by_c.items()
                 if not is_missing_country(c) and who and c not in who]
    unmatched.sort(reverse=True)
    if unmatched:
        out("WHO'da eslesmeyen ulke etiketleri (sample): " + ", ".join(f"{c}={n_}" for n_, c in unmatched[:8]))

    # --- 3. ana kapsama metrikleri + donemsel egilim
    if who:
        out("")
        out(f"=== 2. TEMSIL / KAPSAMA METRIKLERI (representativeness / coverage), koleksiyon yili {P0}-{P1} ===")
        m = metrics(gen_per[(P0, P1)], main_inc)
        out(fmt_metrics(m))
        out("")
        out("=== 3. DONEMSEL EGILIM (kapsama zamanla duzeliyor mu?) ===")
        with open(OUT / "periods.csv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["period", "countries", "genomes", "genomes_per_1000_cases", "pct_burden_in_low_coverage",
                        "gini", "top30_burden_pct", "top30_genome_pct", "zero_genome_countries"])
            for p in PERIODS:
                mp = metrics(gen_per[p], period_inc(who, *p))
                if mp:
                    out(f"{p[0]}-{p[1]}: {fmt_metrics(mp)}")
                    w.writerow([f"{p[0]}-{p[1]}", mp["countries"], mp["genomes"], f"{mp['per1000']:.2f}",
                                f"{mp['burden_low_pct']:.1f}", f"{mp['gini']:.3f}", f"{mp['top30_burden_pct']:.1f}",
                                f"{mp['top30_genome_pct']:.1f}", mp["zero_countries"]])

        # --- 4. duyarlilik
        out("")
        out("=== 4. EKSIK ULKE BILGISI DUYARLILIK SENARYOLARI ===")
        base = dict(gen_per[(P0, P1)])
        nm = sum(1 for s in miss if s["cyear"] and P0 <= s["cyear"] <= P1)
        out(f"Donemde koleksiyon yili gecerli ama ulkesi eksik sample: {nm}")
        ti = sum(main_inc.values())
        tw = sum(base.get(c, 0) for c in main_inc)
        scen = []
        scen.append(("S0 gozlenen (eksikler dislanir)", dict(base)))
        k = (tw + nm) / tw if tw else 1
        scen.append(("S1 eksikler gozlenen dagilima orantili", {c: base.get(c, 0) * k for c in main_inc}))
        scen.append(("S2 eksikler TB yukune orantili (uc, esitligi lehine)",
                     {c: base.get(c, 0) + nm * main_inc[c] / ti for c in main_inc}))
        low = [c for c in main_inc if 1000.0 * base.get(c, 0) / main_inc[c] < THR]
        li = sum(main_inc[c] for c in low)
        g3 = {c: base.get(c, 0) for c in main_inc}
        for c in low:
            g3[c] += nm * main_inc[c] / li
        scen.append(("S3 eksiklerin tamami dusuk kapsamali ulkelere (en uc)", g3))
        with open(OUT / "scenarios.csv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["scenario", "genomes_per_1000", "pct_burden_in_low_coverage", "gini", "zero_genome_countries"])
            for name, g in scen:
                ms = metrics(g, main_inc)
                out(f"{name}: genom/1000={ms['per1000']:.1f}, yukun %{ms['burden_low_pct']:.1f}'i dusuk kapsama, "
                    f"Gini={ms['gini']:.3f}, sifir genomlu ulke={ms['zero_countries']}")
                w.writerow([name, f"{ms['per1000']:.2f}", f"{ms['burden_low_pct']:.1f}",
                            f"{ms['gini']:.3f}", ms["zero_countries"]])

    # --- 5. proje yogunlasmasi
    out("")
    out("=== 5. PROJE (STUDY) YOGUNLASMASI ===")
    st = Counter(s["study"] for s in samples)
    top10 = sum(v for _, v in st.most_common(10))
    out(f"Benzersiz study: {len(st)} | ilk 10 study toplam sample'in %{pct(top10, n)}'ini olusturuyor")
    st_info = {}
    for s in samples:
        e = st_info.setdefault(s["study"], {"countries": Counter(), "center": s["center"], "title": s["title"]})
        e["countries"][s["country"]] += 1
    with open(OUT / "studies_top.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["study_accession", "samples", "top_countries", "center_name", "study_title"])
        for acc, v in st.most_common(30):
            e = st_info[acc]
            w.writerow([acc, v, "; ".join(f"{c}={k}" for c, k in e["countries"].most_common(3)),
                        e["center"][:80], e["title"][:160]])
    big = [c for c, v in sorted(by_c.items(), key=lambda kv: -len(kv[1])) if not is_missing_country(c)][:12]
    out("En cok sample'a sahip ulkelerde proje yogunlasmasi (study sayisi / en buyuk study payi):")
    for c in big:
        out(f"  {c}: {len(by_c[c])} sample | {studies_c[c][0]} study | en buyuk study %{100 * studies_c[c][1]:.1f}")

    # --- 6. metadata / yil
    out("")
    out("=== 6. METADATA DOLULUGU VE ULKE EKSIKLIGI (yayin yilina gore) ===")
    yy = defaultdict(lambda: [0, 0, 0, 0, 0])
    for s in samples:
        if not s["pyear"]:
            continue
        r = yy[s["pyear"]]
        r[0] += 1
        r[1] += is_missing_country(s["country"])
        r[2] += s["has_date"]
        r[3] += s["has_host"]
        r[4] += s["has_iso"]
    with open(OUT / "year_metadata.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["first_public_year", "samples", "pct_missing_country", "pct_valid_collection_date",
                    "pct_host", "pct_isolation_source"])
        for y in sorted(yy):
            r = yy[y]
            w.writerow([y, r[0], pct(r[1], r[0]), pct(r[2], r[0]), pct(r[3], r[0]), pct(r[4], r[0])])
        for y in sorted(yy)[-6:]:
            r = yy[y]
            out(f"  {y}: {r[0]} sample | ulke eksik %{pct(r[1], r[0])} | tarih %{pct(r[2], r[0])} | "
                f"konak %{pct(r[3], r[0])} | kaynak %{pct(r[4], r[0])}")
    lags = [s["pyear"] - s["cyear"] for s in samples if s["cyear"] and s["pyear"] and -1 <= s["pyear"] - s["cyear"] <= 60]
    if lags:
        out(f"Koleksiyon -> yayin gecikmesi (yil): medyan={statistics.median(lags):.1f}, n={len(lags)}")

    # --- 7. Turkiye, komsular, kor noktalar
    out("")
    out("=== 7. TURKIYE VE KOMSULARI ===")
    tbl = {r["country"]: r for r in rows}
    out("(sample / genom 2010-24 / WHO vaka 2010-24 / genom per 1000 vaka / study / en buyuk study %)")
    with open(OUT / "turkey_neighbors.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "wgs_samples_clean", f"wgs_samples_{P0}_{P1}", f"who_incidence_{P0}_{P1}",
                    "genomes_per_1000_cases", "n_studies", "top_study_share_pct"])
        for c in NEIGHBORS:
            r = tbl.get(c)
            if r:
                out(f"  {c}: {r['wgs_samples_clean']} / {r[f'wgs_samples_{P0}_{P1}']} / {r[f'who_incidence_{P0}_{P1}']} / "
                    f"{r['genomes_per_1000_cases_period'] or '-'} / {r['n_studies']} / %{r['top_study_share_pct']}")
                w.writerow([c, r["wgs_samples_clean"], r[f"wgs_samples_{P0}_{P1}"], r[f"who_incidence_{P0}_{P1}"],
                            r["genomes_per_1000_cases_period"], r["n_studies"], r["top_study_share_pct"]])
            elif who and c in who:
                inc_p = int(main_inc.get(c, 0))
                out(f"  {c}: 0 / 0 / {inc_p} / 0.00 / 0 / -")
                w.writerow([c, 0, 0, inc_p, "0.00", 0, ""])
    tr = by_c.get("Turkey", [])
    out(f"Turkey sample listesi: {len(tr)} | koleksiyon yillari: "
        + ", ".join(str(s['cyear'] or 'yok') for s in tr[:20]))

    if who:
        out("")
        out("=== 8. GENOMIK KOR NOKTALAR (WHO 2024 insidansi >= 5.000 ve genom < 1 / 1.000 vaka, 2010-24) ===")
        bs = []
        for c, i24 in inc_2024.items():
            ip = main_inc.get(c, 0)
            if i24 >= 5000 and ip:
                per = 1000.0 * gen_per[(P0, P1)].get(c, 0) / ip
                if per < 1.0:
                    bs.append((i24, c, gen_per[(P0, P1)].get(c, 0), per))
        bs.sort(reverse=True)
        with open(OUT / "blind_spots.csv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["country", "who_incidence_2024", f"genomes_{P0}_{P1}", "genomes_per_1000_cases"])
            for i24, c, g, per in bs:
                w.writerow([c, int(i24), g, f"{per:.3f}"])
        out(f"Kor nokta ulke sayisi: {len(bs)}; ilk 15:")
        for i24, c, g, per in bs[:15]:
            out(f"  {c}: vaka(2024)={i24:.0f}, genom={g}, genom/1000 vaka={per:.2f}")

    text = "\n".join(LINES)
    (OUT / "summary.txt").write_text(text, encoding="utf-8")
    z = OUT / "mycotb_analysis_results.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for name in ["summary.txt", "attrition.csv", "country_table.csv", "periods.csv", "scenarios.csv",
                     "studies_top.csv", "year_metadata.csv", "turkey_neighbors.csv", "blind_spots.csv"]:
            p = OUT / name
            if p.exists():
                zf.write(p, name)
    print(f"\nBITTI. Bu dosyayi Claude'a yukleyin: {z}")


if __name__ == "__main__":
    main()
