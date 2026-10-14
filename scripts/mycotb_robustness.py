#!/usr/bin/env python3
"""
MycoTB dayaniklilik analizi (v4, OFFLINE): "kapsama artiyor ama dagilim bozuluyor" bulgusu saglam mi?

Gereksinim: mycotb_analysis.py ile AYNI klasorde olmali (temizleme hattini ondan alir).
Girdi: ~/Downloads/mycotb_pilot_out/mtb_read_runs_all.tsv ve who_tb_estimates.csv (indirme YOK)
Cikti: ~/Downloads/mycotb_robustness_out/  (+ mycotb_robustness_results.zip)

Testler:
  A. Yayin-yili kohortlari (release-year): yayin gecikmesinden bagimsiz egilim
  B. Gecikme-eslestirilmis (lag<=3 yil) koleksiyon-yili egilimi
  C. Duyarlilik: ABD haric / her ulkenin en buyuk projesi haric / ikisi birden
  D. WHO bolgeleri: yuk payi vs genom payi (representation index)
  E. Belirsizlik: WHO alt/ust insidans tahminleriyle Gini ve kapsama araliklari (Monte Carlo)
  F. Kapsama kategorileri, Lorenz egrisi verisi, ulke bazli temsil indeksi
  G. Turkiye: sifir genomlu ulkeler arasindaki sirasi

Calistirma:  cd ~/Downloads && python3 mycotb_robustness.py
"""
import csv
import datetime
import os
import random
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
MBTC_BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    import mycotb_analysis as A
except ImportError:
    sys.exit("HATA: mycotb_analysis.py bu scriptle ayni klasorde olmali.")

OUT = Path(os.environ.get("MYCO_OUT2", str(MBTC_BASE / "mycotb_robustness_out")))
LINES = []

REL_WINDOWS = [(2013, 2015), (2016, 2018), (2019, 2021), (2022, 2024)]
LAG_PERIODS = [(2010, 2014), (2015, 2019), (2020, 2022)]
MAIN = (2010, 2024)
MC_N = 1000


def out(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


def load_who_full():
    inc, lo, hi, region = defaultdict(dict), defaultdict(dict), defaultdict(dict), {}
    with open(A.WHO, "r", encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            try:
                y = int(r["year"])
                v = float(r["e_inc_num"])
            except (ValueError, KeyError, TypeError):
                continue
            c = A.WHO_MAP.get(r["country"], r["country"])
            inc[c][y] = v
            try:
                lo[c][y] = float(r["e_inc_num_lo"])
                hi[c][y] = float(r["e_inc_num_hi"])
            except (ValueError, KeyError, TypeError):
                lo[c][y], hi[c][y] = v, v
            region[c] = r.get("g_whoregion", "") or ""
    return inc, lo, hi, region


def psum(d, y0, y1):
    res = {}
    for c, yy in d.items():
        s = sum(v for y, v in yy.items() if y0 <= y <= y1)
        if s > 0:
            res[c] = s
    return res


def gen_counts(samples, by, y0, y1, filt=None):
    g = Counter()
    for s in samples:
        if A.is_missing_country(s["country"]):
            continue
        y = s["cyear"] if by == "c" else s["pyear"]
        if not y or y < y0 or y > y1:
            continue
        if filt and not filt(s):
            continue
        g[s["country"]] += 1
    return g


def fm(m):
    if not m:
        return "hesaplanamadi"
    return (f"ulke={m['countries']}, genom/1000 vaka={m['per1000']:.2f}, Gini={m['gini']:.3f}, "
            f"ilk30 yuk %{m['top30_burden_pct']:.1f} / genom %{m['top30_genome_pct']:.1f}, "
            f"sifir genomlu ulke={m['zero_countries']}")


def pctile(vals, q):
    vals = sorted(vals)
    return vals[min(len(vals) - 1, max(0, int(round(q * (len(vals) - 1)))))]


def category(per):
    if per == 0:
        return "0"
    if per < 1:
        return "<1"
    if per < 10:
        return "1-10"
    return ">=10"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out(f"MycoTB dayaniklilik analizi v4 - {datetime.datetime.now():%Y-%m-%d %H:%M}")
    attr, samples, _ = A.load_and_clean()
    out(f"Analiz birimi: {len(samples)} temiz WGS sample")
    if not A.WHO.exists():
        sys.exit("HATA: WHO tablosu yok.")
    inc, lo, hi, region = load_who_full()

    # en buyuk proje (study) ulke basina
    by_c_study = defaultdict(Counter)
    for s in samples:
        if not A.is_missing_country(s["country"]):
            by_c_study[s["country"]][s["study"]] += 1
    dominant = {c: cnt.most_common(1)[0][0] for c, cnt in by_c_study.items() if cnt}

    # ------------------------------------------------ A. yayin-yili kohortlari
    out("")
    out("=== A. YAYIN-YILI KOHORTLARI (genom hangi yil yayinlandi; gecikme etkisinden bagimsiz) ===")
    with open(OUT / "release_cohorts.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["release_window", "countries", "genomes", "genomes_per_1000", "gini",
                    "top30_burden_pct", "top30_genome_pct", "zero_genome_countries"])
        for p in REL_WINDOWS:
            m = A.metrics(gen_counts(samples, "p", *p), psum(inc, *p))
            out(f"{p[0]}-{p[1]}: {fm(m)}")
            if m:
                w.writerow([f"{p[0]}-{p[1]}", m["countries"], m["genomes"], f"{m['per1000']:.3f}",
                            f"{m['gini']:.3f}", f"{m['top30_burden_pct']:.1f}",
                            f"{m['top30_genome_pct']:.1f}", m["zero_countries"]])

    # ------------------------------------------------ B. gecikme-eslestirilmis
    out("")
    out("=== B. GECIKME-ESLESTIRILMIS KOLEKSIYON-YILI EGILIMI (yalniz yayin gecikmesi <=3 yil olanlar) ===")
    lagf = lambda s: s["pyear"] is not None and s["cyear"] is not None and 0 <= s["pyear"] - s["cyear"] <= 3  # noqa: E731
    with open(OUT / "lag_matched.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["collection_period", "countries", "genomes", "genomes_per_1000", "gini",
                    "top30_burden_pct", "top30_genome_pct", "zero_genome_countries"])
        for p in LAG_PERIODS:
            m = A.metrics(gen_counts(samples, "c", *p, filt=lagf), psum(inc, *p))
            out(f"{p[0]}-{p[1]}: {fm(m)}")
            if m:
                w.writerow([f"{p[0]}-{p[1]}", m["countries"], m["genomes"], f"{m['per1000']:.3f}",
                            f"{m['gini']:.3f}", f"{m['top30_burden_pct']:.1f}",
                            f"{m['top30_genome_pct']:.1f}", m["zero_countries"]])

    # ------------------------------------------------ C. duyarlilik
    out("")
    out("=== C. DUYARLILIK ANALIZLERI (koleksiyon yili donemleri) ===")
    scen = [
        ("temel", None),
        ("ABD haric", lambda s: s["country"] != "USA"),
        ("her ulkenin en buyuk projesi haric", lambda s: s["study"] != dominant.get(s["country"])),
        ("ABD ve en buyuk projeler haric",
         lambda s: s["country"] != "USA" and s["study"] != dominant.get(s["country"])),
    ]
    with open(OUT / "sensitivity.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "period", "genomes", "genomes_per_1000", "gini", "top30_genome_pct",
                    "zero_genome_countries"])
        for name, flt in scen:
            out(f"[{name}]")
            for p in [(2010, 2014), (2015, 2019), (2020, 2024)]:
                m = A.metrics(gen_counts(samples, "c", *p, filt=flt), psum(inc, *p))
                out(f"  {p[0]}-{p[1]}: {fm(m)}")
                if m:
                    w.writerow([name, f"{p[0]}-{p[1]}", m["genomes"], f"{m['per1000']:.3f}",
                                f"{m['gini']:.3f}", f"{m['top30_genome_pct']:.1f}", m["zero_countries"]])

    # ------------------------------------------------ E. belirsizlik (Monte Carlo)
    out("")
    out(f"=== E. WHO BELIRSIZLIGI ILE GUVEN ARALIKLARI ({MC_N} Monte Carlo cekilisi, %95) ===")
    random.seed(1)
    with open(OUT / "uncertainty.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["period", "gini_point", "gini_lo", "gini_hi", "per1000_point", "per1000_lo", "per1000_hi"])
        for p in [(2010, 2014), (2015, 2019), (2020, 2024)]:
            gen = gen_counts(samples, "c", *p)
            pt, l_, h_ = psum(inc, *p), psum(lo, *p), psum(hi, *p)
            base = A.metrics(gen, pt)
            gs, ps = [], []
            for _ in range(MC_N):
                draw = {}
                for c, v in pt.items():
                    a, b = min(l_.get(c, v), v), max(h_.get(c, v), v)
                    draw[c] = max(random.triangular(a, b, v), 1e-9)
                m = A.metrics(gen, draw)
                if m:
                    gs.append(m["gini"])
                    ps.append(m["per1000"])
            if base and gs:
                out(f"{p[0]}-{p[1]}: Gini={base['gini']:.3f} [{pctile(gs, .025):.3f}-{pctile(gs, .975):.3f}], "
                    f"genom/1000={base['per1000']:.2f} [{pctile(ps, .025):.2f}-{pctile(ps, .975):.2f}]")
                w.writerow([f"{p[0]}-{p[1]}", f"{base['gini']:.4f}", f"{pctile(gs, .025):.4f}",
                            f"{pctile(gs, .975):.4f}", f"{base['per1000']:.3f}", f"{pctile(ps, .025):.3f}",
                            f"{pctile(ps, .975):.3f}"])

    # ------------------------------------------------ D. WHO bolgeleri
    out("")
    out("=== D. WHO BOLGELERI: yuk payi vs genom payi (RI = genom payi / yuk payi) ===")
    reg_rows = []
    for p in [(2010, 2014), (2015, 2019), (2020, 2024)]:
        gen = gen_counts(samples, "c", *p)
        pt = psum(inc, *p)
        tot_i = sum(pt.values())
        tot_g = sum(gen.get(c, 0) for c in pt)
        agg_i, agg_g = Counter(), Counter()
        for c, v in pt.items():
            r = region.get(c, "?") or "?"
            agg_i[r] += v
            agg_g[r] += gen.get(c, 0)
        out(f"{p[0]}-{p[1]}:")
        for r in sorted(agg_i, key=lambda k: -agg_i[k]):
            bi, gi = 100.0 * agg_i[r] / tot_i, (100.0 * agg_g[r] / tot_g if tot_g else 0)
            ri = gi / bi if bi else 0
            per = 1000.0 * agg_g[r] / agg_i[r]
            out(f"  {r}: yuk %{bi:.1f} | genom %{gi:.1f} | RI={ri:.2f} | genom/1000 vaka={per:.2f}")
            reg_rows.append([f"{p[0]}-{p[1]}", r, f"{bi:.2f}", f"{gi:.2f}", f"{ri:.3f}", f"{per:.3f}"])
    with open(OUT / "regions.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["period", "who_region", "burden_share_pct", "genome_share_pct", "representation_index",
                    "genomes_per_1000"])
        w.writerows(reg_rows)

    # ------------------------------------------------ F. kategoriler, ulke tablosu, Lorenz
    out("")
    out(f"=== F. KAPSAMA KATEGORILERI ({MAIN[0]}-{MAIN[1]}; genom / 1.000 kumulatif vaka) ===")
    gen = gen_counts(samples, "c", *MAIN)
    pt = psum(inc, *MAIN)
    ti = sum(pt.values())
    tg = sum(gen.get(c, 0) for c in pt)
    cat_n, cat_i, cat_g = Counter(), Counter(), Counter()
    rows = []
    for c, v in pt.items():
        g = gen.get(c, 0)
        per = 1000.0 * g / v
        k = category(per)
        cat_n[k] += 1
        cat_i[k] += v
        cat_g[k] += g
        ri = ((100.0 * g / tg) / (100.0 * v / ti)) if tg else 0
        rows.append({"country": c, "who_region": region.get(c, ""), f"who_incidence_{MAIN[0]}_{MAIN[1]}": int(v),
                     "genomes": g, "genomes_per_1000": f"{per:.3f}", "burden_share_pct": f"{100.0 * v / ti:.3f}",
                     "genome_share_pct": f"{(100.0 * g / tg) if tg else 0:.3f}", "representation_index": f"{ri:.3f}",
                     "category": k})
    rows.sort(key=lambda r: -r[f"who_incidence_{MAIN[0]}_{MAIN[1]}"])
    with open(OUT / "country_representation.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for k in ["0", "<1", "1-10", ">=10"]:
        out(f"  {k}: {cat_n[k]} ulke | yukun %{100.0 * cat_i[k] / ti:.1f}'i | genomlarin %{(100.0 * cat_g[k] / tg) if tg else 0:.1f}'i")
    order = sorted(pt, key=lambda c: gen.get(c, 0) / pt[c])
    cw = cb = 0.0
    with open(OUT / "lorenz_main.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "cum_burden_share", "cum_genome_share"])
        w.writerow(["(origin)", 0, 0])
        for c in order:
            cb += pt[c] / ti
            cw += (gen.get(c, 0) / tg) if tg else 0
            w.writerow([c, f"{cb:.5f}", f"{cw:.5f}"])
    over = [r for r in rows if float(r["representation_index"]) > 5][:8]
    out("Asiri temsil edilen ulkeler (RI>5, yuke gore): " + ", ".join(
        f"{r['country']}={float(r['representation_index']):.0f}" for r in sorted(over, key=lambda r: -float(r['representation_index']))))

    # ------------------------------------------------ G. Turkiye
    out("")
    out("=== G. TURKIYE ===")
    zero = sorted([(v, c) for c, v in pt.items() if gen.get(c, 0) == 0], reverse=True)
    big_zero = [x for x in zero if x[0] >= 100000]
    tr = pt.get("Turkey", 0)
    out(f"Turkey: kumulatif tahmini vaka {int(tr)}, genom {gen.get('Turkey', 0)}")
    if tr and gen.get("Turkey", 0) == 0:
        rank = [c for _, c in zero].index("Turkey") + 1
        out(f"Sifir genomlu {len(zero)} ulke icinde yuke gore {rank}. sirada; "
            f"kumulatif vakasi >=100.000 olup genomu sifir olan ulke sayisi: {len(big_zero)}")
        out("Sifir genomlu ilk 10 ulke (kumulatif vaka): " + ", ".join(f"{c}={int(v)}" for v, c in zero[:10]))
    elif tr:
        per = 1000.0 * gen.get("Turkey", 0) / tr
        out(f"Genom/1000 vaka = {per:.3f}")

    text = "\n".join(LINES)
    (OUT / "summary.txt").write_text(text, encoding="utf-8")
    z = OUT / "mycotb_robustness_results.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in ["summary.txt", "release_cohorts.csv", "lag_matched.csv", "sensitivity.csv", "uncertainty.csv",
                  "regions.csv", "country_representation.csv", "lorenz_main.csv"]:
            p = OUT / n
            if p.exists():
                zf.write(p, n)
    print(f"\nBITTI. Bu dosyayi Claude'a yukleyin: {z}")


if __name__ == "__main__":
    main()
