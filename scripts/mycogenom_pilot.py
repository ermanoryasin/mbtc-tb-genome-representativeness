#!/usr/bin/env python3
"""
MycoGenom pilotu: halka acik Mycobacterium metadata'sinin yapilabilirlik testi.

Ne yapar (sadece METADATA indirir, ham okuma/FASTQ indirmez):
  1. NCBI SRA'da MTBC kayit sayisini kontrol eder (erisim testi + capraz kontrol).
  2. ENA'dan (INSDC aynasi = NCBI SRA/BioSample ile ayni kayitlar) Turkiye kaynakli
     tum Mycobacterium run ve sample kayitlarini ceker.
  3. ENA'dan tum MTBC (taxid 77643) read_run metadata'sini yil yil ceker.
  4. Ulke bazli ozet, yil ozeti, Turkiye ozeti ve metadata tamlik oranlarini hesaplar.
  5. (Opsiyonel) WHO TB tahmin tablosunu indirir.
  6. Gonderilecek kucuk dosyalari tek bir zip'e toplar.

Gereksinim: python3 (sadece standart kutuphane) ve curl (macOS'ta hazir).
Calistirma:  cd ~/Downloads && python3 mycogenom_pilot.py
Kesilirse ayni komutla tekrar calistirin; inen parcalar atlanir.
"""
import csv
import datetime
import json
import os
import re
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(10 ** 8)

ENA = "https://www.ebi.ac.uk/ena/portal/api"
OUT = Path(os.environ.get("MYCO_OUT", str(Path.home() / "Downloads" / "mycogenom_pilot_out")))

MISSING = {
    "", "missing", "not applicable", "not collected", "not provided", "na", "n/a",
    "unknown", "none", "restricted access", "not available", "null", "-", "not determined",
    "not recorded", "unspecified", "not known",
}

COUNTRY_FIX = {
    "turkey": "Türkiye", "turkiye": "Türkiye", "türkiye": "Türkiye", "türkei": "Türkiye",
    "usa": "USA", "united states": "USA", "united states of america": "USA", "us": "USA",
    "uk": "United Kingdom", "great britain": "United Kingdom", "england": "United Kingdom",
    "viet nam": "Vietnam", "russian federation": "Russia",
    "republic of korea": "South Korea", "korea": "South Korea",
    "people's republic of china": "China", "czech republic": "Czechia", "burma": "Myanmar",
    "iran, islamic republic of": "Iran", "islamic republic of iran": "Iran",
    "syrian arab republic": "Syria", "republic of georgia": "Georgia",
}

RUN_FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession", "study_accession",
    "experiment_accession", "scientific_name", "tax_id", "country", "collection_date",
    "first_public", "first_created", "host", "host_tax_id", "isolation_source",
    "instrument_platform", "instrument_model", "library_strategy", "library_source",
    "library_layout", "read_count", "base_count", "fastq_bytes", "center_name", "strain",
]
RUN_FIELDS_TURKEY_EXTRA = ["study_title"]
SAMPLE_FIELDS = [
    "sample_accession", "secondary_sample_accession", "scientific_name", "tax_id", "country",
    "collection_date", "host", "isolation_source", "strain", "sample_title", "center_name",
    "first_public",
]


# ---------------------------------------------------------------- yardimcilar
def log(msg):
    line = f"[{datetime.datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    try:
        with open(OUT / "run_log.txt", "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except OSError:
        pass


def curl(url, params=None, dest=None, timeout=900):
    cmd = ["curl", "-sS", "--fail", "--retry", "3", "--retry-delay", "5",
           "--max-time", str(timeout)]
    target = None
    if params:
        cmd.append("-G")
        for k, v in params.items():
            cmd += ["--data-urlencode", f"{k}={v}"]
    if dest:
        target = Path(str(dest) + ".part")
        cmd += ["-o", str(target)]
    cmd.append(url)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        if target and target.exists():
            target.unlink()
        raise RuntimeError(f"curl cikis kodu {r.returncode}: {r.stderr.strip()[:300]}")
    if target:
        target.replace(dest)
    return r.stdout


def valid_fields(result, wanted):
    try:
        txt = curl(f"{ENA}/returnFields", {"result": result, "format": "json"}, timeout=120)
        ok = {d["columnId"] for d in json.loads(txt)}
        use = [f for f in wanted if f in ok]
        dropped = [f for f in wanted if f not in ok]
        if dropped:
            log(f"{result}: ENA'da olmayan alanlar atlandi: {', '.join(dropped)}")
        return use
    except Exception as e:  # noqa: BLE001
        log(f"{result}: alan dogrulamasi yapilamadi ({e}); liste oldugu gibi kullaniliyor")
        return wanted


def ena_download(result, query, fields, dest):
    params = {"result": result, "query": query, "fields": ",".join(fields),
              "format": "tsv", "limit": "0"}
    curl(f"{ENA}/search", params, dest=dest)


def count_data_lines(path):
    if not Path(path).exists():
        return 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)


def present(v):
    return (v or "").strip().lower() not in MISSING


def norm_country(raw):
    c = (raw or "").strip()
    if c.lower() in MISSING:
        return "(missing)"
    c = c.split(":")[0].strip()
    if c.lower() in MISSING:
        return "(missing)"
    return COUNTRY_FIX.get(c.lower(), c)


def year_of(v, n=4):
    s = (v or "").strip()
    return s[:n] if len(s) >= n and s[:n].isdigit() else "(missing)"


def pct(a, b):
    return f"{(100.0 * a / b):.1f}" if b else "0.0"


# ---------------------------------------------------------------- asamalar
def ncbi_check(lines):
    try:
        xml = curl("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                   {"db": "sra", "term": "txid77643[Organism:exp]", "rettype": "count"},
                   timeout=60)
        m = re.search(r"<Count>(\d+)</Count>", xml)
        if m:
            lines.append(f"NCBI SRA (esearch) MTBC kayit sayisi: {m.group(1)}")
            log(f"NCBI erisimi calisiyor; SRA MTBC kayit sayisi: {m.group(1)}")
        else:
            lines.append("NCBI SRA esearch yaniti beklenen formatta degil.")
    except Exception as e:  # noqa: BLE001
        lines.append(f"NCBI E-utilities erisimi basarisiz: {e}")
        log(f"NCBI E-utilities erisimi basarisiz (kritik degil): {e}")


def download_turkey(lines):
    run_fields = valid_fields("read_run", RUN_FIELDS + RUN_FIELDS_TURKEY_EXTRA)
    sample_fields = valid_fields("sample", SAMPLE_FIELDS)
    cond = '(country="Turkey*" OR country="Turkiye*" OR country="Türkiye*" OR country="Turkei*")'
    q = f"tax_tree(1763) AND {cond}"
    runs = OUT / "turkey_myco_runs.tsv"
    samples = OUT / "turkey_myco_samples.tsv"
    for result, fields, dest in (("read_run", run_fields, runs), ("sample", sample_fields, samples)):
        if dest.exists() and dest.stat().st_size > 0:
            log(f"{dest.name} zaten var, atlandi")
        else:
            log(f"Turkiye {result} kayitlari cekiliyor...")
            try:
                ena_download(result, q, fields, dest)
            except Exception as e:  # noqa: BLE001
                log(f"HATA: Turkiye {result} indirilemedi: {e}")
                lines.append(f"Turkiye {result} indirilemedi: {e}")
                continue
        log(f"  {dest.name}: {count_data_lines(dest)} satir")
    return runs, samples


def download_mtbc(lines):
    fields = valid_fields("read_run", RUN_FIELDS)
    (OUT / "raw").mkdir(exist_ok=True)
    parts = []
    for y in range(2008, datetime.date.today().year + 1):
        dest = OUT / "raw" / f"mtbc_runs_{y}.tsv"
        if not (dest.exists() and dest.stat().st_size > 0):
            q = f"tax_tree(77643) AND first_public>={y}-01-01 AND first_public<={y}-12-31"
            log(f"MTBC {y} cekiliyor...")
            try:
                ena_download("read_run", q, fields, dest)
            except Exception as e:  # noqa: BLE001
                log(f"HATA: MTBC {y} indirilemedi: {e}")
                lines.append(f"MTBC {y} indirilemedi: {e}")
                continue
        log(f"  {y}: {count_data_lines(dest)} run")
        parts.append(dest)
    merged = OUT / "mtbc_read_runs_all.tsv"
    header_written = False
    with open(merged, "w", encoding="utf-8", newline="") as out:
        for p in parts:
            with open(p, "r", encoding="utf-8", errors="replace") as f:
                first = True
                for line in f:
                    if first:
                        first = False
                        if not header_written:
                            out.write(line)
                            header_written = True
                        continue
                    out.write(line)
    log(f"Birlestirildi: {merged.name} ({count_data_lines(merged)} run)")
    return merged


def analyze_global(merged, lines):
    if not merged.exists() or merged.stat().st_size == 0:
        lines.append("Global MTBC dosyasi bos/yok; analiz yapilamadi.")
        return
    runs = Counter()
    wgs_runs = Counter()
    samples = defaultdict(set)
    wgs_samples = defaultdict(set)
    d_ok, h_ok, i_ok = Counter(), Counter(), Counter()
    years = Counter()
    years_missing_country = Counter()
    strategy = Counter()
    platform = Counter()
    total = 0
    all_samples = set()
    with open(merged, "r", encoding="utf-8", errors="replace", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            c = norm_country(row.get("country"))
            runs[c] += 1
            sa = row.get("sample_accession") or row.get("run_accession")
            samples[c].add(sa)
            all_samples.add(sa)
            is_wgs = (row.get("library_strategy") or "").upper() == "WGS"
            if is_wgs:
                wgs_runs[c] += 1
                wgs_samples[c].add(sa)
            if present(row.get("collection_date")):
                d_ok[c] += 1
            if present(row.get("host")):
                h_ok[c] += 1
            if present(row.get("isolation_source")):
                i_ok[c] += 1
            y = year_of(row.get("first_public"))
            years[y] += 1
            if c == "(missing)":
                years_missing_country[y] += 1
            strategy[(row.get("library_strategy") or "(missing)")] += 1
            platform[(row.get("instrument_platform") or "(missing)")] += 1

    with open(OUT / "country_summary.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "runs_all", "samples_all", "runs_wgs", "samples_wgs",
                    "pct_runs_with_collection_date", "pct_runs_with_host",
                    "pct_runs_with_isolation_source"])
        for c, n in runs.most_common():
            w.writerow([c, n, len(samples[c]), wgs_runs[c], len(wgs_samples[c]),
                        pct(d_ok[c], n), pct(h_ok[c], n), pct(i_ok[c], n)])
    with open(OUT / "year_summary.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["first_public_year", "runs", "runs_missing_country", "pct_missing_country"])
        for y in sorted(years):
            w.writerow([y, years[y], years_missing_country[y], pct(years_missing_country[y], years[y])])

    miss = runs.get("(missing)", 0)
    lines.append("")
    lines.append("=== GLOBAL MTBC (ENA read_run, taxid 77643 agaci) ===")
    lines.append(f"Toplam run: {total} | benzersiz sample: {len(all_samples)}")
    lines.append(f"Ulke bilgisi eksik/gecersiz run: {miss} (%{pct(miss, total)})")
    lines.append(f"Ulke sayisi (eksik haric): {len([c for c in runs if c != '(missing)'])}")
    lines.append("library_strategy: " + ", ".join(f"{k}={v}" for k, v in strategy.most_common(6)))
    lines.append("platform: " + ", ".join(f"{k}={v}" for k, v in platform.most_common(5)))
    lines.append("En cok run'a sahip 15 ulke (run / sample / WGS sample):")
    for c, n in runs.most_common(16):
        lines.append(f"  {c}: {n} / {len(samples[c])} / {len(wgs_samples[c])}")
    lines.append("Turkiye ve cevresi (run / sample / WGS sample / %tarih / %konak):")
    wanted = ["türkiye", "georgia", "greece", "iran", "iraq", "syria", "bulgaria", "armenia",
              "azerbaijan", "russia", "egypt", "romania", "lebanon", "jordan", "saudi arabia"]
    for c in runs:
        if c.lower() in wanted:
            lines.append(f"  {c}: {runs[c]} / {len(samples[c])} / {len(wgs_samples[c])} / "
                         f"%{pct(d_ok[c], runs[c])} / %{pct(h_ok[c], runs[c])}")


def analyze_turkey(runs_path, samples_path, lines):
    lines.append("")
    lines.append("=== TURKIYE KAYNAKLI MYCOBACTERIUM (genus, taxid 1763 agaci) ===")
    if runs_path.exists() and runs_path.stat().st_size > 0:
        rows = []
        with open(runs_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        n = len(rows)
        lines.append(f"Run sayisi: {n} | benzersiz sample: "
                     f"{len({r.get('sample_accession') for r in rows})} | benzersiz study: "
                     f"{len({r.get('study_accession') for r in rows})}")
        sp = Counter(r.get("scientific_name") or "(missing)" for r in rows)
        lines.append("Tur dagilimi (run): " + ", ".join(f"{k}={v}" for k, v in sp.most_common(12)))
        st = Counter(r.get("library_strategy") or "(missing)" for r in rows)
        lines.append("library_strategy: " + ", ".join(f"{k}={v}" for k, v in st.most_common(6)))
        pf = Counter(r.get("instrument_platform") or "(missing)" for r in rows)
        lines.append("platform: " + ", ".join(f"{k}={v}" for k, v in pf.most_common(5)))
        cy = Counter(year_of(r.get("collection_date")) for r in rows)
        lines.append("collection_date yili: " + ", ".join(f"{k}={cy[k]}" for k in sorted(cy)))
        py = Counter(year_of(r.get("first_public")) for r in rows)
        lines.append("first_public yili: " + ", ".join(f"{k}={py[k]}" for k in sorted(py)))
        for fld in ("collection_date", "host", "isolation_source", "strain"):
            ok = sum(1 for r in rows if present(r.get(fld)))
            lines.append(f"  {fld} dolu: %{pct(ok, n)}")
        mtbc_like = sum(v for k, v in sp.items() if "tuberculosis" in k.lower()
                        or "bovis" in k.lower() or "africanum" in k.lower()
                        or "caprae" in k.lower() or "complex" in k.lower())
        lines.append(f"MTBC benzeri run (ad eslesmesi): {mtbc_like} | digerleri (NTM vb.): {n - mtbc_like}")
        studies = defaultdict(lambda: {"runs": 0, "title": "", "center": "", "species": Counter()})
        for r in rows:
            s = studies[r.get("study_accession") or "(missing)"]
            s["runs"] += 1
            s["title"] = s["title"] or (r.get("study_title") or "")
            s["center"] = s["center"] or (r.get("center_name") or "")
            s["species"][r.get("scientific_name") or "(missing)"] += 1
        with open(OUT / "turkey_studies.csv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["study_accession", "runs", "center_name", "species", "study_title"])
            for acc, s in sorted(studies.items(), key=lambda kv: -kv[1]["runs"]):
                w.writerow([acc, s["runs"], s["center"],
                            "; ".join(f"{k}={v}" for k, v in s["species"].most_common(3)),
                            s["title"][:200]])
        lines.append("En buyuk 8 study: " + "; ".join(
            f"{a} ({s['runs']} run)" for a, s in
            sorted(studies.items(), key=lambda kv: -kv[1]["runs"])[:8]))
    else:
        lines.append("Turkiye run dosyasi yok/bos.")
    if samples_path.exists() and samples_path.stat().st_size > 0:
        with open(samples_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            srows = list(csv.DictReader(f, delimiter="\t"))
        ssp = Counter(r.get("scientific_name") or "(missing)" for r in srows)
        lines.append(f"Sample (BioSample) kaydi: {len(srows)} | tur: "
                     + ", ".join(f"{k}={v}" for k, v in ssp.most_common(10)))
    else:
        lines.append("Turkiye sample dosyasi yok/bos.")


def download_who(lines):
    dest = OUT / "who_tb_estimates.csv"
    if dest.exists() and dest.stat().st_size > 0:
        return
    try:
        curl("https://extranet.who.int/tme/generateCSV.asp?ds=estimates", dest=dest, timeout=180)
        with open(dest, "r", encoding="utf-8", errors="replace") as f:
            head = f.read(300)
        if "<html" in head.lower() or head.lstrip().startswith("<"):
            dest.unlink()
            raise RuntimeError("CSV yerine HTML dondu")
        log("WHO TB tahmin tablosu indirildi")
    except Exception as e:  # noqa: BLE001
        log(f"WHO tablosu indirilemedi (kritik degil, elle indirilebilir): {e}")
        lines.append(f"WHO tablosu indirilemedi: {e}")


def make_zip():
    names = ["summary.txt", "country_summary.csv", "year_summary.csv", "turkey_myco_runs.tsv",
             "turkey_myco_samples.tsv", "turkey_studies.csv", "who_tb_estimates.csv", "run_log.txt"]
    z = OUT / "mycogenom_pilot_results.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in names:
            p = OUT / n
            if p.exists():
                zf.write(p, n)
    return z


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"Cikti klasoru: {OUT}")
    lines = [f"MycoGenom pilotu - {datetime.datetime.now():%Y-%m-%d %H:%M}"]
    ncbi_check(lines)
    t_runs, t_samples = download_turkey(lines)
    merged = download_mtbc(lines)
    analyze_global(merged, lines)
    analyze_turkey(t_runs, t_samples, lines)
    download_who(lines)
    text = "\n".join(lines)
    (OUT / "summary.txt").write_text(text, encoding="utf-8")
    print("\n" + text + "\n")
    z = make_zip()
    log(f"BITTI. Bu dosyayi Claude'a yukleyin: {z}")


if __name__ == "__main__":
    sys.exit(main())
