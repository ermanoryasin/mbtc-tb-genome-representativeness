#!/usr/bin/env python3
"""
MycoTB pilotu (v2): yalnizca Mycobacterium tuberculosis odakli metadata taramasi.

Ne yapar (sadece METADATA indirir, FASTQ indirmez):
  1. NCBI SRA erisim testi (capraz kontrol sayisi).
  2. Turkiye kaynakli MTBC kayitlari: country alani Turkey/TURKEY/Turkiye/Türkiye vb.
     tum yazim varyantlariyla (run + BioSample).
  3. Global M. tuberculosis (taxid 1773 agaci) read_run metadata'si, yil yil.
  4. Ulke bilgisi EKSIK kayitlarda merkez adi / study basligi uzerinden Turk kurumlari taramasi
     ("gizli Turkiye adaylari").
  5. Ulke bazli ozet (siki etiket filtresi ile), yil ozeti, WHO TB yukune gore
     genom/vaka adalet metrikleri.
  6. Gonderilecek kucuk dosyalari tek zip'e toplar.

Gereksinim: python3 (standart kutuphane) ve curl (macOS'ta hazir).
Calistirma:  cd ~/Downloads && python3 mycotb_pilot.py
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
MBTC_BASE = Path(os.environ.get("MBTC_BASE", str(Path(__file__).resolve().parent.parent / "data")))

csv.field_size_limit(10 ** 8)

ENA = "https://www.ebi.ac.uk/ena/portal/api"
OUT = Path(os.environ.get("MYCO_OUT", str(MBTC_BASE / "mycotb_pilot_out")))
V1_MERGED = MBTC_BASE / "mycogenom_pilot_out" / "mtbc_read_runs_all.tsv"

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
    "syrian arab republic": "Syria", "republic of georgia": "Georgia",
    "swaziland": "Eswatini",
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
}

RUN_FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession", "study_accession",
    "experiment_accession", "scientific_name", "tax_id", "country", "collection_date",
    "first_public", "first_created", "host", "host_tax_id", "isolation_source",
    "instrument_platform", "instrument_model", "library_strategy", "library_source",
    "library_layout", "read_count", "base_count", "fastq_bytes", "center_name", "strain",
    "study_title",
]
SAMPLE_FIELDS = [
    "sample_accession", "secondary_sample_accession", "scientific_name", "tax_id", "country",
    "collection_date", "host", "isolation_source", "strain", "sample_title", "center_name",
    "first_public",
]

# Turkiye'deki kurum/sehir anahtar kelimeleri (merkez adi + study basligi taramasi)
TR_KW = [
    "turkey", "turkiye", "türkiye", "ankara", "istanbul", "stanbul", "izmir", "bursa", "konya",
    "samsun", "antalya", "adana", "trabzon", "erzurum", "diyarbakir", "diyarbakır", "gaziantep",
    "kayseri", "hacettepe", "marmara", "ondokuz", "dokuz eylul", "dokuz eylül", "gazi university",
    "erciyes", "ataturk", "atatürk", "bogazici", "boğaziçi", "sabanci", "sabancı",
    "middle east technical", "yeditepe", "bezmialem", "medeniyet", "cerrahpasa", "cerrahpaşa",
    "cukurova", "çukurova", "akdeniz university", "pamukkale", "sakarya", "kocaeli", "inonu",
    "inönü", "firat university", "fırat university", "dicle", "uludag", "uludağ", "trakya",
    "eskisehir", "eskişehir", "anadolu", "tubitak", "tübitak", "refik saydam", "koc university",
    "koç university", "sureyyapasa", "süreyyapaşa", "ege university", "ege üniversitesi",
]
TR_RE = re.compile(r"(?<![a-z])(" + "|".join(re.escape(k) for k in TR_KW) + r")(?![a-z])",
                   re.IGNORECASE)

TURKEY_COUNTRY_Q = (
    '(country="Turkey*" OR country="TURKEY*" OR country="turkey*" OR country="Türkiye*" '
    'OR country="TÜRKİYE*" OR country="Turkiye*" OR country="TURKIYE*" '
    'OR country="Republic of Turkey*")'
)

NON_STRICT_TOKENS = ("variant", "complex", "bovis", "caprae", "africanum", "microti", "canettii",
                     "orygis", "pinnipedii", "mungi", "suricattae", "bcg")


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


def year_of(v, n=4):
    s = (v or "").strip()
    return s[:n] if len(s) >= n and s[:n].isdigit() else "(missing)"


def pct(a, b):
    return f"{(100.0 * a / b):.1f}" if b else "0.0"


# ---------------------------------------------------------------- indirme
def ncbi_check(lines):
    try:
        xml = curl("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                   {"db": "sra", "term": "txid1773[Organism:exp]", "rettype": "count"}, timeout=60)
        m = re.search(r"<Count>(\d+)</Count>", xml)
        if m:
            lines.append(f"NCBI SRA (esearch) M. tuberculosis (txid1773 agaci) kayit sayisi: {m.group(1)}")
            log(f"NCBI erisimi calisiyor; SRA M. tuberculosis kayit sayisi: {m.group(1)}")
    except Exception as e:  # noqa: BLE001
        lines.append(f"NCBI E-utilities erisimi basarisiz: {e}")
        log(f"NCBI E-utilities erisimi basarisiz (kritik degil): {e}")


def download_turkey(lines):
    run_fields = valid_fields("read_run", RUN_FIELDS)
    sample_fields = valid_fields("sample", SAMPLE_FIELDS)
    q = f"tax_tree(77643) AND {TURKEY_COUNTRY_Q}"
    runs = OUT / "turkey_mtbc_runs.tsv"
    samples = OUT / "turkey_mtbc_samples.tsv"
    for result, fields, dest in (("read_run", run_fields, runs), ("sample", sample_fields, samples)):
        if dest.exists() and dest.stat().st_size > 0:
            log(f"{dest.name} zaten var, atlandi")
        else:
            log(f"Turkey/TURKEY varyantlariyla MTBC {result} kayitlari cekiliyor...")
            try:
                ena_download(result, q, fields, dest)
            except Exception as e:  # noqa: BLE001
                log(f"HATA: Turkiye {result} indirilemedi: {e}")
                lines.append(f"Turkiye {result} indirilemedi: {e}")
                continue
        log(f"  {dest.name}: {count_data_lines(dest)} satir")
    return runs, samples


def download_tb(lines):
    fields = valid_fields("read_run", RUN_FIELDS)
    (OUT / "raw").mkdir(exist_ok=True)
    parts = []
    for y in range(2008, datetime.date.today().year + 1):
        dest = OUT / "raw" / f"mtb_runs_{y}.tsv"
        if not (dest.exists() and dest.stat().st_size > 0):
            q = f"tax_tree(1773) AND first_public>={y}-01-01 AND first_public<={y}-12-31"
            log(f"M. tuberculosis {y} cekiliyor...")
            try:
                ena_download("read_run", q, fields, dest)
            except Exception as e:  # noqa: BLE001
                log(f"HATA: {y} indirilemedi: {e}")
                lines.append(f"M. tuberculosis {y} indirilemedi: {e}")
                continue
        log(f"  {y}: {count_data_lines(dest)} run")
        parts.append(dest)
    merged = OUT / "mtb_read_runs_all.tsv"
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


def download_who(lines):
    dest = OUT / "who_tb_estimates.csv"
    old = MBTC_BASE / "mycogenom_pilot_out" / "who_tb_estimates.csv"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        if old.exists() and old.stat().st_size > 0:
            dest.write_bytes(old.read_bytes())
            log("WHO tablosu onceki pilottan kopyalandi")
            return dest
        curl("https://extranet.who.int/tme/generateCSV.asp?ds=estimates", dest=dest, timeout=180)
        with open(dest, "r", encoding="utf-8", errors="replace") as f:
            head = f.read(300)
        if head.lstrip().startswith("<"):
            dest.unlink()
            raise RuntimeError("CSV yerine HTML dondu")
        log("WHO TB tahmin tablosu indirildi")
        return dest
    except Exception as e:  # noqa: BLE001
        log(f"WHO tablosu indirilemedi (kritik degil): {e}")
        lines.append(f"WHO tablosu indirilemedi: {e}")
        return None


# ---------------------------------------------------------------- analiz
class Acc:
    def __init__(self):
        self.runs = Counter()
        self.samples = defaultdict(set)
        self.wgs_samples = defaultdict(set)
        self.d = Counter()
        self.h = Counter()
        self.i = Counter()
        self.years = Counter()
        self.years_missing = Counter()

    def add(self, c, sa, is_wgs, row, y):
        self.runs[c] += 1
        self.samples[c].add(sa)
        if is_wgs:
            self.wgs_samples[c].add(sa)
        if present(row.get("collection_date")):
            self.d[c] += 1
        if present(row.get("host")):
            self.h[c] += 1
        if present(row.get("isolation_source")):
            self.i[c] += 1
        self.years[y] += 1
        if c == "(missing)":
            self.years_missing[y] += 1


def analyze_global(merged, lines):
    if not merged.exists() or merged.stat().st_size == 0:
        lines.append("Global M. tuberculosis dosyasi bos/yok; analiz yapilamadi.")
        return {}, {}
    tree, strict = Acc(), Acc()
    names = Counter()
    strategy, platform = Counter(), Counter()
    hidden = []
    total = 0
    with open(merged, "r", encoding="utf-8", errors="replace", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            c = norm_country(row.get("country"))
            name = row.get("scientific_name") or "(missing)"
            names[name] += 1
            sa = row.get("sample_accession") or row.get("run_accession")
            is_wgs = (row.get("library_strategy") or "").upper() == "WGS"
            y = year_of(row.get("first_public"))
            tree.add(c, sa, is_wgs, row, y)
            if is_strict_mtb(name):
                strict.add(c, sa, is_wgs, row, y)
            strategy[(row.get("library_strategy") or "(missing)")] += 1
            platform[(row.get("instrument_platform") or "(missing)")] += 1
            if c != "Turkey":
                txt = f"{row.get('center_name') or ''} | {row.get('study_title') or ''}"
                m = TR_RE.search(txt)
                if m:
                    cat = ("missing_country" if c in ("(missing)", "(multiple)")
                           else f"other_country:{c}")
                    hidden.append({
                        "category": cat, "keyword": m.group(1).lower(),
                        "run_accession": row.get("run_accession", ""),
                        "sample_accession": row.get("sample_accession", ""),
                        "study_accession": row.get("study_accession", ""),
                        "country_raw": row.get("country", ""),
                        "scientific_name": name,
                        "library_strategy": row.get("library_strategy", ""),
                        "collection_date": row.get("collection_date", ""),
                        "first_public": row.get("first_public", ""),
                        "center_name": row.get("center_name", ""),
                        "study_title": (row.get("study_title") or "")[:200],
                    })

    # country_summary
    countries = sorted(tree.runs, key=lambda k: -len(strict.wgs_samples[k]))
    with open(OUT / "country_summary.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country", "runs_tree1773", "samples_tree1773", "wgs_samples_tree1773",
                    "runs_strict", "samples_strict", "wgs_samples_strict",
                    "pct_strict_runs_with_collection_date", "pct_strict_runs_with_host",
                    "pct_strict_runs_with_isolation_source"])
        for c in countries:
            n = strict.runs[c]
            w.writerow([c, tree.runs[c], len(tree.samples[c]), len(tree.wgs_samples[c]),
                        n, len(strict.samples[c]), len(strict.wgs_samples[c]),
                        pct(strict.d[c], n), pct(strict.h[c], n), pct(strict.i[c], n)])
    with open(OUT / "year_summary.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["first_public_year", "runs_tree", "pct_missing_country_tree",
                    "runs_strict", "pct_missing_country_strict"])
        for y in sorted(set(tree.years) | set(strict.years)):
            w.writerow([y, tree.years[y], pct(tree.years_missing[y], tree.years[y]),
                        strict.years[y], pct(strict.years_missing[y], strict.years[y])])
    if hidden:
        with open(OUT / "turkey_hidden_candidates.csv", "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(hidden[0].keys()))
            w.writeheader()
            w.writerows(hidden)

    tr = sum(tree.runs.values())
    sr = sum(strict.runs.values())
    lines.append("")
    lines.append("=== GLOBAL M. tuberculosis (ENA read_run, taxid 1773 agaci) ===")
    lines.append(f"Toplam run: {tr} | benzersiz sample: {len({s for c in tree.samples for s in tree.samples[c]})}")
    lines.append(f"Siki etiket filtresi (adi 'Mycobacterium tuberculosis', varyant/complex vb. haric): "
                 f"{sr} run (%{pct(sr, tr)})")
    lines.append("scientific_name dagilimi: " + ", ".join(f"{k}={v}" for k, v in names.most_common(8)))
    if V1_MERGED.exists():
        lines.append(f"Karsilastirma: onceki pilot tum MTBC (taxid 77643) run = {count_data_lines(V1_MERGED)}")
    for label, acc in (("tree", tree), ("strict", strict)):
        tot = sum(acc.runs.values())
        miss = acc.runs.get("(missing)", 0)
        mult = acc.runs.get("(multiple)", 0)
        lines.append(f"[{label}] ulke eksik: {miss} run (%{pct(miss, tot)}), cok-ulkeli: {mult}; "
                     f"ulke sayisi: {len([c for c in acc.runs if c not in ('(missing)', '(multiple)')])}")
    lines.append("library_strategy: " + ", ".join(f"{k}={v}" for k, v in strategy.most_common(6)))
    lines.append("platform: " + ", ".join(f"{k}={v}" for k, v in platform.most_common(5)))
    lines.append("En cok WGS sample'a sahip 15 ulke (siki filtre: run / sample / WGS sample):")
    for c in countries[:16]:
        lines.append(f"  {c}: {strict.runs[c]} / {len(strict.samples[c])} / {len(strict.wgs_samples[c])}")
    lines.append("Turkiye ve cevresi (siki filtre: run / sample / WGS sample / %tarih / %konak):")
    wanted = {"turkey", "georgia", "greece", "iran", "iraq", "syria", "bulgaria", "armenia",
              "azerbaijan", "russia", "egypt", "romania", "lebanon"}
    for c in countries:
        if c.lower() in wanted:
            n = strict.runs[c]
            lines.append(f"  {c}: {n} / {len(strict.samples[c])} / {len(strict.wgs_samples[c])} / "
                         f"%{pct(strict.d[c], n)} / %{pct(strict.h[c], n)}")
    if "Turkey" not in tree.runs:
        lines.append("  Turkey: global dosyada ulke etiketi 'Turkey' olan kayit yok")

    # gizli Turkiye adaylari
    lines.append("")
    lines.append("=== GIZLI TURKIYE ADAYLARI (merkez adi / study basligi taramasi) ===")
    if hidden:
        cats = Counter("missing_country" if h["category"] == "missing_country" else "other_country"
                       for h in hidden)
        lines.append(f"Aday run: {len(hidden)} | ulkesi eksik: {cats['missing_country']} | "
                     f"baska ulke etiketli: {cats['other_country']}")
        miss_h = [h for h in hidden if h["category"] == "missing_country"]
        lines.append(f"Ulkesi eksik adaylarda benzersiz sample: {len({h['sample_accession'] for h in miss_h})}, "
                     f"study: {len({h['study_accession'] for h in miss_h})}")
        kw = Counter(h["keyword"] for h in miss_h)
        lines.append("Anahtar kelimeler (ulkesi eksik): " + ", ".join(f"{k}={v}" for k, v in kw.most_common(10)))
        st = defaultdict(lambda: [0, "", ""])
        for h in miss_h:
            s = st[h["study_accession"]]
            s[0] += 1
            s[1] = s[1] or h["center_name"]
            s[2] = s[2] or h["study_title"]
        for acc, (n, cen, ti) in sorted(st.items(), key=lambda kv: -kv[1][0])[:8]:
            lines.append(f"  {acc}: {n} run | {cen[:50]} | {ti[:90]}")
        lines.append("(Not: bunlar otomatik aday. Turkiye kaynakli oldugu elle dogrulanmalidir.)")
    else:
        lines.append("Aday bulunamadi.")
    hid_miss_samples = len({h["sample_accession"] for h in hidden if h["category"] == "missing_country"})
    return ({c: len(strict.wgs_samples[c]) for c in strict.wgs_samples},
            {"hidden_missing_samples": hid_miss_samples})


def analyze_turkey(runs_path, samples_path, lines):
    lines.append("")
    lines.append("=== TURKEY ETIKETLI MTBC KAYITLARI (taxid 77643 agaci, tum yazim varyantlari) ===")
    if runs_path.exists() and runs_path.stat().st_size > 0:
        with open(runs_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        n = len(rows)
        raw_vals = Counter((r.get("country") or "").strip() for r in rows)
        lines.append("Ham country yazimlari: " + ", ".join(f"'{k}'={v}" for k, v in raw_vals.most_common(8)))
        lines.append(f"Run: {n} | sample: {len({r.get('sample_accession') for r in rows})} | "
                     f"study: {len({r.get('study_accession') for r in rows})}")
        sp = Counter(r.get("scientific_name") or "(missing)" for r in rows)
        lines.append("Tur dagilimi (run): " + ", ".join(f"{k}={v}" for k, v in sp.most_common(10)))
        strict_n = sum(1 for r in rows if is_strict_mtb(r.get("scientific_name")))
        strict_s = len({r.get("sample_accession") for r in rows if is_strict_mtb(r.get("scientific_name"))})
        lines.append(f"Siki M. tuberculosis: {strict_n} run, {strict_s} sample")
        st = Counter(r.get("library_strategy") or "(missing)" for r in rows)
        lines.append("library_strategy: " + ", ".join(f"{k}={v}" for k, v in st.most_common(6)))
        cy = Counter(year_of(r.get("collection_date")) for r in rows)
        lines.append("collection_date yili: " + ", ".join(f"{k}={cy[k]}" for k in sorted(cy)))
        py = Counter(year_of(r.get("first_public")) for r in rows)
        lines.append("first_public yili: " + ", ".join(f"{k}={py[k]}" for k in sorted(py)))
        for fld in ("collection_date", "host", "isolation_source"):
            ok = sum(1 for r in rows if present(r.get(fld)))
            lines.append(f"  {fld} dolu: %{pct(ok, n)}")
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
        lines.append("Study'ler: " + "; ".join(
            f"{a} ({s['runs']} run)" for a, s in
            sorted(studies.items(), key=lambda kv: -kv[1]["runs"])[:8]))
    else:
        lines.append("Turkiye run dosyasi yok/bos.")
    if samples_path.exists() and samples_path.stat().st_size > 0:
        with open(samples_path, "r", encoding="utf-8", errors="replace", newline="") as f:
            srows = list(csv.DictReader(f, delimiter="\t"))
        ssp = Counter(r.get("scientific_name") or "(missing)" for r in srows)
        lines.append(f"BioSample kaydi: {len(srows)} | " + ", ".join(f"{k}={v}" for k, v in ssp.most_common(8)))


def equity(wgs_by_country, who_path, extra, lines):
    if not who_path or not Path(who_path).exists():
        lines.append("WHO tablosu yok; adalet metrikleri hesaplanamadi.")
        return
    with open(who_path, "r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    year = max(int(r["year"]) for r in rows if r.get("e_inc_num"))
    inc = {r["country"]: float(r["e_inc_num"]) for r in rows
           if int(r["year"]) == year and r.get("e_inc_num") and float(r["e_inc_num"]) > 0}
    table = []
    for c, i in inc.items():
        k = WHO_MAP.get(c, c)
        w = wgs_by_country.get(k, 0)
        table.append((c, k, i, w))
    ti = sum(t[2] for t in table)
    tw = sum(t[3] for t in table)
    with open(OUT / "equity_table.csv", "w", encoding="utf-8", newline="") as f:
        w_ = csv.writer(f)
        w_.writerow(["who_country", "ena_country", f"who_incidence_{year}", "wgs_samples_strict",
                     "genomes_per_1000_cases"])
        for c, k, i, w in sorted(table, key=lambda t: -t[2]):
            w_.writerow([c, k, int(i), w, f"{1000 * w / i:.2f}"])
    lo = sum(i for c, k, i, w in table if 1000 * w / i < 10)
    rs = sorted(table, key=lambda t: t[3] / t[2])
    ci = cw = area = 0.0
    for c, k, i, w in rs:
        pi, pw = i / ti, (w / tw if tw else 0)
        area += pi * (cw + pw / 2)
        ci += pi
        cw += pw
    top = sorted(table, key=lambda t: -t[2])[:30]
    lines.append("")
    lines.append(f"=== ADALET METRIKLERI (WHO {year} insidansi; siki M. tuberculosis WGS sample) ===")
    lines.append(f"Ulke: {len(table)} | genom/1000 yillik vaka (kuresel): {1000 * tw / ti:.1f}")
    lines.append(f"Yukun %{100 * lo / ti:.1f}'i yillik 1000 vakaya <10 genom dusen ulkelerde")
    lines.append(f"Gini (genom dagilimi vs yuk): {1 - 2 * area:.3f}")
    lines.append(f"Ilk 30 yuk ulkesi: yukun %{100 * sum(t[2] for t in top) / ti:.1f}'i, "
                 f"WGS'nin %{100 * sum(t[3] for t in top) / tw:.1f}'i" if tw else "WGS yok")
    names = ["Türkiye", "Georgia", "Greece", "Iran (Islamic Republic of)", "Iraq",
             "Syrian Arab Republic", "Bulgaria", "Armenia", "Azerbaijan", "Russian Federation",
             "Romania", "Lebanon", "Egypt"]
    lines.append("Komsular (insidans / WGS sample / genom per 1000 vaka):")
    for c, k, i, w in table:
        if c in names:
            lines.append(f"  {c}: {i:.0f} / {w} / {1000 * w / i:.1f}")
    hid = extra.get("hidden_missing_samples", 0)
    for c, k, i, w in table:
        if k == "Turkey":
            lines.append(f"Turkey: etiketli WGS={w}; gizli adaylar (ulkesi eksik) dahil ust sinir="
                         f"{w + hid} -> {1000 * (w + hid) / i:.1f} genom/1000 vaka")


def make_zip():
    names = ["summary.txt", "country_summary.csv", "year_summary.csv", "equity_table.csv",
             "turkey_mtbc_runs.tsv", "turkey_mtbc_samples.tsv", "turkey_studies.csv",
             "turkey_hidden_candidates.csv", "who_tb_estimates.csv", "run_log.txt"]
    z = OUT / "mycotb_pilot_results.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in names:
            p = OUT / n
            if p.exists():
                zf.write(p, n)
    return z


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"Cikti klasoru: {OUT}")
    lines = [f"MycoTB pilotu v2 - {datetime.datetime.now():%Y-%m-%d %H:%M}"]
    ncbi_check(lines)
    t_runs, t_samples = download_turkey(lines)
    merged = download_tb(lines)
    wgs_by_country, extra = analyze_global(merged, lines)
    analyze_turkey(t_runs, t_samples, lines)
    who = download_who(lines)
    equity(wgs_by_country, who, extra, lines)
    text = "\n".join(lines)
    (OUT / "summary.txt").write_text(text, encoding="utf-8")
    print("\n" + text + "\n")
    z = make_zip()
    log(f"BITTI. Bu dosyayi Claude'a yukleyin: {z}")


if __name__ == "__main__":
    sys.exit(main())
