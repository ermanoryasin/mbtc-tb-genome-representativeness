#!/usr/bin/env python3
"""Builds Supplementary_Tables.xlsx with sheets numbered by order of first citation in the manuscript (supp_map.json) and readable headers."""
import json, re
import pandas as pd
from pathlib import Path
from openpyxl import Workbook
HERE = Path(__file__).resolve().parent
D = HERE.parent / "data" / "mycotb_extended_out"; T = HERE / "tables"
SMAP = {int(k): v for k, v in json.load(open(HERE / "supp_map.json")).items()}
FIX = {"RI_LIC": "RI low", "RI_LMC": "RI lower-middle", "RI_UMC": "RI upper-middle", "RI_HIC": "RI high", "gini": "Gini coefficient", "zero_countries": "Countries with no genome",
 "zero_genome_countries": "Countries with no genome", "genomes_per_1000": "Genomes per 1000 cases", "top30_genome_pct": "Top-30 genome share (%)", "rr": "Rate ratio", "ci_lo": "CI lower", "ci_hi": "CI upper",
 "p": "p", "n": "n", "who_region": "WHO region", "release_window": "Release window", "burden_share_pct": "Burden share (%)", "genome_share_pct": "Genome share (%)", "representation_index": "RI",
 "RI_year_specific": "RI, year-specific classification", "RI_current_classification": "RI, current classification", "burden_share_year_specific_pct": "Burden share, year-specific (%)",
 "genome_share_year_specific_pct": "Genome share, year-specific (%)", "lag_threshold_years": "Lag threshold (years)", "p_obs_le_under_nb2": "P(count <= observed), NB2", "p_obs_le_averaged": "P(count <= observed), averaged over draws",
 "share_without_country_pct": "Share without country (%)", "study_income_group": "Income group of study", "without_country": "Samples without country", "depositor_class": "Depositor class",
 "step": "Step", "count": "Count", "reason_or_note": "Reason or note", "expected_poisson": "Expected, Poisson", "expected_nb2": "Expected, NB2", "expected_nb2_lo": "NB2 interval lower", "expected_nb2_hi": "NB2 interval upper",
 "shortfall_rank_nb2": "Shortfall rank, NB2", "shortfall_rank_poisson": "Shortfall rank, Poisson", "observed": "Observed genomes", "expected_poisson_P3": "Expected, Poisson P3", "expected_nb2_P3": "Expected, NB2 P3"}
for g in ("LIC", "LMC", "UMC", "HIC"): FIX[f"share_{g}_pct"] = f"Share of genomes, {dict(LIC='low',LMC='lower-middle',UMC='upper-middle',HIC='high')[g]} income (%)"
def pretty(c): return FIX.get(c) or c.replace("_", " ").capitalize()
def load(p): df = pd.read_csv(p); df.columns = [pretty(c) for c in df.columns]; return df
m = pd.read_csv(T / "Table3_models.csv")
mm = {'P1 gelir grubu (ref=LIC)': 'P1 income groups (ref = low income)', 'P2 GSYH + insidans hizi': 'P2 GDP + incidence rate', 'P3 ANA MODEL: GSYH + insidans hizi + makale/milyon': 'P3 main model: GDP + incidence rate + articles per million', 'P4 P3 + log2 kumulatif vaka (olcek)': 'P4 P3 + log2 cumulative cases (case volume)', 'P5 P3 + WHO bolgesi (ref=AFR)': 'P5 P3 + WHO region (ref = Africa)', 'P6 saglik harcamasi ile (GSYH yerine)': 'P6 health expenditure instead of GDP', 'S1 ABD haric': 'S1 excluding the USA', 'S2 yalniz dusuk/orta gelirli': 'S2 low- and middle-income countries only', 'L1 >=1 genom (vaka olcegi dahil)': 'L1 logistic: at least one genome (with case volume)'}
tm = {'inc_LMC': 'income lower-middle', 'inc_UMC': 'income upper-middle', 'inc_HIC': 'income high', 'log2gdp': 'log2 GDP per capita', 'log2rate': 'log2 incidence rate', 'log2artpc': 'log2 articles per million', 'log2cases': 'log2 cumulative cases', 'reg_AMR': 'region Americas', 'reg_EMR': 'region Eastern Mediterranean', 'reg_EUR': 'region Europe', 'reg_SEA': 'region South-East Asia', 'reg_WPR': 'region Western Pacific', 'log2health': 'log2 health expenditure per capita'}
m["model"] = m.model.map(lambda x: mm.get(x, x)); m["term"] = m.term.map(lambda x: tm.get(x, x))
m = m.rename(columns={"ratio_exp_coef": "Rate ratio (or odds ratio)"}); m.columns = [c if c == "Rate ratio (or odds ratio)" else pretty(c) for c in m.columns]
ls = load(D / "literature_search_log.csv")
# (old number, suffix, descriptor, dataframe)
SHEETS = [(1, "", "study cap", load(D / "study_cap_sensitivity.csv")), (2, "", "region by release", load(D / "regions_release.csv")),
 (3, "", "model specification", load(D / "study_dependence_models.csv")), (4, "", "study bootstrap", load(D / "study_bootstrap.csv")),
 (5, "", "regression models", m), (6, "", "Turkey provenance", load(D / "turkey_provenance.csv")),
 (7, "", "NB and Poisson models", load(D / "negbin_models.csv")), (8, "", "Turkey NB expectations", load(D / "negbin_turkey.csv")),
 (9, "", "expected counts P3", load(D / "negbin_expected_P3.csv")), (10, "", "analysis set flow", load(D / "flow_table.csv")),
 (11, "", "missing scenarios", load(D / "missing_country_scenarios.csv")), (11, "b", "missing probe", load(D / "missing_country_probe.csv")),
 (12, "", "lag thresholds", load(D / "lag_threshold_windows.csv")), (12, "b", "fast vs slow", load(D / "lag_fast_slow_composition.csv")),
 (12, "c", "annual lag series", load(D / "lag_annual_k3.csv")), (13, "", "income by year", load(D / "income_year_specific.csv")), (14, "", "literature search", ls)]
SHEETS.sort(key=lambda s: (SMAP[s[0]], s[1]))
wb = Workbook(); wb.remove(wb.active); index = []
for old, suf, desc, df in SHEETS:
    name = f"S{SMAP[old]}{suf} {desc}"[:31]; ws = wb.create_sheet(name); ws.append(list(df.columns))
    for r in df.itertuples(index=False): ws.append([None if (isinstance(v, float) and pd.isna(v)) else v for v in r])
    for col in ws.columns: ws.column_dimensions[col[0].column_letter].width = min(60, max(10, max(len(str(c.value or "")) for c in col) + 2))
    index.append((name, desc))
wb.save(HERE / "submission_EI" / "Supplementary_Tables.xlsx"); json.dump(index, open(HERE / "supp_index.json", "w"))
print([n for n, _ in index])
