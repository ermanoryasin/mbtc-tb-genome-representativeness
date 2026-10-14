"""Exports manuscript/supplementary tables from the analysis outputs to paper/tables/."""
import os, pandas as pd
BASE=os.path.dirname(os.path.abspath(__file__)); D=os.path.join(BASE,"..","data")+"/"; T=os.path.join(BASE,"tables")+"/"
os.makedirs(T,exist_ok=True)
pd.read_csv(D+"mycotb_robustness_out/country_representation.csv").head(20).to_csv(T+"Table1_top20_burden.csv",index=False)
pd.read_csv(D+"mycotb_regression_out/income_groups.csv").to_csv(T+"Table2_income_groups.csv",index=False)
pd.read_csv(D+"mycotb_regression_out/models.csv").to_csv(T+"Table3_models.csv",index=False)
for n in ["periods","scenarios"]: pd.read_csv(D+f"mycotb_analysis_out/{n}.csv").to_csv(T+f"TableS_{n}.csv",index=False)
for n in ["sensitivity","uncertainty","release_cohorts","lag_matched"]: pd.read_csv(D+f"mycotb_robustness_out/{n}.csv").to_csv(T+f"TableS_{n}.csv",index=False)
for n in ["regions_release","study_cap_sensitivity","study_bootstrap","study_dependence_models"]: pd.read_csv(D+f"mycotb_extended_out/{n}.csv").to_csv(T+f"TableS_{n}.csv",index=False)
print("tables ok")
