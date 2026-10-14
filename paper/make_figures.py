import pandas as pd, numpy as np, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os
BASE=os.path.dirname(os.path.abspath(__file__))
S=os.path.join(BASE,"..","data")+"/"
OUT=os.path.join(BASE,"figures")+"/"
mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.linewidth":.7,"axes.spines.top":False,"axes.spines.right":False,
 "axes.titlesize":9,"axes.titleweight":"bold","axes.titlelocation":"left","legend.frameon":False,"figure.dpi":150,"savefig.dpi":600,"pdf.fonttype":42})
C=dict(bur="#0072B2",gen="#D55E00",gray="#8c8c8c",tr="#CC79A7",green="#009E73",yel="#E69F00")
def lab(ax,t): ax.text(-0.14,1.08,t,transform=ax.transAxes,fontsize=11,fontweight="bold")
def save(f,n):
    f.savefig(OUT+n+".png",bbox_inches="tight"); f.savefig(OUT+n+".pdf",bbox_inches="tight")
    f.savefig(OUT+n+".tiff",bbox_inches="tight",dpi=600,pil_kwargs={"compression":"tiff_lzw"}); plt.close(f)
cr=pd.read_csv(S+"mycotb_robustness_out/country_representation.csv")
lz=pd.read_csv(S+"mycotb_robustness_out/lorenz_main.csv")
ig=pd.read_csv(S+"mycotb_regression_out/income_groups.csv")
ds=pd.read_csv(S+"mycotb_regression_out/analysis_dataset.csv")
ev=pd.read_csv(S+"mycotb_regression_out/expected_vs_observed.csv")
mo=pd.read_csv(S+"mycotb_regression_out/models.csv")
rr=pd.read_csv(S+"mycotb_extended_out/regions_release.csv"); bs=pd.read_csv(S+"mycotb_extended_out/study_bootstrap.csv").set_index("metric"); cap=pd.read_csv(S+"mycotb_extended_out/study_cap_sensitivity.csv"); sd=pd.read_csv(S+"mycotb_extended_out/study_dependence_models.csv")
rc=pd.read_csv(S+"mycotb_robustness_out/release_cohorts.csv"); lm=pd.read_csv(S+"mycotb_robustness_out/lag_matched.csv")
un=pd.read_csv(S+"mycotb_robustness_out/uncertainty.csv"); tn=pd.read_csv(S+"mycotb_analysis_out/turkey_neighbors.csv")

INC=dict(LIC="#D55E00",LMC="#E69F00",UMC="#009E73",HIC="#0072B2")
INCN=dict(LIC="Low",LMC="Lower-middle",UMC="Upper-middle",HIC="High")
d=ds[ds.income.notna()].copy()
# ---------- Fig 1
f,ax=plt.subplots(1,3,figsize=(7.4,2.9),gridspec_kw=dict(width_ratios=[.9,1.35,.95],wspace=.5))
a=ax[0]; a.plot([0,1],[0,1],ls="--",c=C["gray"],lw=.8); a.plot(lz.cum_burden_share,lz.cum_genome_share,c=C["gen"],lw=1.8)
a.fill_between(lz.cum_burden_share,lz.cum_genome_share,lz.cum_burden_share,color=C["gen"],alpha=.12)
a.set_xlabel("Cumulative share of TB burden"); a.set_ylabel("Cumulative share of genomes"); a.set_aspect("equal"); a.set_title("Lorenz curve"); lab(a,"a")
a.text(.04,.9,f"Gini = {bs.loc['gini','point']:.3f}\n(95% CI {bs.loc['gini','ci_lo']:.2f}-{bs.loc['gini','ci_hi']:.2f})",fontsize=6.5,va="top")
a=ax[1]; y=d.obs.clip(lower=.3)
for g in ["LIC","LMC","UMC","HIC"]:
    m=d.income==g; a.scatter(d.cases[m]/1000,y[m],s=9,c=INC[g],alpha=.75,lw=0,label=INCN[g])
xs=np.array([3,6e4])
for r_,l_ in [(.01,"0.01"),(.1,"0.1"),(1,"1"),(10,"10"),(100,"100")]:
    a.plot(xs,xs*r_,c=C["gray"],lw=.5,ls=":"); 
    if r_ in (.1,1): a.text(xs[1]*.9,xs[1]*r_*.9,l_,fontsize=5.5,color=C["gray"],ha="right",va="top")
a.axhline(.3,c=C["gray"],lw=.4); a.text(5e4,.36,"zero genomes",fontsize=5.5,color=C["gray"],ha="right")
for n in ["India","China","South Africa","Indonesia","Philippines","Pakistan","Nigeria","United States of America","Georgia","Turkey","Viet Nam"]:
    r=d[d.country==n]
    if len(r): r=r.iloc[0]; a.annotate(n.replace("United States of America","USA"),(r.cases/1000,max(r.obs,.3)),xytext=(3,3),textcoords="offset points",fontsize=5.5)
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(3,6e4); a.set_ylim(.2,1e5)
a.set_xlabel("Incident cases 2010-2024 (thousands)"); a.set_ylabel("Public genomes"); a.set_title("Countries (n = 207)"); lab(a,"b")
a.legend(fontsize=5.8,loc="upper center",bbox_to_anchor=(.5,-.2),ncol=4,handletextpad=.1,columnspacing=.7,markerscale=1.2);
a=ax[2]; keys=["LIC","LMC","UMC","HIC"]; pt=[bs.loc[f"RI_{k}","point"] for k in keys]; lo=[bs.loc[f"RI_{k}","ci_lo"] for k in keys]; hi=[bs.loc[f"RI_{k}","ci_hi"] for k in keys]
for i,k in enumerate(keys): a.errorbar(i,pt[i],yerr=[[pt[i]-lo[i]],[hi[i]-pt[i]]],fmt="o",c=INC[k],capsize=2.5,ms=5,lw=1.2)
a.axhline(1,c="k",lw=.6,ls=":"); a.set_yscale("log"); a.set_yticks([.1,.3,1,3,10,30]); a.set_yticklabels(["0.1","0.3","1","3","10","30"]); a.set_ylim(.05,45); a.minorticks_off()
a.set_xticks(range(4)); a.set_xticklabels(keys); a.set_xlim(-.5,3.5); a.set_ylabel("Representation index (log)"); a.set_title("Income group"); lab(a,"c")
a.text(.02,.97,"RI = genome share /\nburden share;\nbars: study-level\nbootstrap 95% CI",transform=a.transAxes,fontsize=5.5,color=C["gray"],va="top")
save(f,"Fig1_concentration")
# ---------- Fig 2
cr["category"]=cr.category.astype(str)
t=cr.groupby("category").agg(n=("country","count"),b=("burden_share_pct","sum"),g=("genome_share_pct","sum")).reindex(["0","<1","1-10",">=10"]).fillna(0)
f,ax=plt.subplots(1,3,figsize=(7.4,3.0),gridspec_kw=dict(width_ratios=[.8,1.0,1.25],wspace=.7))
a=ax[0]; cols=["#999999","#D55E00","#E69F00","#0072B2"]; bot=np.zeros(2)
for c_,col in zip(t.index,cols):
    v=np.array([t.loc[c_,"b"],t.loc[c_,"g"]]); a.bar([0,1],v,bottom=bot,color=col,width=.6,label=f"{c_} (n={int(t.loc[c_,'n'])})"); bot+=v
a.set_xticks([0,1]); a.set_xticklabels(["Burden","Genomes"]); a.set_ylabel("Share (%)"); a.legend(title="Genomes / 1,000 cases",fontsize=6,title_fontsize=6.5,loc="upper center",bbox_to_anchor=(.5,-.12),ncol=2)
a.set_title("Coverage category"); lab(a,"a")
a=ax[1]; g=cr.groupby("who_region").agg(b=("burden_share_pct","sum"),g=("genome_share_pct","sum")).sort_values("b"); y=np.arange(len(g))
a.barh(y+.2,g.b,.38,color=C["bur"],label="TB burden"); a.barh(y-.2,g.g,.38,color=C["gen"],label="Genomes")
a.set_yticks(y); a.set_yticklabels(g.index); a.set_xlabel("Share (%)"); a.set_title("WHO region"); lab(a,"b"); a.legend(fontsize=6,loc="upper center",bbox_to_anchor=(.5,-.2),ncol=2)
a=ax[2]; top=cr.sort_values("who_incidence_2010_2024",ascending=False).head(15).iloc[::-1]
v=top.genomes_per_1000.clip(lower=.003); a.barh(top.country.str.replace("Democratic Republic of the Congo","DR Congo"),v,color=[C["bur"] if x>=1 else C["gen"] for x in top.genomes_per_1000])
a.set_xscale("log"); a.axvline(1,c="k",lw=.6,ls=":"); a.set_xlabel("Genomes per 1,000 cases (log)"); a.set_title("15 highest-burden countries"); lab(a,"c")
save(f,"Fig2_coverage")
# ---------- Fig 3 (2x2)
f,ax=plt.subplots(2,2,figsize=(6.6,5.2),gridspec_kw=dict(wspace=.38,hspace=.55)); ax=ax.ravel()
a=ax[0]; xs=np.arange(len(un)); a.errorbar(xs,un.gini_point,yerr=[un.gini_point-un.gini_lo,un.gini_hi-un.gini_point],fmt="o-",c=C["gray"],capsize=2)
a.set_xticks(xs); a.set_xticklabels(un.period); a.set_ylim(.8,1); a.set_ylabel("Gini coefficient"); a.set_title("Naive, by collection period"); lab(a,"a")
a.text(.5,.06,"confounded by publication lag",transform=a.transAxes,ha="center",fontsize=6.5,color=C["gray"])
a=ax[1]; a.plot(np.arange(len(lm)),lm.gini,"s-",c=C["bur"]); a.set_xticks(np.arange(len(lm))); a.set_xticklabels(lm.collection_period); a.set_ylim(.8,1); a.set_ylabel("Gini coefficient"); a.set_title("Lag-matched (lag <= 3 y)"); lab(a,"b")
a=ax[2]; xs=np.arange(len(rc)); a.bar(xs,rc.genomes_per_1000,color=C["green"]); a.set_xticks(xs); a.set_xticklabels([w_.replace("20","") for w_ in rc.release_window],fontsize=7)
a.set_ylabel("Genomes per 1,000 cases"); a.set_xlabel("Release window (20xx)"); a.set_title("By release window"); lab(a,"c")
for i,r in rc.iterrows(): a.text(i,r.genomes_per_1000+.05,f"{int(r.zero_genome_countries)}",ha="center",fontsize=6.5)
a.text(.02,.95,"numbers: countries with\nzero genomes",transform=a.transAxes,fontsize=6,color=C["gray"],va="top")
a=ax[3]; rcol=dict(AFR=C["yel"],AMR=C["green"],EMR=C["tr"],EUR=C["bur"],SEA=C["gen"],WPR="#444444")
for reg,gdf in rr.groupby("who_region"):
    gdf=gdf.sort_values("release_window"); yv=gdf.representation_index.clip(lower=.03); a.plot(range(len(gdf)),yv,"o-",c=rcol[reg],ms=3,lw=1.2); a.text(3.08,yv.iloc[-1]*{"AFR":1.18,"SEA":.85}.get(reg,1),reg,fontsize=6.5,color=rcol[reg],va="center")
a.axhline(1,c="k",lw=.6,ls=":"); a.set_yscale("log"); a.set_yticks([.03,.1,.3,1,3,10]); a.set_yticklabels(["0","0.1","0.3","1","3","10"]); a.minorticks_off(); a.set_xlim(-.2,3.6)
a.set_xticks(range(4)); a.set_xticklabels(sorted(rr.release_window.unique()),fontsize=7); a.set_ylabel("Representation index (log)"); a.set_title("WHO region, by release window"); lab(a,"d")
save(f,"Fig3_temporal")
# ---------- Fig 4 forest
want=[("P3 ANA","Main model (P3)"),("P4","+ case volume (P4)"),("P5","+ WHO region (P5)"),("S1","Excluding USA"),("S2","Low/middle income only"),("P6","Health expenditure (P6)")]
terms=[("log2gdp","GDP per capita"),("log2rate","TB incidence rate"),("log2artpc","Articles per million")]
f,ax=plt.subplots(1,3,figsize=(7.2,2.6),sharey=True,gridspec_kw=dict(wspace=.08))
for a,(tm,ttl) in zip(ax,terms):
    for i,(k,nm) in enumerate(want):
        r=mo[mo.model.str.startswith(k)&(mo.term==tm)|(mo.model.str.startswith(k)&(k=="P6")&(mo.term==tm.replace("gdp","health")))]
        if r.empty: continue
        r=r.iloc[0]; a.errorbar(r.ratio_exp_coef,-i,xerr=[[r.ratio_exp_coef-r.ci_lo],[r.ci_hi-r.ratio_exp_coef]],fmt="o",c=C["bur"] if r.ci_lo>1 or r.ci_hi<1 else C["gray"],ms=4,lw=1,capsize=2)
    a.axvline(1,c="k",lw=.6,ls=":"); a.set_xscale("log"); a.set_xticks([.25,.5,1,2,4,8]); a.set_xticklabels(["0.25","0.5","1","2","4","8"]); a.minorticks_off(); a.set_xlim(.3,8); a.set_title(ttl,fontsize=7,loc="center"); a.set_xlabel("Rate ratio per doubling (95% CI)")
ax[0].set_yticks([-i for i in range(len(want))]); ax[0].set_yticklabels([w[1] for w in want])
save(f,"Fig4_regression")
# ---------- Fig 5
f,ax=plt.subplots(1,3,figsize=(7.4,2.9),gridspec_kw=dict(width_ratios=[1.1,1,1],wspace=.65))
a=ax[0]; e=ev[ev.observed_genomes>0]; a.scatter(e.expected_genomes_P3,e.observed_genomes,s=8,c=C["gray"],alpha=.6,lw=0)
m=max(ev.expected_genomes_P3.max(),ev.observed_genomes.max()); a.plot([1,m],[1,m],"--",c="k",lw=.7)
t_=ev[ev.country=="Turkey"].iloc[0]; a.scatter([t_.expected_genomes_P3],[.7],s=40,c=C["tr"],zorder=5,marker="v"); a.annotate("Turkey\n(0 observed)",(t_.expected_genomes_P3,.7),xytext=(-62,16),textcoords="offset points",arrowprops=dict(arrowstyle="-",color=C["tr"],lw=.6),fontsize=7,color=C["tr"],fontweight="bold")
a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(.5,None); a.set_xlabel("Expected genomes (model P3)"); a.set_ylabel("Observed genomes"); a.set_title("Observed vs expected"); lab(a,"a")
a=ax[1]; s=ev.sort_values("expected_minus_observed_P3",ascending=False).head(12).iloc[::-1]
a.barh(s.country,s.expected_minus_observed_P3/1000,color=[C["tr"] if c=="Turkey" else C["gen"] for c in s.country]); a.set_xlabel("Shortfall (thousand genomes)"); a.set_title("Largest shortfalls"); lab(a,"b")
a=ax[2]; n=tn.iloc[::-1].reset_index(drop=True); nz=n.genomes_per_1000_cases>0
a.barh(n.index[nz],n.genomes_per_1000_cases[nz],color=[C["tr"] if c=="Turkey" else C["bur"] for c in n.country[nz]]); a.set_xscale("log"); a.set_xlim(.004,400)
a.scatter([.006]*int((~nz).sum()),n.index[~nz],marker="x",s=14,c=[C["tr"] if c=="Turkey" else C["gray"] for c in n.country[~nz]],lw=1)
a.set_yticks(n.index); a.set_yticklabels(n.country); a.set_xlabel("Genomes per 1,000 cases (log)"); a.set_title("Turkey and neighbours"); a.text(.97,.02,"x = zero genomes",transform=a.transAxes,fontsize=6,color=C["gray"],ha="right"); lab(a,"c")
save(f,"Fig5_turkey")
# ---------- Fig 6 study dependence
f,ax=plt.subplots(1,3,figsize=(7.4,2.8),gridspec_kw=dict(width_ratios=[1.15,.9,1.2],wspace=.55))
a=ax[0]; lbl=["No cap","Cap 100\nper study","Cap 25\nper study"]; w=.2
for j,k in enumerate(["LIC","LMC","UMC","HIC"]):
    a.bar(np.arange(3)+(j-1.5)*w,cap[f"RI_{k}"],w,color=INC[k],label=k)
a.axhline(1,c="k",lw=.6,ls=":"); a.set_yscale("log"); a.set_ylim(.07,40); a.set_yticks([.1,.3,1,3,10,30]); a.set_yticklabels(["0.1","0.3","1","3","10","30"]); a.minorticks_off()
a.set_xticks(range(3)); a.set_xticklabels(lbl,fontsize=6.5); a.set_ylabel("Representation index (log)"); a.legend(ncol=4,fontsize=5.5,loc="upper center",bbox_to_anchor=(.5,-.28),columnspacing=.6,handlelength=.8); a.set_title("Cap per study"); lab(a,"a")
a=ax[1]; mt=["gini","top30_genome_pct"]
a.errorbar([bs.loc["gini","point"]],[0],xerr=[[bs.loc["gini","point"]-bs.loc["gini","ci_lo"]],[bs.loc["gini","ci_hi"]-bs.loc["gini","point"]]],fmt="o",c=C["gen"],capsize=2.5)
for i,(k,nm) in enumerate([("no cap","No cap"),("cap 100 per study","Cap 100"),("cap 25 per study","Cap 25")]):
    gv=cap[cap.scenario==k].gini.iloc[0]; a.scatter([gv],[-(i+1)*.6],c=C["bur"],s=14,zorder=3)
a.set_yticks([0,-.6,-1.2,-1.8]); a.set_yticklabels(["Bootstrap\n(studies)","No cap","Cap 100","Cap 25"],fontsize=6.5); a.set_ylim(-2.3,.5); a.set_xlim(.74,.95); a.set_xlabel("Gini coefficient"); a.set_title("Gini"); lab(a,"b")
a=ax[2]; sel=[("P3 HC1 (published)","Robust SE (P3)"),("P3 region-clustered CR1 (t, G-1 df)","Region-clustered SE"),("Outcome: genomes capped at 25/study, HC1","Capped at 25/study"),("Outcome: number of distinct studies, HC1","Number of studies")]
for i,(k,nm) in enumerate(sel):
    for tm,col,off_ in [("log2 GDP",C["bur"],.12),("log2 rate",C["gen"],-.12)]:
        r=sd[(sd.model==k)&(sd.term==tm)].iloc[0]; a.errorbar(r.rr,-i+off_,xerr=[[r.rr-r.ci_lo],[r.ci_hi-r.rr]],fmt="o",c=col,ms=3.5,lw=1,capsize=1.8,label=("GDP per capita" if tm=="log2 GDP" else "Incidence rate") if i==0 else None)
a.axvline(1,c="k",lw=.6,ls=":"); a.set_xscale("log"); a.set_xticks([.3,.5,1,2,4,8]); a.set_xticklabels(["0.3","0.5","1","2","4","8"]); a.minorticks_off(); a.set_xlim(.28,9)
a.set_yticks([-i for i in range(4)]); a.set_yticklabels([s_[1] for s_ in sel],fontsize=6.5); a.set_xlabel("Rate ratio per doubling (95% CI)"); a.legend(fontsize=5.8,loc="upper center",bbox_to_anchor=(.5,-.3),ncol=2); a.set_title("Model specification"); lab(a,"c")
save(f,"Fig6_study_dependence")
print("ok")
