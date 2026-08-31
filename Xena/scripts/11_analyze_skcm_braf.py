"""Analyze real SKCM-BRAF patient table and emit aggregate result files."""
from pathlib import Path
import json, math
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; P=ROOT/"processed"; df=pd.read_csv(P/"skcm_braf_patient_table.tsv",sep="\t")
def perm(a,b,n=10000):
    import numpy as np
    rng=np.random.default_rng(20260830); z=np.r_[a,b]; k=len(a); obs=abs(np.mean(a)-np.mean(b)); hits=0
    for _ in range(n):
        rng.shuffle(z); hits += abs(np.mean(z[:k])-np.mean(z[k:]))>=obs
    return (hits+1)/(n+1)
def compare(col,out):
    x=df[df.braf_v600_hotspot.eq(1)][col].dropna().to_numpy(); y=df[df.braf_v600_hotspot.eq(0)][col].dropna().to_numpy(); pooled=((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2)
    rows={"braf_mutated_n":len(x),"braf_wild_type_n":len(y),"braf_mutated_mean":float(x.mean()),"braf_wild_type_mean":float(y.mean()),"mean_difference_mutated_minus_wild_type":float(x.mean()-y.mean()),"cohens_d":float((x.mean()-y.mean())/math.sqrt(pooled)),"two_sided_permutation_p_value":float(perm(x,y))}
    pd.DataFrame(rows.items(),columns=["metric","value"]).to_csv(P/out,sep="\t",index=False); return rows
expression=compare("braf_expression","braf_expression_comparison.tsv")
groups=[]
for g,label in [(1,"BRAF hotspot"),(0,"BRAF wild type")]:
    q=df[df.braf_v600_hotspot.eq(g)]; groups.append({"group":label,"n":len(q),"female_n":int(q.gender.astype(str).str.upper().eq("FEMALE").sum()),"female_percent":float(100*q.gender.astype(str).str.upper().eq("FEMALE").mean()),"median_age":float(pd.to_numeric(q.age_at_initial_pathologic_diagnosis,errors="coerce").median()),"primary_tumour_n":int(q.sample_type.eq("Primary Tumor").sum())})
json.dump({"groups":groups,"note":"Descriptive, unadjusted comparisons"},open(P/"skcm_clinical_characteristics.json","w"),indent=2)
def cox(rows):
    beta=0.; times=sorted(set(rows.os_days[rows.os_event.eq(1)]));
    for _ in range(40):
        score=info=0.
        for t in times:
            risk=rows[rows.os_days>=t]; ev=rows[(rows.os_days==t)&rows.os_event.eq(1)]; w=(beta*risk.braf_v600_hotspot).map(math.exp); s0=w.sum(); s1=(w*risk.braf_v600_hotspot).sum(); s2=(w*risk.braf_v600_hotspot**2).sum(); m=len(ev)
            score += ev.braf_v600_hotspot.sum()-m*s1/s0; info += m*(s2/s0-(s1/s0)**2)
        step=score/info; beta+=step
        if abs(step)<1e-10: break
    se=1/math.sqrt(info); hr=math.exp(beta); return {"analysis_n":len(rows),"deaths":int(rows.os_event.sum()),"hazard_ratio":hr,"95ci_lower":math.exp(beta-1.96*se),"95ci_upper":math.exp(beta+1.96*se),"p_value":math.erfc(abs(beta/se)/math.sqrt(2))}
surv=df.dropna(subset=["os_days","os_event"]); mut=surv[surv.braf_v600_hotspot.eq(1)]; wt=surv[surv.braf_v600_hotspot.eq(0)]
surv_result={"braf_mutated_n":len(mut),"braf_wild_type_n":len(wt),"braf_mutated_deaths":int(mut.os_event.sum()),"braf_wild_type_deaths":int(wt.os_event.sum()),**cox(surv)}
json.dump(surv_result,open(P/"braf_mutation_survival.json","w"),indent=2)
expr_surv=df.dropna(subset=["braf_expression","os_days","os_event"]).copy(); expr_surv["braf_v600_hotspot"]=pd.qcut(expr_surv.braf_expression,2,labels=[0,1]).astype(int); ex=cox(expr_surv); ex["expression_median_split_n"]=len(expr_surv); json.dump(ex,open(P/"braf_expression_survival.json","w"),indent=2)
summary={"site":"Xena_TCGA_SKCM","analysis_type":"local_summary_statistics","cohort":{"clinical_patients":len(df),"mutation_evaluable_patients":int(df.mutation_evaluable.sum()),"braf_v600_hotspot_patients":int(df.braf_v600_hotspot.sum()),"braf_wild_type_patients":int((df.braf_v600_hotspot==0).sum()),"primary_expression_patients":int(df.expression_available.sum())},"braf_expression_by_mutation_status":expression,"clinical_characteristics_by_mutation_status":groups,"braf_mutation_overall_survival":surv_result,"braf_expression_overall_survival":ex,"analysis_definitions":{"hotspot":"GDC transcript annotation p.V640E or p.V640G, corresponding to the canonical BRAF V600 activating hotspot","expression":"BRAF row from Xena TCGA.SKCM.sampleMap/HiSeqV2; primary-tumour sample when available","survival":"Dead=event; living patients use latest follow-up"},"source_files":["SKCM/HiSeqV2.gz","SKCM/SKCM_clinicalMatrix","SKCM/TCGA-SKCM.somaticmutation_wxs.tsv.gz"],"limitations":["Observational TCGA data; survival models are unadjusted","Only 104 patients have primary-tumour expression in the legacy Xena matrix","This is aggregate-only output; patient rows are not included"]}
json.dump(summary,open(P/"skcm_braf_site_summary.json","w"),indent=2); print(json.dumps({"expression":expression,"survival":surv_result,"cohort":summary["cohort"]},indent=2))

