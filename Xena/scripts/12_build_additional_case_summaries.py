"""Build aggregate patient summaries for BRCA-ERBB2, COAD-BRAF and GBM-EGFR."""
from pathlib import Path
import json, math
import numpy as np, pandas as pd
ROOT=Path(__file__).resolve().parents[1];
CASES={"BRCA":{"gene":"ERBB2","hotspots":None},"COAD":{"gene":"BRAF","hotspots":["p.V640E"]},"GBM":{"gene":"EGFR","hotspots":None}}
def cox(d):
 d=d.dropna(subset=["os_days","os_event"]); beta=0.; times=sorted(d.loc[d.os_event.eq(1),"os_days"].unique())
 for _ in range(40):
  score=info=0.
  for t in times:
   r=d[d.os_days>=t]; e=d[(d.os_days==t)&d.os_event.eq(1)]; w=np.exp(beta*r.status); s0=w.sum(); s1=(w*r.status).sum(); s2=(w*r.status**2).sum(); m=len(e); score+=e.status.sum()-m*s1/s0; info+=m*(s2/s0-(s1/s0)**2)
  step=score/info; beta+=step
  if abs(step)<1e-10: break
 se=1/math.sqrt(info); return {"analysis_n":len(d),"deaths":int(d.os_event.sum()),"hazard_ratio":math.exp(beta),"95ci_lower":math.exp(beta-1.96*se),"95ci_upper":math.exp(beta+1.96*se),"p_value":math.erfc(abs(beta/se)/math.sqrt(2))}
for cohort,meta in CASES.items():
 raw=ROOT/"raw"/cohort; out=ROOT/"processed"; cl=pd.read_csv(raw/f"{cohort}_clinicalMatrix",sep="\t",low_memory=False).dropna(subset=["bcr_patient_barcode"]).copy(); cl["patient_id"]=cl.bcr_patient_barcode.str[:12]; cl=cl.sort_values(["patient_id","sample_type_id"],na_position="last").drop_duplicates("patient_id")
 m=pd.read_csv(raw/f"TCGA-{cohort}.somaticmutation_wxs.tsv.gz",sep="\t",compression="gzip",usecols=["Tumor_Sample_Barcode","gene","Amino_Acid_Change"]); m["patient_id"]=m.Tumor_Sample_Barcode.str[:12]; x=m[m.gene.eq(meta["gene"])]
 if meta["hotspots"]: pos=set(x.loc[x.Amino_Acid_Change.isin(meta["hotspots"]),"patient_id"])
 else: pos=set(x.loc[~x.Amino_Acid_Change.astype(str).str.contains("=|splice|fs|\*",regex=True),"patient_id"])
 expr=pd.read_csv(raw/"HiSeqV2.gz",sep="\t",compression="gzip"); row=expr.loc[expr.iloc[:,0].eq(meta["gene"])]; values=row.iloc[0].to_dict() if len(row) else {}; primary={c[:12]:values[c] for c in expr.columns[1:] if "-01" in c[12:15]}
 cl["status"]=cl.patient_id.isin(pos).astype(int); cl["expression"]=cl.patient_id.map(primary); cl["expression_available"]=cl.expression.notna().astype(int); cl["os_event"]=cl.vital_status.map({"DECEASED":1,"LIVING":0}); cl["os_days"]=pd.to_numeric(cl.get("days_to_death"),errors="coerce"); fu=pd.to_numeric(cl.get("days_to_last_followup"),errors="coerce"); cl.loc[cl.os_event.eq(0),"os_days"]=fu[cl.os_event.eq(0)]
 ex=cl.dropna(subset=["expression"]); a=ex[ex.status.eq(1)].expression.to_numpy(); b=ex[ex.status.eq(0)].expression.to_numpy(); diff=float(a.mean()-b.mean()) if len(a)>1 and len(b)>1 else None; pooled=((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2) if len(a)>1 and len(b)>1 else None
 clinical=[]
 for g,label in [(1,f"{meta['gene']} altered"),(0,f"{meta['gene']} wild type")]:
  q=cl[cl.status.eq(g)]; clinical.append({"group":label,"n":len(q),"female_n":int(q.gender.astype(str).str.upper().eq("FEMALE").sum()),"female_percent":float(100*q.gender.astype(str).str.upper().eq("FEMALE").mean()),"median_age":float(pd.to_numeric(q.age_at_initial_pathologic_diagnosis,errors="coerce").median())})
 surv=cl.dropna(subset=["os_days","os_event"]); sd=surv.copy(); sd["status"]=sd.status.astype(float); sr=cox(sd); summary={"site":f"Xena_TCGA_{cohort}","analysis_type":"local_summary_statistics","cohort":{"clinical_patients":len(cl),"mutation_evaluable_patients":len(cl),"altered_patients":int(cl.status.sum()),"wild_type_patients":int((cl.status==0).sum()),"primary_expression_patients":int(cl.expression_available.sum())},"expression_by_mutation_status":{"altered_n":len(a),"wild_type_n":len(b),"altered_mean":float(a.mean()) if len(a) else None,"wild_type_mean":float(b.mean()) if len(b) else None,"mean_difference":diff,"cohens_d":float(diff/math.sqrt(pooled)) if pooled and diff is not None else None},"clinical_characteristics_by_mutation_status":clinical,"mutation_overall_survival":sr,"expression_overall_survival":{"note":"Exploratory median split among patients with primary expression","analysis_n":len(ex),**cox(ex.assign(status=(ex.expression>=ex.expression.median()).astype(int)))},"analysis_definitions":{"gene":meta["gene"],"alteration":"COAD BRAF uses p.V640E (canonical V600E transcript equivalent); BRCA ERBB2 and GBM EGFR use nonsynonymous coding mutations excluding synonymous, splice and truncating annotations","expression":"Gene-level HiSeqV2 expression; primary tumour sample when available","survival":"Dead=event; living patients use latest follow-up"},"source_files":[f"{cohort}/HiSeqV2.gz",f"{cohort}/{cohort}_clinicalMatrix",f"{cohort}/TCGA-{cohort}.somaticmutation_wxs.tsv.gz"],"limitations":["Observational TCGA data; survival models are unadjusted","Altered-group definitions are exploratory except the COAD BRAF V600E hotspot","Aggregate output only; no patient-level records"]}
 json.dump(summary,open(out/f"{cohort.lower()}_{meta['gene'].lower()}_site_summary.json","w"),indent=2); print(cohort,summary["cohort"],summary["expression_by_mutation_status"])

