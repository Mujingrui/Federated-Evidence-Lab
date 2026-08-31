"""Build real Xena SKCM-BRAF patient table."""
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; RAW=ROOT/"raw"/"SKCM"; OUT=ROOT/"processed"/"skcm_braf_patient_table.tsv"
cl=pd.read_csv(RAW/"SKCM_clinicalMatrix",sep="\t",low_memory=False).dropna(subset=["bcr_patient_barcode"]).copy()
cl["patient_id"]=cl["bcr_patient_barcode"].str[:12]
cl=cl.sort_values(["patient_id","sample_type_id"],na_position="last").drop_duplicates("patient_id")
mut=pd.read_csv(RAW/"TCGA-SKCM.somaticmutation_wxs.tsv.gz",sep="\t",compression="gzip",usecols=["Tumor_Sample_Barcode","gene","Amino_Acid_Change"])
mut["patient_id"]=mut["Tumor_Sample_Barcode"].str[:12]
# GDC's transcript annotation reports canonical V600E/G as V640E/G.
hot=set(mut.loc[(mut.gene=="BRAF")&mut.Amino_Acid_Change.isin(["p.V640E","p.V640G"]),"patient_id"])
expr=pd.read_csv(RAW/"HiSeqV2.gz",sep="\t",compression="gzip")
braf=expr.loc[expr.iloc[:,0].eq("BRAF")].iloc[0]
primary={c[:12]:braf[c] for c in expr.columns[1:] if "-01" in c[12:15]}
cl["mutation_evaluable"]=1; cl["braf_v600_hotspot"]=cl.patient_id.isin(hot).astype(int)
cl["braf_expression"]=cl.patient_id.map(primary); cl["expression_available"]=cl.braf_expression.notna().astype(int)
cl["os_event"]=cl.vital_status.map({"DECEASED":1,"LIVING":0})
cl["os_days"]=pd.to_numeric(cl.get("days_to_death"),errors="coerce")
follow=pd.to_numeric(cl.get("days_to_last_followup"),errors="coerce"); cl.loc[cl.os_event.eq(0),"os_days"]=follow[cl.os_event.eq(0)]
keep=["patient_id","bcr_patient_barcode","sample_type","mutation_evaluable","braf_v600_hotspot","braf_expression","expression_available","vital_status","os_event","os_days","age_at_initial_pathologic_diagnosis","gender","pathologic_stage"]
cl[keep].to_csv(OUT,sep="\t",index=False)
print("Output:",OUT); print("Patients:",len(cl)); print("BRAF V600 hotspot:",int(cl.braf_v600_hotspot.sum())); print("Primary expression:",int(cl.expression_available.sum()))

