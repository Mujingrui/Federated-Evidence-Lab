from pathlib import Path
import json

ROOT=Path('/Users/jingruimu/Documents/ChatGPT/federated AI')
sites=[('SKCM','BRAF','skcm_braf','skcm_braf'),('BRCA','ERBB2','brca_erbb2','brca_erbb2'),('COAD','BRAF','coad_braf','coad_braf'),('GBM','EGFR','gbm_egfr','gb_egfr')]
for disease,gene,xkey,dkey in sites:
    xf=ROOT/'Xena/processed'/f'{xkey}_site_summary.json'
    df=ROOT/'DepMap/processed'/f'{dkey}_drug_summary.json'
    x=json.loads(xf.read_text()) if xf.exists() else {'status':'not_available'}
    d=json.loads(df.read_text()) if df.exists() else {'status':'not_available'}
    out={'site':f'{disease}-{gene}','patient_site':x,'cell_drug_site':d,'federated_interpretation':{'patient_evidence':'Xena provides disease-cohort mutation, expression, clinical, and survival aggregates.','cell_drug_evidence':'DepMap provides disease-lineage cell-line mutation groups and PRISM AUC aggregates.','cross_site_rule':'Evidence is considered concordant only when the patient biomarker signal and lower-AUC cell-line drug signal point in the same direction; this is hypothesis generation, not patient-level linkage.','privacy_boundary':'Only aggregate statistics are exchanged; raw patient and cell-line rows remain local.'}}
    (ROOT/'DepMap/processed'/f'{disease.lower()}_{gene.lower()}_federated_evidence.json').write_text(json.dumps(out,indent=2))
    print('wrote',out['site'])

