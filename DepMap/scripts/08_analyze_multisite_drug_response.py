from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path('/Users/jingruimu/Documents/ChatGPT/federated AI/DepMap')
prism=pd.read_csv(ROOT/'Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv',low_memory=False)
idcol=prism.columns[0]
drug_cols=list(prism.columns[7:])
rng=np.random.default_rng(42)
pairs=[('SKCM','BRAF'),('BRCA','ERBB2'),('COAD','BRAF'),('GB','EGFR')]

def perm_p(a,b,n=500):
    if len(a)<2 or len(b)<2:return None
    obs=float(np.mean(a)-np.mean(b)); pooled=np.r_[a,b]; k=len(a); hits=0
    for _ in range(n):
        z=rng.permutation(pooled); hits += abs(z[:k].mean()-z[k:].mean()) >= abs(obs)
    return (hits+1)/(n+1)

for disease,gene in pairs:
    tab=pd.read_csv(ROOT/'processed'/f'{disease.lower()}_{gene.lower()}_cell_line_table.tsv',sep='\t')
    x=prism[prism[idcol].isin(set(tab.DepMap_ID))].copy().set_index(idcol)
    status=tab.set_index('DepMap_ID').mutated
    rows=[]
    for drug in drug_cols:
        vals=pd.to_numeric(x[drug],errors='coerce').dropna()
        vals=vals[vals.index.isin(status.index)]
        a=vals[status.loc[vals.index].eq(1)].to_numpy(); b=vals[status.loc[vals.index].eq(0)].to_numpy()
        if len(a)>=2 and len(b)>=2:
            rows.append({'drug':drug,'mutated_n':len(a),'wild_type_n':len(b),'mutated_mean_auc':float(a.mean()),'wild_type_mean_auc':float(b.mean()),'auc_difference_mut_minus_wt':float(a.mean()-b.mean()),'permutation_p':perm_p(a,b)})
    out=pd.DataFrame(rows, columns=['drug','mutated_n','wild_type_n','mutated_mean_auc','wild_type_mean_auc','auc_difference_mut_minus_wt','permutation_p'])
    if len(out): out=out.sort_values('auc_difference_mut_minus_wt')
    out.to_csv(ROOT/'processed'/f'{disease.lower()}_{gene.lower()}_drug_screen.tsv',sep='\t',index=False)
    controls=['OSIMERTINIB','GEFITINIB','ERLOTINIB','AFATINIB','VEMURAFENIB','DABRAFENIB','TRAMETINIB','TUCATINIB','LAPATINIB']
    control=out[out.drug.str.upper().apply(lambda s:any(c in s for c in controls))].head(30).to_dict('records')
    summary={'site':f'{disease}-{gene}','cell_lines_total':int(len(tab)),'mutated_lines':int(tab.mutated.sum()),'wild_type_lines':int((tab.mutated==0).sum()),'drugs_tested':int(len(out)),'known_control_results':control,'top_exploratory_hits':out.head(20).to_dict('records'),'warning':'Exploratory: small mutated group; AUC associations are not causal.'}
    (ROOT/'processed'/f'{disease.lower()}_{gene.lower()}_drug_summary.json').write_text(json.dumps(summary,indent=2))
    print(f'{disease}-{gene}: drugs={len(out)}, controls={len(control)}, top={out.iloc[0].drug if len(out) else "none"}')

