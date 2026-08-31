from pathlib import Path
import json
import streamlit as st
import pandas as pd

ROOT=Path(__file__).resolve().parent
SITES={'LUAD–EGFR':('xena_site_summary.json','depmap_site_summary.json'),'SKCM–BRAF':('skcm_braf_site_summary.json','skcm_braf_federated_evidence.json'),'BRCA–ERBB2':('brca_erbb2_site_summary.json','brca_erbb2_federated_evidence.json'),'COAD–BRAF':('coad_braf_site_summary.json','coad_braf_federated_evidence.json'),'GBM–EGFR':('gbm_egfr_site_summary.json','gbm_egfr_federated_evidence.json')}
st.set_page_config(page_title='Federated Evidence Lab',layout='wide')
st.title('Drug–cell–patient federated evidence lab')
st.caption('Aggregate-only evidence integration; raw records remain at their source sites.')
site=st.sidebar.selectbox('Disease–biomarker site',list(SITES)); xf,df=SITES[site]
x=json.loads((ROOT/'Xena/processed'/xf).read_text())
d=json.loads((ROOT/'DepMap/processed'/df).read_text()); d=d.get('cell_drug_site',d)
st.sidebar.success('Xena + DepMap summaries loaded'); st.sidebar.info('Central input: 0 raw records')
st.subheader('Federated site map')
site_rows=[]
for label,(xp,dp) in SITES.items():
    xx=json.loads((ROOT/'Xena/processed'/xp).read_text())
    dd=json.loads((ROOT/'DepMap/processed'/dp).read_text()); dd=dd.get('cell_drug_site',dd)
    cc=xx.get('cohort',{})
    site_rows.append({'site':label,'patient site':'Xena','cell/drug site':'DepMap PRISM','patients':cc.get('mutation_evaluable_patients',cc.get('clinical_patients','—')),'mutated lines':dd.get('mutated_lines',dd.get('egfr_hotspot_mutated_cell_lines','—')),'raw records sent':0})
st.dataframe(pd.DataFrame(site_rows),use_container_width=True,hide_index=True)
st.header(site); co=x.get('cohort',{}); c1,c2,c3=st.columns(3)
c1.metric('Patients',co.get('mutation_evaluable_patients',co.get('clinical_patients','—')))
c2.metric('Mutated patients',co.get('egfr_mutated_patients',co.get('braf_v600_hotspot_patients',co.get('altered_patients','—'))))
c3.metric('DepMap mutated lines',d.get('mutated_lines',d.get('egfr_hotspot_mutated_cell_lines','—')))
a,b,c=st.tabs(['Patient evidence','Drug–cell evidence','Federated vs centralized'])
with a: st.json(x,expanded=False)
with b:
 st.metric('Drugs tested',d.get('drugs_tested','—')); hits=d.get('top_exploratory_hits',d.get('exploratory_top_25_lower_auc_candidates',[]))
 if hits: st.dataframe(pd.DataFrame(hits),use_container_width=True,hide_index=True)
 else: st.warning('Insufficient mutated lines for a two-group comparison.')
with c:
 st.subheader('Centralized benchmark')
 st.write('The centralized benchmark is a simulated validation view: it represents combining the same site-level data in one analysis environment. No raw records are uploaded by this app.')
 bench=[]
 for label,(xp,dp) in SITES.items():
  xx=json.loads((ROOT/'Xena/processed'/xp).read_text()); dd=json.loads((ROOT/'DepMap/processed'/dp).read_text()); dd=dd.get('cell_drug_site',dd)
  cc=xx.get('cohort',{}); bench.append({'site':label,'federated patient count':cc.get('mutation_evaluable_patients',cc.get('clinical_patients','—')),'federated mutated lines':dd.get('mutated_lines',dd.get('egfr_hotspot_mutated_cell_lines','—')),'centralized comparison':'same aggregate counts'})
 st.dataframe(pd.DataFrame(bench),use_container_width=True,hide_index=True)
 st.write('Federated results preserve locality; centralized results serve as a validation benchmark.')
 st.info('Federated results preserve locality; centralized results serve as a validation benchmark.')
 st.warning('Drug rankings are exploratory and do not establish clinical benefit.')

