# Federated Evidence Lab

Privacy-preserving biomarker–drug evidence integration across independent hospital and laboratory sites.

## Demonstration cases
- LUAD–EGFR
- SKCM–BRAF
- BRCA–ERBB2
- COAD–BRAF
- GBM–EGFR

## Federated design
Hospital A/B analyze local patient data. Lab A/B analyze local cell-line and PRISM drug-response data. Only aggregate counts, moments, effect sizes, p-values, survival summaries, and quality flags are transferred to the central agent. Raw identifiers and individual measurements remain local.

## Architecture (what each component does)

| Component | Role |
|---|---|
| Xena/raw/ | Local TCGA/Xena clinical, mutation, expression, and survival inputs |
| Xena/scripts/01–09 | LUAD–EGFR patient audit, table construction, validation, and survival analyses |
| Xena/scripts/10_build_skcm_braf_patient_table.py | Build real SKCM–BRAF patient table |
| Xena/scripts/11_analyze_skcm_braf.py | SKCM–BRAF expression, clinical, and survival summaries |
| Xena/scripts/12_build_additional_case_summaries.py | BRCA–ERBB2, COAD–BRAF, and GBM–EGFR patient summaries |
| DepMap/Model.csv | Cell-line metadata and OncoTree disease lineage |
| DepMap/OmicsSomaticMutationsMAF.maf | Local somatic mutation calls used to assign biomarker status |
| DepMap/scripts/07_build_multisite_cell_line_tables.py | Intersect disease lineage, mutation status, metadata, and PRISM availability |
| DepMap/scripts/08_analyze_multisite_drug_response.py | AUC comparisons, controls, permutation tests, and exploratory rankings |
| DepMap/scripts/09_combine_patient_cell_drug_evidence.py | Combine patient and cell–drug aggregate contracts |
| central/01_build_central_agent_results.py | Original LUAD–EGFR central validation workflow |
| central/02_build_multisite_federated_summary.py | Build multisite aggregate federated summary |
| federated_streamlit_app.py | Interactive local dashboard for site selection and benchmark viewing |
| demo-app/ | Federated Evidence Lab website source |
| Xena/processed/, DepMap/processed/, central/ JSON | Aggregate-only site contracts and evidence outputs |

## Outputs
The repository contains reproducible Xena and DepMap analysis scripts, aggregate-only site summaries, federated evidence packages, and presentation materials. Results are research prioritization evidence, not clinical treatment guidance.
