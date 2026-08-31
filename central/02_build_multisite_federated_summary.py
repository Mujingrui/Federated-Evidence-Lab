from pathlib import Path
import json

root = Path('/Users/jingruimu/Documents/ChatGPT/federated AI')
items = []
for p in sorted((root / 'DepMap/processed').glob('*_federated_evidence.json')):
    d = json.loads(p.read_text())
    items.append({'site': d['site'], 'patient_site': d['patient_site'], 'cell_drug_site': d['cell_drug_site'], 'privacy_boundary': d['federated_interpretation']['privacy_boundary']})
out = {'analysis': 'multi-disease drug-cell-patient federated evidence', 'sites': items, 'centralized_benchmark': {'status': 'simulated benchmark', 'definition': 'Centralized benchmark combines raw rows in one environment; federated results use aggregate site contracts.', 'comparison_metrics': ['mutation prevalence', 'expression effect', 'survival hazard ratio', 'mutated-versus-wild-type AUC difference', 'permutation p-value'], 'warning': 'Patient and cell-line records are not directly linked; conclusions are hypothesis-generating.'}}
(root / 'central/multisite_federated_summary.json').write_text(json.dumps(out, indent=2))
print('sites:', [x['site'] for x in items])

