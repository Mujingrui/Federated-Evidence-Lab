 # Project progress

## Completed
- Downloaded and audited TCGA/Xena clinical, mutation, expression, and survival data.
- Built LUAD–EGFR patient and DepMap workflows.
- Added real Xena SKCM–BRAF analysis.
- Added real Xena summaries for BRCA–ERBB2, COAD–BRAF, and GBM–EGFR.
- Downloaded full DepMap somatic mutation MAF and matched mutations to Model.csv and PRISM.
- Built disease-specific cell-line tables for SKCM–BRAF, BRCA–ERBB2, COAD–BRAF, and GBM–EGFR.
- Screened PRISM compounds and generated aggregate drug summaries.
- Combined patient and cell–drug summaries into federated evidence packages.

## Validation
- Positive and negative drug controls.
- Permutation tests for AUC differences.
- Leave-one-cell-line-out robustness for supported cases.
- Federated-versus-centralized parity benchmark.
- Explicit inconclusive status for small mutated groups.

## Current limitations
- Patient and cell-line records are not directly linked.
- BRCA–ERBB2 and GBM–EGFR have very small mutated cell-line groups.
- Results are exploratory research prioritization, not clinical treatment guidance.

## Reproducibility
The local project contains Python scripts under Xena/scripts, DepMap/scripts, and central, plus aggregate JSON contracts. Raw data should remain local in a production federated deployment.
