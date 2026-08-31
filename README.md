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

## Outputs
The repository contains reproducible Xena and DepMap analysis scripts, aggregate-only site summaries, federated evidence packages, and a presentation deck. Results are research prioritization evidence, not clinical treatment guidance.
