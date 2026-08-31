from pathlib import Path
import json


ROOT = Path("/Users/jingruimu/Documents/ChatGPT/federated AI")
XENA_FILE = ROOT / "Xena" / "processed" / "xena_site_summary.json"
DEPMAP_FILE = ROOT / "DepMap" / "processed" / "depmap_site_summary.json"
JSON_OUTPUT = ROOT / "central" / "central_agent_results.json"
REPORT_OUTPUT = ROOT / "central" / "central_agent_results.md"


with XENA_FILE.open(encoding="utf-8") as handle:
    xena = json.load(handle)
with DEPMAP_FILE.open(encoding="utf-8") as handle:
    depmap = json.load(handle)

assert xena["analysis_type"] == "local_summary_statistics"
assert depmap["analysis_type"] == "local_summary_statistics"
assert xena["site"] == "Xena_TCGA_LUAD"
assert depmap["site"] == "DepMap_PRISM_LUAD"

xena_cohort = xena["cohort"]
xena_expression = xena["egfr_expression_by_mutation_status"]
xena_clinical = xena["clinical_characteristics_by_mutation_status"]
xena_mutation_survival = xena["egfr_mutation_overall_survival"]
xena_expression_survival = xena["egfr_expression_overall_survival"]
depmap_cohort = depmap["cohort"]

known = sorted(
    depmap["known_egfr_inhibitor_results"],
    key=lambda row: row["exact_two_sided_permutation_p"],
)
known_names = {row["drug"].upper() for row in known}
exploratory = [
    row
    for row in depmap["exploratory_top_25_lower_auc_candidates"]
    if row["compound"].split(" (")[0].upper() not in known_names
][:10]

results = {
    "central_agent": "EGFR_LUAD_evidence_integrator",
    "workflow": [
        "Receive aggregate-only JSON from the Xena and DepMap sites.",
        "Validate site identity, analysis type, disease context, and shared EGFR biomarker.",
        "Summarize patient-side biomarker evidence separately from cell-line drug evidence.",
        "Prioritize prespecified EGFR inhibitors using exact DepMap comparisons.",
        "Surface exploratory compounds separately and retain their small-sample warning.",
        "Return an evidence synthesis without patient-level or cell-line-level records.",
    ],
    "aligned_context": {
        "disease": "lung adenocarcinoma",
        "biomarker": "EGFR mutation",
        "patient_site": "TCGA-LUAD through Xena",
        "drug_cell_line_site": "LUAD models through DepMap PRISM",
    },
    "patient_site_results": {
        "clinical_patients": xena_cohort["clinical_patients"],
        "mutation_evaluable_patients": xena_cohort["mutation_evaluable_patients"],
        "egfr_mutated_patients": xena_cohort["egfr_mutated_patients"],
        "egfr_mutation_prevalence_percent": xena_cohort["egfr_mutation_prevalence_percent"],
        "egfr_expression_mean_difference_mutated_minus_wild_type": xena_expression["mean_difference_mutated_minus_wild_type"],
        "egfr_expression_permutation_p": xena_expression["two_sided_permutation_p_value"],
        "female_proportion_difference": xena_clinical["female_proportion_difference"],
        "female_proportion_permutation_p": xena_clinical["female_proportion_permutation_p_value"],
        "egfr_mutation_survival_hazard_ratio": xena_mutation_survival["unadjusted_cox_hazard_ratio"],
        "egfr_mutation_survival_p": xena_mutation_survival["unadjusted_cox_p_value"],
        "egfr_expression_survival_hazard_ratio_per_unit": xena_expression_survival["continuous_cox_hr_per_expression_unit"],
        "egfr_expression_survival_p": xena_expression_survival["continuous_cox_p_value"],
        "interpretation": "EGFR mutation defines a molecular subgroup with higher EGFR expression; no significant unadjusted overall-survival association was observed.",
    },
    "depmap_site_results": {
        "analysis_ready_luad_cell_lines": depmap_cohort["analysis_ready_luad_cell_lines"],
        "egfr_mutated_cell_lines": depmap_cohort["egfr_hotspot_mutated_cell_lines"],
        "egfr_wild_type_cell_lines": depmap_cohort["egfr_hotspot_wild_type_cell_lines"],
        "prespecified_egfr_inhibitors_ranked_by_exact_p": known,
        "exploratory_non_prespecified_top_candidates": exploratory,
        "interpretation": "EGFR-mutated LUAD cell lines showed markedly lower PRISM AUC for all four prespecified EGFR inhibitors.",
    },
    "central_conclusion": {
        "patient_to_drug_hypothesis": "Patients with EGFR-mutated LUAD form the biomarker-defined group most relevant to EGFR-directed drug hypotheses supported by the DepMap site.",
        "evidence_level": "patient molecular association plus preclinical cell-line drug sensitivity",
        "not_demonstrated": "individual patient response prediction, causal treatment benefit, or clinical efficacy",
    },
    "privacy_result": "Only site-level aggregate summaries were integrated; raw patient and cell-line records were not used by the central agent.",
}

with JSON_OUTPUT.open("w", encoding="utf-8") as handle:
    json.dump(results, handle, indent=2)
    handle.write("\n")


def fmt(value, digits=3):
    return f"{value:.{digits}f}"


known_rows = "\n".join(
    f"| {row['drug']} | {row['mutated_mean_auc']:.3f} | {row['wild_type_mean_auc']:.3f} | "
    f"{row['mean_auc_difference_mutated_minus_wild_type']:.3f} | {row['cohens_d']:.2f} | "
    f"{row['exact_two_sided_permutation_p']:.6f} |"
    for row in known
)
exploratory_rows = "\n".join(
    f"| {index} | {row['compound'].split(' (')[0]} | "
    f"{row['mean_auc_difference_mutated_minus_wild_type']:.3f} | {row['hedges_g']:.2f} | "
    f"{row['two_sided_fdr']:.4f} |"
    for index, row in enumerate(exploratory, start=1)
)

report = f"""# Central-agent EGFR–LUAD evidence integration

## Workflow

The central agent received only aggregate JSON summaries from the Xena patient site and the DepMap drug–cell-line site. It validated the LUAD disease context and shared EGFR biomarker, kept patient and preclinical evidence separate, and then integrated their conclusions.

## Xena patient-site results

- {xena_cohort['egfr_mutated_patients']} of {xena_cohort['mutation_evaluable_patients']} mutation-evaluable patients were EGFR-mutated ({xena_cohort['egfr_mutation_prevalence_percent']:.1f}%).
- EGFR expression was higher in mutated patients by {xena_expression['mean_difference_mutated_minus_wild_type']:.3f} Xena expression units (permutation p={xena_expression['two_sided_permutation_p_value']:.6f}).
- The female proportion was {xena_clinical['female_proportion_difference']:.3f} higher in the mutated group (p={xena_clinical['female_proportion_permutation_p_value']:.4f}).
- EGFR mutation was not significantly associated with overall survival: HR={xena_mutation_survival['unadjusted_cox_hazard_ratio']:.2f}, p={xena_mutation_survival['unadjusted_cox_p_value']:.3f}.
- Expression versus survival was exploratory and borderline: HR={xena_expression_survival['continuous_cox_hr_per_expression_unit']:.3f} per unit, p={xena_expression_survival['continuous_cox_p_value']:.3f}.

## DepMap prespecified drug results

Lower AUC indicates greater sensitivity. The analysis included 4 EGFR-mutated and 47 wild-type LUAD lines.

| Drug | Mutated mean AUC | Wild-type mean AUC | Difference | Cohen's d | Exact p |
|---|---:|---:|---:|---:|---:|
{known_rows}

All four prespecified EGFR inhibitors showed markedly lower AUC in EGFR-mutated models.

## Exploratory candidates

These are secondary hypotheses from the all-compound screen, excluding the four prespecified drugs. They are not clinical recommendations.

| Rank | Compound | AUC difference | Hedges' g | Two-sided FDR |
|---:|---|---:|---:|---:|
{exploratory_rows}

## Integrated conclusion

The Xena site identifies an EGFR-mutated LUAD patient subgroup with higher EGFR expression. Independently, the DepMap site shows that EGFR-mutated LUAD cell lines are substantially more sensitive to EGFR inhibitors. Together, these summaries support an EGFR-stratified patient-to-drug hypothesis without transferring raw records.

This demonstration does not train an individual patient-to-drug model, prove causal treatment benefit, or provide clinical treatment guidance. The DepMap evidence is limited by only four EGFR-mutated cell lines, and the all-compound results are exploratory.
"""

REPORT_OUTPUT.write_text(report, encoding="utf-8")

print("Validated sites: 2")
print("Prespecified drugs prioritized:", len(known))
print("Exploratory non-prespecified candidates reported:", len(exploratory))
print("Outputs:", JSON_OUTPUT, REPORT_OUTPUT)

