# ST001639 renal-related signal sensitivity — fixed analysis plan

Recorded 2 October 2026 before fitting the new models. This is a reviewer-motivated, post hoc sensitivity analysis, not a prespecified confirmatory study.

The public metadata contain no age, calibrated serum creatinine, cystatin C, eGFR, renal diagnosis, diet or medication fields usable in the existing model. The metabolomic matrix contains an annotation named `creatinine`, positive for all 648 matrix samples. This is an uncalibrated relative assay signal. It must not be used to calculate eGFR or called clinical renal adjustment.

1. Reconstruct the previously reported five ST001639 models using the original parser, positive log2 transformation, within-study standardisation, covariates (COPD status, sex, race, BMI and smoking category), complete-case handling and OLS/HC3 intervals.
2. Require all five reconstructed coefficients, confidence bounds, P values and complete-case counts to match the saved original results before any new result is accepted.
3. For exactly the same five candidates (indolepropionate, indoleacetate, indolelactate, p-cresol sulfate and hippurate), add the positive log2-standardised `creatinine` relative signal to the same covariate set. Keep the outcome standardisation fixed and verify the model samples remain unchanged and the designs are full rank.
4. Report all five new models, with BH correction across the five additional tests. Do not select a different renal proxy or specification in response to the results. Retain the original 985-feature multiplicity result unchanged. No claim of overall project-wide error control is made.
5. Save aggregate estimates, source-file hashes, package versions and a five-candidate forest figure. Do not export participant records, residuals, source measurements, sample identifiers or the clinical table.
6. Interpretation is restricted to the observed association's sensitivity to this measured signal. Persistence cannot exclude renal confounding; attenuation cannot prove renal causation. Muscle mass, diet, disease-related metabolism and assay effects may affect creatinine. Age and clinical renal function remain unmeasured. Hippurate remains unsuitable as a validated COPD biomarker on these data.
