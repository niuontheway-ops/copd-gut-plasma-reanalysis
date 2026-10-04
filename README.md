# Exploratory COPD gut–plasma reanalysis: code and derived-output snapshot

This archive contains analysis/audit/plotting Python scripts, pinned package lists, aggregate result CSVs, and audit JSONs. It deliberately excludes all `sources/` files and participant-level data; obtain original inputs directly from their custodians using the accessions and DOI links in the manuscript's data-availability statement. A GitHub URL is not a repository DOI or evidence of journal acceptance.

Run on Python 3.12 after obtaining the public source files under the relative paths expected by each script. Install `requirements-pilot.txt` and `requirements-hummanet.txt` into the project-local `.deps` directory, then run in this order: `analyze_pilot.py`, `analyze_copd_stool.py`, `analyze_spiromics_fev1.py`, `analyze_aecopd_paired.py`, `analyze_xuq_gut_serum.py`, `audit_xuq_original_source.py`, `audit_published_ltrc.py`, `analyze_bowerman_stool.py`, `plot_cohort_flow.py`. The separate `analyze_mtbls9119.py` output is an eligibility audit, not a manuscript result. `fetch_public.py` and `audit_hummanet.py` help inspect/download public deposits but do not override their reuse terms.

The scripts expect a project layout with `scripts/`, `sources/`, and `results/`. Some source files have reuse or attribution conditions. The authors have confirmed source-data reuse for this analysis; that confirmation does not grant permission to redistribute third-party participant-level files. The Arivale repository has no visible redistribution licence; the LTRC publisher supplement was used only for aggregate literature checking. Do not add any third-party source file to a public deposit solely because it was downloaded locally. A successful rerun is not an independent scientific validation.

This is an exploratory analysis snapshot, not a validated biomarker or a same-person COPD gut–plasma–lung analysis. The authors have approved the manuscript and confirmed data reuse. These author confirmations do not establish that every analysis has been independently rerun in a clean environment, and this repository does not claim such validation. A GitHub commit is a versioned code record, not a repository DOI or evidence of journal acceptance.

## Reviewer addendum published 4 October 2026

`reviewer_addendum_2026-10-02/` contains the unchanged analysis script and fixed plan recorded on 2 October, all ten aggregate model estimates in S21, the S5 figure in PNG, TIFF and vector PDF formats, the run audit, and a file manifest. The earlier 26 September ZIP remains unchanged. No source matrices, participant identifiers, residuals, authors' private review materials, ethics documents or credentials are included in this addendum.

The five original ST001639 candidate models were reproduced before fitting the five fixed models with an additional uncalibrated creatinine relative signal. The hippurate coefficient changes from 0.219943 to 0.220825. This is a post hoc sensitivity analysis, not adjustment for calibrated clinical creatinine or eGFR; renal-function and age confounding remain unresolved. BH correction is across five candidates separately within each specification, not project-wide error control. The earlier 985-feature correction is not replaced or recalculated here.

### Reproduce the addendum

First extract the original ZIP into a separate project directory, preserving its `scripts/` and `results/` layout. Obtain the four ST001639 assay matrices and its sample metadata directly from Metabolomics Workbench under the filenames recorded in the audit, subject to the applicable terms. Do not upload those source files to this repository. The new audit records the tested Python and package versions, source hashes and the exact original helper/results hashes.

The old `requirements-pilot.txt` is not, by itself, a complete clean-environment guarantee: `pandas.read_html()` also requires a compatible HTML parser, such as `lxml`. Install the pinned original analysis dependencies and an appropriate parser in the environment you actually use. A clean-environment installation was not newly tested for this release.

```text
python reviewer_addendum_2026-10-02/analyze_renal_signal_sensitivity.py --source-project PATH_TO_EXTRACTED_ORIGINAL_PROJECT --output-dir PATH_TO_NEW_AGGREGATE_OUTPUTS
```

The script verifies original-model agreement and exports only aggregate estimates, a forest figure and run metadata. The two original helper/result files needed by this addendum are already in the old ZIP; their hashes are recorded in the new audit and were checked before release. No Arivale or additional discovery-cohort analysis is required for this addendum. Plotting uses Codex-assisted Python code applied to observed aggregate estimates, not generated participant measurements. File checksums are listed in `SHA256SUMS.txt`.
