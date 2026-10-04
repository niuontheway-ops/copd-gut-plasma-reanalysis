"""Fixed, post hoc renal-related metabolomics-signal sensitivity; aggregate output only."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-project", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.source_project.resolve(strict=True)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(root / ".deps"))
    spec = importlib.util.spec_from_file_location("original_pilot", root / "scripts/analyze_pilot.py")
    helper = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(helper)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import statsmodels
    from statsmodels.stats.multitest import multipletests

    files = [root / "sources" / f"ST001639-{name}.txt"
             for name in ("AN002681", "AN002682", "AN002683", "AN002684")]
    meta_path = root / "sources/ST001639-samples.html"
    original_path = root / "results/copd_metabolites_adjusted.csv"
    matrices = [helper.read_workbench_matrix(path) for path in files]
    keys = matrices[0][0].index
    assert keys.is_unique and len(keys) == 648
    for matrix, factors in matrices[1:]:
        assert set(matrix.index) == set(keys)
        assert factors.sort_index().equals(matrices[0][1].sort_index())
    features = pd.concat([matrix.reindex(keys) for matrix, _ in matrices], axis=1)
    assert features.columns.is_unique
    meta = next(table for table in pd.read_html(meta_path) if "Sample name" in table.columns)
    metadata_fields = list(meta.columns)
    meta = meta.set_index("Sample name")
    assert meta.index.is_unique and set(keys).issubset(meta.index)
    factors = matrices[0][1].str.split(":", n=1).str[-1]
    assert factors.loc[keys].fillna("NA").eq(meta.loc[keys, "COPD"].fillna("NA")).all()
    clinical = meta.loc[keys, ["COPD", "Sample_Data:Gender", "Race", "BMI", "Smoking Status"]].rename(
        columns={"Sample_Data:Gender": "sex", "Race": "race", "BMI": "bmi", "Smoking Status": "smoking"})
    clinical["bmi"] = pd.to_numeric(clinical["bmi"], errors="coerce")
    clinical = clinical.loc[clinical["COPD"].isin(["case", "control"])].copy()
    assert len(clinical) == 636
    clinical["copd"] = clinical["COPD"].eq("case").astype(int)
    assert "creatinine" in features.columns
    renal_signal = pd.to_numeric(features["creatinine"], errors="coerce")
    assert renal_signal.gt(0).all() and renal_signal.notna().all()
    clinical["creatinine_signal_zlog2"] = helper.zlog(renal_signal.loc[clinical.index])
    fields = ["copd", "sex", "race", "bmi", "smoking"]
    designs = {"original_covariates": helper.design(clinical, fields),
               "plus_creatinine_relative_signal": helper.design(clinical, fields + ["creatinine_signal_zlog2"])}
    expected = pd.read_csv(original_path).set_index("metabolite")
    rows = []
    matched = []
    for name in helper.TARGETS:
        outcome = helper.zlog(features.loc[clinical.index, name])
        complete_keys = []
        for label, design in designs.items():
            complete = design.assign(outcome=outcome).replace([np.inf, -np.inf], np.nan).dropna()
            assert np.linalg.matrix_rank(complete[design.columns].to_numpy()) == len(design.columns)
            complete_keys.append(set(complete.index))
            result = helper.model(outcome, design, "copd")
            assert result["n_complete"] == len(complete)
            if label == "original_covariates":
                assert result["n_complete"] == int(expected.loc[name, "n_complete"])
                for field in ("beta", "lo", "hi", "p"):
                    assert np.isclose(result[field], expected.loc[name, field], rtol=1e-10, atol=1e-12), (name, field)
                matched.append(name)
            rows.append({"metabolite": name, "specification": label, "n_complete": result["n_complete"],
                         "beta": result["beta"], "lo": result["lo"], "hi": result["hi"], "p": result["p"]})
        assert complete_keys[0] == complete_keys[1]
    results = pd.DataFrame(rows)
    for label in designs:
        use = results["specification"].eq(label)
        results.loc[use, "q_bh_5_within_specification"] = multipletests(results.loc[use, "p"], method="fdr_bh")[1]
    original_q = results.loc[results.specification.eq("original_covariates")].set_index("metabolite")
    assert np.allclose(original_q.loc[list(helper.TARGETS), "q_bh_5_within_specification"],
                       expected.loc[list(helper.TARGETS), "q_bh_5"], rtol=1e-10, atol=1e-12)
    results.to_csv(out / "supplementary_table_s21_renal_signal_sensitivity.csv", index=False)
    hippurate = results.loc[results.metabolite.eq("hippurate")].copy()
    baseline_beta, adjusted_beta = hippurate.beta.to_numpy()
    summary = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "study": "ST001639", "analysis": "post_hoc_creatinine_relative_signal_sensitivity",
        "plan": "RENAL_SIGNAL_SENSITIVITY_PLAN.md", "matrix_samples": len(keys),
        "copd_labelled_samples": len(clinical), "metadata_fields": metadata_fields,
        "creatinine_signal_positive_samples": int(renal_signal.gt(0).sum()),
        "creatinine_units": "uncalibrated relative metabolomics signal, not clinical concentration or eGFR",
        "baseline_match_targets": matched, "baseline_match_tolerances": {"rtol": 1e-10, "atol": 1e-12},
        "same_complete_case_samples_between_specifications": True,
        "all_designs_full_rank": True, "multiplicity": "BH across 5 candidates separately for each specification; not project-wide error control",
        "hippurate": hippurate.to_dict(orient="records"),
        "hippurate_coefficient_relative_attenuation": float(1 - adjusted_beta / baseline_beta),
        "versions": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
                     "statsmodels": statsmodels.__version__, "matplotlib": matplotlib.__version__},
        "sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in [*files, meta_path, original_path, root / "scripts/analyze_pilot.py", Path(__file__)]},
        "limitations": ["No age or clinical renal function", "Cannot estimate eGFR from relative signal",
                        "No proof of renal cause, no exclusion of renal confounding", "Post hoc analysis; no biomarker validation"],
        "participant_level_output": False,
    }
    (out / "renal_signal_sensitivity_audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    labels = ["Indolepropionate", "Indoleacetate", "Indolelactate", "p-Cresol sulfate", "Hippurate"]
    y = np.arange(len(labels))[::-1]
    for (label, delta, colour, display) in [("original_covariates", .12, "#315f86", "Original covariates"),
                                          ("plus_creatinine_relative_signal", -.12, "#b66b33", "+ relative creatinine signal")]:
        part = results.loc[results.specification.eq(label)].set_index("metabolite").loc[list(helper.TARGETS)]
        ax.errorbar(part.beta, y + delta, xerr=np.array([part.beta - part.lo, part.hi - part.beta]),
                    fmt="o", ms=5, capsize=3, color=colour, label=display)
    ax.axvline(0, color="#666666", lw=.8, ls="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel("COPD-minus-control difference (log2 signal SD; HC3 95% CI)")
    ax.set_title("Post hoc sensitivity to an uncalibrated creatinine signal", fontsize=11)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.19), ncol=2, fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "supplementary_figure_s5_renal_signal_sensitivity.png", dpi=300)
    fig.savefig(out / "supplementary_figure_s5_renal_signal_sensitivity.tiff", dpi=300,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(out / "supplementary_figure_s5_renal_signal_sensitivity.pdf")
    plt.close(fig)
    print(json.dumps({"baseline_match": matched, "hippurate": summary["hippurate"],
                      "output": str(out), "exported_participant_records": 0}, indent=2))


if __name__ == "__main__":
    main()
