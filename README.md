# Reproducibility package — tele-oncology conformance and performance

This package accompanies the manuscript *"Evaluating conformance and
performance in oncology pathways via telemedicine using process mining: a study
in a Chilean hospital network."* It provides an **executable, publicly
shareable** implementation of the full analysis pipeline.

The real clinical event log and the institutional source code cannot be
released, because the data are protected patient records and the code is
institutional property. To make the methodology auditable and reproducible
nonetheless, this package ships:

1. a **synthetic, de-identified event log** of fifty illustrative patients that
   reproduces the *structure* of the study data (activity vocabulary,
   telemedicine delivery modality, out-of-model activities, trace repetition and
   right-censoring), and
2. a **self-contained pipeline** that runs the complete conformance and
   performance analysis on that log and writes every result type reported in the
   paper and in the response to reviewers.

The synthetic figures are **illustrative** and are not the study's numbers; the
purpose of the package is to expose the exact procedure, not to restate the
clinical results.

## Contents

| File | Purpose |
|------|---------|
| `make_synthetic_log.py` | Generates `synthetic_event_log.csv` (50 patients, deterministic seed). |
| `synthetic_event_log.csv` | The synthetic event log (schema-identical to the processed study log). |
| `pathway_model.py` | Builds the guideline reference model (acyclic Petri net) in code and classifies alignment moves. |
| `analyses.py` | Modality, exposure, patient-level, transition, regression/survival and sensitivity analyses. |
| `run_pipeline.py` | Orchestrator: runs the whole pipeline end-to-end. |
| `outputs/` | Result tables (CSV) and the sensitivity figure produced by a run. |

## Requirements

```bash
pip install pm4py pandas numpy scipy statsmodels lifelines matplotlib
```

Only the conformance stage requires PM4Py (the same library used in the study);
the remaining analyses use pandas, scipy, statsmodels and lifelines.

## How to run

```bash
python make_synthetic_log.py     # (re)generate the synthetic event log
python run_pipeline.py           # run the analysis; results land in ./outputs/
```

## What the pipeline reproduces

| Stage | Output | Corresponds to |
|-------|--------|----------------|
| 1–2 | `conformance_by_case.csv`, `decomposition_summary.csv`, `decomposition_by_cohort.csv` | Censoring-aware conformance (prefix-alignments) and the decomposition of non-conformance into right-censoring vs genuine deviation, overall and per cohort (Figure 9b). |
| 3 | `modality_association.csv`, `modality_by_activity.csv` | Whether structural deviations concentrate in the remote modality (Cramer's V). |
| 4 | `exposure_adjusted.csv` | Out-of-model activity rate per 100 patient-days and the Poisson incidence-rate ratio for observation time. |
| 5 | `patient_level_tests.csv` | Cohort-vs-rest activity-presence tests at the patient level, with FDR correction. |
| 6 | `transition_times.csv` | Robust transition times (median and IQR alongside the mean). |
| 7 | `regression_conformance.csv`, `regression_duration.csv`, `survival_completion.csv` | Regression of conformance and duration, and Kaplan–Meier time-to-completion (right-censored). |
| 8 | `sensitivity_sweep.csv`, `sensitivity_bootstrap.csv`, `fig_sensitivity.png` | Sensitivity of the decomposition to the definition of the censored tail. |

## Method notes

* **Reference model.** The normative oncology pathway is encoded as an acyclic
  Petri net (`pathway_model.py`): mandatory consultation → oncology consultation
  → oncology committee → (exclusive) treatment → nursing follow-up → discharge,
  with psychology and non-medical care as optional supportive steps. Because the
  model is built in code, no proprietary BPMN file is required.
* **Censoring vs deviation.** For each trace, an optimal alignment is computed
  against the reference model. Synchronous moves are matches; log-moves are
  observed activities not permitted by the model (insertions, i.e. genuine
  deviation); model-moves are required steps not observed — those at the **end**
  of the alignment are attributed to right-censoring (the pathway has not yet
  been completed), while **interior** ones are genuine skipped steps.
* **Abstraction level.** The pipeline runs at the category level, the primary
  analytical level of the study. The identical procedure applies at the finer
  activity level by pointing `run_conformance` at the `ACTIVITY_MODEL` column.
* **Determinism.** The generator and the bootstrap use fixed seeds, so a run is
  fully reproducible.
