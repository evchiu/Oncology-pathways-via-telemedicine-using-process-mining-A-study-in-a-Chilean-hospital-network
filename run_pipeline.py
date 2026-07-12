#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Conformance-and-performance pipeline for tele-oncology care pathways.

Self-contained, publicly shareable implementation of the analysis reported in
the manuscript. It runs end-to-end on the synthetic event log produced by
``make_synthetic_log.py`` and reproduces every result type used in the paper and
the response to reviewers:

  1-2. Censoring-aware conformance via optimal alignments against a
       guideline-derived reference model, and decomposition of non-conformance
       into right-censoring versus genuine structural deviation.
  3.   Delivery-modality analysis (telemedicine vs in-person), Cramer's V.
  4.   Exposure-adjusted performance (out-of-model rate per 100 patient-days).
  5.   Patient-level independence tests with Cramer's V and FDR correction.
  6.   Robust transition times (median and IQR alongside the mean).
  7.   Regression and survival models (conformance, duration, completion).
  8.   Sensitivity analysis of the censoring/deviation decomposition.

Alignments are computed with PM4Py, the same library used in the study. The
reference model is built in code (see ``pathway_model.py``), so no proprietary
BPMN file is needed.

Usage:
    python make_synthetic_log.py        # writes synthetic_event_log.csv
    python run_pipeline.py               # writes outputs/*.csv and a figure
"""
from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from pm4py.objects.conversion.log import converter as log_converter
from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments

from pathway_model import build_reference_net, classify
import analyses

LOG_FILE = "synthetic_event_log.csv"
OUT = "outputs"
os.makedirs(OUT, exist_ok=True)


def run_conformance(df):
    """Align every trace against the reference model and classify the moves."""
    net, im, fm = build_reference_net()
    work = df.rename(columns={"ANON.ID": "case:concept:name",
                              "CATEGORY_MODEL": "concept:name",
                              "ts": "time:timestamp"})
    work = work.dropna(subset=["concept:name"]).sort_values(
        ["case:concept:name", "time:timestamp"])
    cohort = df.groupby("ANON.ID")["cohorte_C"].first().to_dict()
    event_log = log_converter.apply(
        work, variant=log_converter.Variants.TO_EVENT_LOG)

    records = []
    for trace in event_log:
        cid = trace.attributes["concept:name"]
        result = alignments.apply(trace, net, im, fm)
        rec = classify(result["alignment"])
        rec["case"] = cid
        rec["cohorte"] = cohort.get(cid, "?")
        rec["std_fitness"] = round(result["fitness"], 4)
        rec["trace_len"] = len(trace)
        records.append(rec)
    out = pd.DataFrame(records)
    out.to_csv(f"{OUT}/conformance_by_case.csv", index=False)
    return out


def decomposition_summary(R):
    total = R["genuine_dev"].sum() + R["censoring"].sum()
    summary = pd.DataFrame([dict(
        n=len(R),
        prefix_conform_pct=round(100 * R["prefix_conform"].mean(), 2),
        mean_prefix_fitness=round(R["prefix_fitness"].mean(), 4),
        mean_strict_fitness=round(R["std_fitness"].mean(), 4),
        genuine_units=int(R["genuine_dev"].sum()),
        censoring_units=int(R["censoring"].sum()),
        pct_genuine=round(100 * R["genuine_dev"].sum() / total, 1) if total else np.nan,
        pct_censoring=round(100 * R["censoring"].sum() / total, 1) if total else np.nan)])
    summary.to_csv(f"{OUT}/decomposition_summary.csv", index=False)
    return summary


def decomposition_by_cohort(R):
    """Per-cohort split of non-conformance into genuine deviation vs censoring
    (the breakdown shown in Figure 9b of the manuscript)."""
    rows = []
    for coh, g in R.groupby("cohorte"):
        total = g["genuine_dev"].sum() + g["censoring"].sum()
        rows.append(dict(
            cohort=coh, n=len(g),
            genuine_units=int(g["genuine_dev"].sum()),
            censoring_units=int(g["censoring"].sum()),
            pct_genuine=round(100 * g["genuine_dev"].sum() / total, 2) if total else np.nan,
            pct_censoring=round(100 * g["censoring"].sum() / total, 2) if total else np.nan))
    out = pd.DataFrame(rows).sort_values("n", ascending=False)
    out.to_csv(f"{OUT}/decomposition_by_cohort.csv", index=False)
    return out


def main():
    df = pd.read_csv(LOG_FILE, sep=";", encoding="utf-8-sig")
    df["ts"] = pd.to_datetime(df["FEC_CREACION"], utc=True, errors="coerce")
    df = df.dropna(subset=["ts"]).sort_values(["ANON.ID", "ts"])

    print("=" * 68)
    print(" Tele-oncology conformance-and-performance pipeline (synthetic run)")
    print("=" * 68)
    print(f" patients={df['ANON.ID'].nunique()}  events={len(df)}")

    print("\n[1-2] Censoring-aware conformance and decomposition")
    conformance = run_conformance(df)
    print(decomposition_summary(conformance).to_string(index=False))
    decomposition_by_cohort(conformance)

    print("\n[3] Delivery-modality analysis")
    analyses.modality_analysis(df)

    print("[4] Exposure-adjusted performance")
    analyses.exposure_analysis(df)

    print("[5] Patient-level tests (Cramer's V, FDR)")
    analyses.patient_level_tests(df)

    print("[6] Robust transition times (median / IQR)")
    analyses.transition_times(df)

    print("[7] Regression and survival")
    analyses.regression_survival(df, conformance)

    print("[8] Sensitivity of the decomposition")
    analyses.sensitivity_analysis(conformance)

    print("\nDone. All results written to ./outputs/")


if __name__ == "__main__":
    main()
