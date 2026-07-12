#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Synthetic event-log generator for the tele-oncology conformance study.

This script produces a fully synthetic, de-identified event log of fifty
illustrative patients that reproduces the *structure* of the clinical log
analysed in the manuscript -- activity vocabulary, telemedicine delivery
modality, out-of-model activities, trace repetition and right-censoring --
without containing any real patient information.

The resulting file (``synthetic_event_log.csv``) is schema-compatible with
``run_pipeline.py``, so the complete conformance-and-performance analysis can be
run end-to-end on data that may be shared publicly.

Output columns (identical to the processed study log):
    ANON.ID, ACTIVITY, ACTIVITY_MODEL, CATEGORY, CATEGORY_MODEL,
    FEC_CREACION, SUBSISTEMA, MEDICO, ESPECIALIDAD, COD_CIE_10,
    COD_CIE_10_clean, DES_CIE_10, CENTRO, cohorte_A, cohorte_C, quarter_caso
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SEED = 20260710
OUTFILE = "synthetic_event_log.csv"
N_PATIENTS = 50

# --------------------------------------------------------------------------- #
# Reference (guideline) vocabulary                                            #
# --------------------------------------------------------------------------- #
# Categories permitted by the normative BPMN model, in guideline order. The
# first three steps and the discharge are mandatory; the treatment step is an
# exclusive choice; psychology and non-medical care are optional supportive
# activities.
MANDATORY_HEAD = ["SPECIALIST MEDICAL CONSULTATION",
                  "ONCOLOGY CONSULTATION",
                  "ONCOLOGY COMMITTEE"]
TREATMENT_CHOICE = ["SYSTEMIC TREATMENT", "SURGERY", "RADIOTHERAPY COMMITTEE"]
OPTIONAL_SUPPORT = ["PSYCHOLOGY", "NON-MEDICAL CARE"]
FOLLOW_UP = ["NURSING CARE"]
TERMINAL = "HOSPITAL DISCHARGE"

# Activities present in real logs but absent from the reference model.
OUT_OF_MODEL = ["EMERGENCY CARE", "GENERAL MEDICINE CONSULTATION", "PROCEDURES"]

# Categories with a remote (telemedicine) delivery variant. The variant label
# is what appears in the raw CATEGORY column; the base label is what the
# preprocessing step maps it to (CATEGORY_MODEL).
TELE_VARIANT = {
    "SPECIALIST MEDICAL CONSULTATION": "SPECIALIST TELEMEDICAL CONSULTATION",
    "ONCOLOGY CONSULTATION":           "TELEMEDICINE ONCOLOGY CONSULTATION",
    "NURSING CARE":                    "TELEMEDICAL NURSING CARE",
    "PSYCHOLOGY":                      "TELEMEDICINE PSYCHOLOGY",
    "NON-MEDICAL CARE":                "NON-MEDICAL TELEMEDICINE CARE",
    "GENERAL MEDICINE CONSULTATION":   "GENERAL TELEMEDICINE CONSULTATION",
}

# Fine-grained activity labels (ACTIVITY is more specific than CATEGORY).
ACTIVITY_LABEL = {
    "SPECIALIST MEDICAL CONSULTATION": ["Digestive Surgery Consultation",
                                        "Endocrinology Consultation",
                                        "Bronchopulmonary Consultation"],
    "ONCOLOGY CONSULTATION":           ["Oncology Consultation"],
    "ONCOLOGY COMMITTEE":              ["Oncology Committee"],
    "SYSTEMIC TREATMENT":              ["Chemotherapy Administration"],
    "SURGERY":                         ["Oncological Surgery"],
    "RADIOTHERAPY COMMITTEE":          ["Radiotherapy Committee"],
    "NURSING CARE":                    ["Nursing Care"],
    "PSYCHOLOGY":                      ["Psychology Session"],
    "NON-MEDICAL CARE":                ["Nutritionist Consultation"],
    "HOSPITAL DISCHARGE":              ["Hospital Discharge"],
    "EMERGENCY CARE":                  ["Emergency Care"],
    "GENERAL MEDICINE CONSULTATION":   ["General Medicine Consultation"],
    "PROCEDURES":                      ["Diagnostic Procedure"],
}

SPECIALTY = {
    "SPECIALIST MEDICAL CONSULTATION": "SPECIALIST PHYSICIAN",
    "ONCOLOGY CONSULTATION":           "ADULT ONCOLOGY",
    "ONCOLOGY COMMITTEE":              "ADULT ONCOLOGY",
    "SYSTEMIC TREATMENT":              "ADULT ONCOLOGY",
    "SURGERY":                         "ONCOLOGICAL SURGERY",
    "RADIOTHERAPY COMMITTEE":          "RADIOTHERAPY",
    "NURSING CARE":                    "NURSE",
    "PSYCHOLOGY":                      "PSYCHOLOGY",
    "NON-MEDICAL CARE":                "NUTRITION",
    "HOSPITAL DISCHARGE":              "ADULT ONCOLOGY",
    "EMERGENCY CARE":                  "EMERGENCY",
    "GENERAL MEDICINE CONSULTATION":   "GENERAL MEDICINE",
    "PROCEDURES":                      "DIAGNOSTIC IMAGING",
}

CENTERS = ["CENTER_A", "CENTER_B", "TELEMEDICINE", "CANCER_CENTER", "HOSPITAL"]

# Cohorts: ICD-10 primary code, description and the two cohort groupings used in
# the study (cohorte_C is the six-way split; cohorte_A the three-way one), with
# the relative frequency used to assign patients (large colon/breast/other
# cohorts and small exploratory lung/pancreas/rectum cohorts, as in the study).
COHORTS = {
    "Colon":    ("C18", "MALIGNANT NEOPLASM OF COLON",    "Colon", 0.26),
    "Mama":     ("C50", "MALIGNANT NEOPLASM OF BREAST",   "Mama",  0.20),
    "Otros":    ("C73", "MALIGNANT NEOPLASM OF THYROID",  "Otros", 0.24),
    "Recto":    ("C20", "MALIGNANT NEOPLASM OF RECTUM",   "Otros", 0.10),
    "Pulmon":   ("C34", "MALIGNANT NEOPLASM OF BRONCHUS", "Otros", 0.10),
    "Pancreas": ("C25", "MALIGNANT NEOPLASM OF PANCREAS", "Otros", 0.10),
}

TARGET_TELE_SHARE = 0.185   # overall fraction of remote events in the study


def _quarter(ts: pd.Timestamp) -> str:
    return f"{ts.year}Q{(ts.month - 1) // 3 + 1}"


def build_trace(rng: np.random.Generator, is_completer: bool) -> list[str]:
    """Return the ordered list of CATEGORY_MODEL values for one patient."""
    seq: list[str] = []

    # Entry activity: ~45% of patients enter at the specialist consultation;
    # the remainder enter through another node (structural entry variability).
    if rng.random() < 0.45:
        seq.append("SPECIALIST MEDICAL CONSULTATION")
    else:
        seq.append(str(rng.choice(["ONCOLOGY CONSULTATION",
                                   "GENERAL MEDICINE CONSULTATION"])))

    # Mandatory head, each step occasionally skipped (interior model-move).
    for act in MANDATORY_HEAD:
        if act in seq:
            continue
        if rng.random() < 0.85:
            seq.append(act)
        if act in ("SPECIALIST MEDICAL CONSULTATION", "ONCOLOGY CONSULTATION"):
            for _ in range(int(rng.integers(0, 3))):
                seq.append(act)

    # Optional supportive activities.
    for act in OPTIONAL_SUPPORT:
        if rng.random() < 0.45:
            seq.append(act)

    # Treatment: exclusive choice, with realistic repetition (e.g. chemo
    # cycles) -- the dominant source of trace repetition in oncology logs.
    treatment = str(rng.choice(TREATMENT_CHOICE, p=[0.7, 0.2, 0.1]))
    n_cycles = int(rng.integers(3, 9)) if treatment == "SYSTEMIC TREATMENT" else 1
    for _ in range(n_cycles):
        seq.append(treatment)
        if rng.random() < 0.5:
            seq.append("NURSING CARE")

    # Sparse out-of-model events (emergency visits, procedures).
    for _ in range(int(rng.poisson(0.4))):
        pos = int(rng.integers(1, len(seq) + 1))
        seq.insert(pos, str(rng.choice(OUT_OF_MODEL, p=[0.5, 0.3, 0.2])))

    # Follow-up nursing care.
    if rng.random() < 0.7:
        seq.append("NURSING CARE")

    # Terminal discharge only for completers (otherwise right-censored).
    if is_completer:
        seq.append(TERMINAL)

    return seq


def main() -> None:
    rng = np.random.default_rng(SEED)

    names = list(COHORTS)
    weights = np.array([COHORTS[c][3] for c in names])
    weights = weights / weights.sum()
    patient_cohorts = rng.choice(names, size=N_PATIENTS, p=weights)

    # A small minority of trajectories reach a terminal discharge, reproducing
    # the heavy right-censoring reported in the study (about 1-4%).
    completers = set(rng.choice(range(1, N_PATIENTS + 1), size=2, replace=False))

    rows = []
    for pid in range(1, N_PATIENTS + 1):
        cohort = str(patient_cohorts[pid - 1])
        code, desc, cohort_a, _ = COHORTS[cohort]
        seq = build_trace(rng, pid in completers)

        start = pd.Timestamp("2020-10-01", tz="UTC") + pd.Timedelta(
            days=int(rng.integers(0, 540)))
        ts = start
        quarter_caso = _quarter(start)

        for i, cat_model in enumerate(seq):
            if i > 0:
                gap = 21 if cat_model == "SYSTEMIC TREATMENT" else \
                    float(rng.exponential(18) + 1)
                ts = ts + pd.Timedelta(days=float(gap),
                                       hours=float(rng.integers(8, 18)))

            remote = (cat_model in TELE_VARIANT
                      and rng.random() < TARGET_TELE_SHARE / 0.75)

            base_label = str(rng.choice(ACTIVITY_LABEL[cat_model]))
            if remote:
                category = TELE_VARIANT[cat_model]
                activity = f"Telemedicine {base_label}"
                center = "TELEMEDICINE"
            else:
                category = cat_model
                activity = base_label
                center = str(rng.choice(CENTERS[:2] + CENTERS[3:]))

            rows.append({
                "ANON.ID": pid,
                "ACTIVITY": activity,
                "ACTIVITY_MODEL": base_label,
                "CATEGORY": category,
                "CATEGORY_MODEL": cat_model,
                "FEC_CREACION": ts.strftime("%Y-%m-%d %H:%M:%S+00:00"),
                "SUBSISTEMA": str(rng.choice(["FCE2", "ISUC", "ALERT"],
                                             p=[0.7, 0.2, 0.1])),
                "MEDICO": f"PHYSICIAN_{int(rng.integers(1, 60)):03d}",
                "ESPECIALIDAD": SPECIALTY[cat_model],
                "COD_CIE_10": code,
                "COD_CIE_10_clean": code,
                "DES_CIE_10": desc,
                "CENTRO": center,
                "cohorte_A": cohort_a,
                "cohorte_C": cohort,
                "quarter_caso": quarter_caso,
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUTFILE, sep=";", index=False, encoding="utf-8-sig")

    tele = df["CATEGORY"].str.contains("TELEMEDIC", case=False, na=False)
    model_cats = set(MANDATORY_HEAD) | set(TREATMENT_CHOICE) | \
        set(OPTIONAL_SUPPORT) | set(FOLLOW_UP) | {TERMINAL}
    oom = ~df["CATEGORY_MODEL"].isin(model_cats)
    last = df.groupby("ANON.ID")["CATEGORY_MODEL"].last()
    print(f"Wrote {OUTFILE}")
    print(f"  patients ............ {df['ANON.ID'].nunique()}")
    print(f"  events .............. {len(df)}")
    print(f"  telemedicine share .. {100 * tele.mean():.1f}%")
    print(f"  out-of-model share .. {100 * oom.mean():.1f}%")
    print(f"  reach discharge ..... {(last == TERMINAL).sum()}/{df['ANON.ID'].nunique()}")
    print(f"  events per patient .. {df.groupby('ANON.ID').size().min()}"
          f"-{df.groupby('ANON.ID').size().max()}")
    print("  cohort sizes ........ "
          + ", ".join(f"{k}:{int(v)}" for k, v in
                      df.groupby('ANON.ID')['cohorte_C'].first().value_counts().items()))


if __name__ == "__main__":
    main()
