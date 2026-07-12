#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Robustness and performance analyses requested during peer review.

Each function reproduces one result type reported in the manuscript and its
response letter: delivery-modality effect sizes, exposure-adjusted deviation
rates, patient-level independence tests, robust transition times, regression
and survival models, and a sensitivity analysis of the censoring/deviation
decomposition.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from pathway_model import MODEL_ACTIVITIES, COHORT_EN

OUT = "outputs"


def cramers_v(table):
    chi2 = stats.chi2_contingency(table)[0]
    n = table.sum()
    return float(np.sqrt(chi2 / (n * (min(table.shape) - 1)))) if n else np.nan


def modality_analysis(df):
    """Whether structural deviations concentrate in the remote modality."""
    df = df.copy()
    df["tele"] = df["CATEGORY"].str.contains("TELEMEDIC", case=False, na=False)
    df["rep_idx"] = df.groupby(["ANON.ID", "CATEGORY_MODEL"]).cumcount()
    df["repeated"] = df["rep_idx"] > 0

    share = df.groupby("CATEGORY_MODEL")["tele"].agg(["sum", "count"])
    share["pct_tele"] = (100 * share["sum"] / share["count"]).round(1)
    share.rename(columns={"sum": "n_tele", "count": "n_total"}, inplace=True)
    share.sort_values("n_total", ascending=False).to_csv(
        f"{OUT}/modality_by_activity.csv")

    table = pd.crosstab(df["repeated"], df["tele"]).values
    pd.DataFrame([dict(test="repetition_x_modality",
                       cramers_v=round(cramers_v(table), 3),
                       tele_event_share=round(100 * df["tele"].mean(), 1))]
                 ).to_csv(f"{OUT}/modality_association.csv", index=False)


def exposure_analysis(df):
    """Out-of-model activity rate per 100 patient-days and a Poisson IRR."""
    df = df.copy()
    df["oom"] = ~df["CATEGORY_MODEL"].isin(MODEL_ACTIVITIES)
    g = df.groupby("ANON.ID").agg(
        n_events=("CATEGORY_MODEL", "count"),
        n_oom=("oom", "sum"),
        t0=("ts", "min"), t1=("ts", "max"))
    g["person_days"] = (g["t1"] - g["t0"]).dt.total_seconds() / 86400
    g["has_oom"] = g["n_oom"] > 0
    g = g[g["person_days"] > 0].copy()
    g["oom_rate_100pd"] = 100 * g["n_oom"] / g["person_days"]

    rho, p_rho = stats.spearmanr(g["person_days"], g["oom_rate_100pd"])
    irr = irr_lo = irr_hi = p_irr = np.nan
    try:
        import statsmodels.api as sm
        g["ld"] = np.log(g["person_days"])
        model = sm.GLM(g["n_oom"], sm.add_constant(g["ld"]),
                       family=sm.families.Poisson(), offset=g["ld"]).fit()
        irr = float(np.exp(model.params["ld"]))
        irr_lo, irr_hi = np.exp(model.conf_int().loc["ld"]).tolist()
        p_irr = float(model.pvalues["ld"])
    except Exception as exc:                       # pragma: no cover
        print("  [exposure] statsmodels unavailable:", exc)

    pd.DataFrame([dict(
        mean_days_with_oom=round(g[g.has_oom].person_days.mean(), 1),
        mean_days_without_oom=round(g[~g.has_oom].person_days.mean(), 1),
        global_oom_rate_100pd=round(100 * g.n_oom.sum() / g.person_days.sum(), 3),
        spearman_rho_rate=round(rho, 3), spearman_p=round(p_rho, 3),
        poisson_irr_duration=round(irr, 3) if np.isfinite(irr) else np.nan,
        irr_ci_low=round(irr_lo, 3) if np.isfinite(irr_lo) else np.nan,
        irr_ci_high=round(irr_hi, 3) if np.isfinite(irr_hi) else np.nan,
        irr_p=round(p_irr, 3) if np.isfinite(p_irr) else np.nan)]
    ).to_csv(f"{OUT}/exposure_adjusted.csv", index=False)


def patient_level_tests(df):
    """Cohort-vs-rest activity-presence tests at the patient level (Cramer's V)."""
    presence = (df.groupby(["ANON.ID", "CATEGORY_MODEL"]).size()
                .unstack(fill_value=0) > 0).astype(int)
    cohort = df.groupby("ANON.ID")["cohorte_C"].first()
    acts = list(presence.columns)
    rows = []
    for c in cohort.unique():
        inside = cohort == c
        if inside.sum() < 2 or (~inside).sum() < 2:
            continue
        ids_in = inside.index[inside]
        ids_out = inside.index[~inside]
        T = np.vstack([presence.loc[ids_in, acts].sum().values,
                       presence.loc[ids_out, acts].sum().values])
        T = T[:, T.sum(0) > 0]
        if T.shape[1] < 2:
            continue
        chi2, p, _, _ = stats.chi2_contingency(T)
        rows.append(dict(cohort=COHORT_EN.get(c, c), n_patients=int(inside.sum()),
                         chi2=round(chi2, 2), p=float(f"{p:.3e}"),
                         cramers_v=round(cramers_v(T), 3)))
    res = pd.DataFrame(rows)
    if len(res):
        try:
            from statsmodels.stats.multitest import multipletests
            res["p_fdr"] = multipletests(res["p"], method="fdr_bh")[1]
            res["sig_fdr"] = res["p_fdr"] < 0.05
        except Exception:
            pass
    res.to_csv(f"{OUT}/patient_level_tests.csv", index=False)


def transition_times(df):
    """Median and IQR of transition times alongside the (outlier-sensitive) mean."""
    d = df.sort_values(["ANON.ID", "ts"]).copy()
    d["next"] = d.groupby("ANON.ID")["CATEGORY_MODEL"].shift(-1)
    d["nts"] = d.groupby("ANON.ID")["ts"].shift(-1)
    d["delta_days"] = (d["nts"] - d["ts"]).dt.total_seconds() / 86400
    tr = d.dropna(subset=["next", "delta_days"])
    agg = tr.groupby(["CATEGORY_MODEL", "next"])["delta_days"].agg(
        n="count", mean="mean", median="median",
        q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75)).reset_index()
    agg = agg[agg["n"] >= 2].copy()
    agg["iqr"] = (agg["q3"] - agg["q1"]).round(1)
    for c in ["mean", "median", "q1", "q3"]:
        agg[c] = agg[c].round(1)
    agg.rename(columns={"CATEGORY_MODEL": "from", "next": "to"}, inplace=True)
    agg.sort_values("mean", ascending=False).to_csv(
        f"{OUT}/transition_times.csv", index=False)


def regression_survival(df, conformance):
    """Regression of conformance and duration, plus time-to-completion (KM)."""
    df = df.copy()
    df["tele"] = df["CATEGORY"].str.contains("TELEMEDIC", case=False, na=False)
    g = df.groupby("ANON.ID").agg(
        n_events=("CATEGORY_MODEL", "count"),
        n_tele=("tele", "sum"),
        t0=("ts", "min"), t1=("ts", "max")).reset_index()
    g["person_days"] = (g["t1"] - g["t0"]).dt.total_seconds() / 86400
    g["tele_share"] = g["n_tele"] / g["n_events"]
    last = df.sort_values(["ANON.ID", "ts"]).groupby("ANON.ID").tail(1)
    g["completion"] = g["ANON.ID"].map(
        last.set_index("ANON.ID")["CATEGORY_MODEL"].eq("HOSPITAL DISCHARGE").astype(int))

    d = g.merge(conformance[["case", "std_fitness"]],
                left_on="ANON.ID", right_on="case", how="inner")
    d = d[d["person_days"] > 0].copy()
    d["log_events"] = np.log(d["n_events"])

    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
        m1 = smf.ols("std_fitness ~ log_events + tele_share", data=d).fit()
        pd.DataFrame({"coef": m1.params, "ci_low": m1.conf_int()[0],
                      "ci_high": m1.conf_int()[1], "p": m1.pvalues}
                     ).round(4).to_csv(f"{OUT}/regression_conformance.csv")
        m2 = smf.glm("person_days ~ tele_share + n_events", data=d,
                     family=sm.families.Gamma(sm.families.links.log())).fit()
        pd.DataFrame({"coef": m2.params, "expB": np.exp(m2.params),
                      "p": m2.pvalues}).round(4).to_csv(
            f"{OUT}/regression_duration.csv")
    except Exception as exc:                       # pragma: no cover
        print("  [regression] statsmodels unavailable:", exc)

    try:
        from lifelines import KaplanMeierFitter
        kmf = KaplanMeierFitter()
        kmf.fit(d["person_days"], event_observed=d["completion"],
                label="Pathway completion")
        med = kmf.median_survival_time_
        cum_inc = 1 - kmf.survival_function_.iloc[-1, 0]
        pd.DataFrame([dict(
            completions=int(d.completion.sum()), n=len(d),
            median_time_to_completion=("not reached" if not np.isfinite(med)
                                       else round(float(med), 1)),
            cumulative_incidence_pct=round(100 * cum_inc, 1))]
        ).to_csv(f"{OUT}/survival_completion.csv", index=False)
    except Exception as exc:                       # pragma: no cover
        print("  [survival] lifelines unavailable:", exc)


def sensitivity_analysis(conformance, n_boot=2000):
    """Robustness of the decomposition to how the censored tail is defined.

    A parameter alpha in [0, 1] reclassifies a fraction of interior model-moves
    as censoring: alpha = 0 is the strict baseline, alpha = 1 the worst case
    (every model-move treated as censoring). Insertions (log-moves) are events
    that did occur and can never be censoring, so genuine deviation dominates if
    it still dominates at alpha = 1.
    """
    rng = np.random.default_rng(20260710)
    alphas = np.round(np.linspace(0, 1, 11), 3)
    logmove = conformance["logmove"].to_numpy(float)
    interior = conformance["interior_mm"].to_numpy(float)
    tail = conformance["tail_mm"].to_numpy(float)
    n = len(conformance)

    rows, brows = [], []
    for a in alphas:
        genuine = logmove + (1 - a) * interior
        censoring = tail + a * interior
        total = genuine.sum() + censoring.sum()
        conform = ((logmove == 0) & ((1 - a) * interior < 1e-9)).astype(int)
        rows.append(dict(
            alpha=a,
            pct_genuine=round(100 * genuine.sum() / total, 2) if total else np.nan,
            pct_censoring=round(100 * censoring.sum() / total, 2) if total else np.nan,
            pct_prefix_conform=round(100 * conform.mean(), 2)))
    for a, tag in [(0.0, "baseline"), (1.0, "worst_case")]:
        genuine = logmove + (1 - a) * interior
        censoring = tail + a * interior
        pg = np.empty(n_boot)
        for b in range(n_boot):
            idx = rng.integers(0, n, n)
            num = genuine[idx].sum()
            den = num + censoring[idx].sum()
            pg[b] = 100 * num / den if den else np.nan
        total = genuine.sum() + censoring.sum()
        brows.append(dict(
            scenario=tag, alpha=a,
            pct_genuine=round(100 * genuine.sum() / total, 2) if total else np.nan,
            ci_low=round(np.nanpercentile(pg, 2.5), 2),
            ci_high=round(np.nanpercentile(pg, 97.5), 2)))
    sweep = pd.DataFrame(rows)
    sweep.to_csv(f"{OUT}/sensitivity_sweep.csv", index=False)
    pd.DataFrame(brows).to_csv(f"{OUT}/sensitivity_bootstrap.csv", index=False)
    _sensitivity_figure(sweep)


def _sensitivity_figure(sweep):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(sweep.alpha, sweep.pct_censoring, "-o", ms=4, color="#1b485e")
        ax[0].set(xlabel="Tail generosity α", ylim=(0, 100),
                  ylabel="% deviation attributed to censoring",
                  title="(a) Censoring share vs tail definition")
        ax[1].plot(sweep.alpha, sweep.pct_prefix_conform, "-o", ms=4, color="#568b87")
        ax[1].set(xlabel="Tail generosity α", ylim=(0, 40),
                  ylabel="% patients on a conformant prefix",
                  title="(b) Conformance is not rescued by censoring")
        for a in ax:
            a.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{OUT}/fig_sensitivity.png", dpi=200)
        plt.close(fig)
    except Exception as exc:                       # pragma: no cover
        print("  [figure] matplotlib unavailable:", exc)
