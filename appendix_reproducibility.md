# Appendix A. Reproducibility: pseudocode and synthetic data

## A.1. Data and code availability

The clinical event log analysed in this study contains protected patient
information and cannot be released, and the institutional implementation of the
analysis code is not publicly distributable. To make the methodology fully
auditable notwithstanding these constraints, we provide three complementary
artefacts:

1. the complete pseudocode of the analysis pipeline — preprocessing and
   vocabulary mapping (Algorithm A1), censoring-aware conformance and its
   decomposition (Algorithm A2), exposure-adjusted performance and
   modality analysis (Algorithm A3), and the sensitivity analysis of the
   decomposition (Algorithm A4);

2. a synthetic, de-identified event log that reproduces the *structure* of the
   study data — the activity vocabulary, the telemedicine delivery modality,
   out-of-model activities, trace repetition and right-censoring — while
   containing no real patient records; and

3. an executable reference implementation that runs the entire pipeline on the
   synthetic log and regenerates every result type reported in the paper
   (conformance decomposition, modality effect sizes, exposure-adjusted rates,
   patient-level tests, robust transition times, regression and survival models,
   and the sensitivity analysis).

These materials are available as Supplementary Material and from the authors on
request. The synthetic log carries no real patient data; the numerical values it
produces are illustrative and are not the study's clinical results.

## A.2. Reference model

Conformance is assessed against the guideline-derived reference model, encoded as
an acyclic workflow net. The normative pathway comprises a mandatory core —
specialist medical consultation, oncology consultation, oncology committee, an
exclusive treatment step (systemic treatment, surgery or radiotherapy committee),
nursing follow-up and hospital discharge — together with optional supportive
activities (psychology and non-medical care) that may be omitted without penalty.
Delivery modality (in-person versus telemedicine) is treated as an attribute of
each event rather than as a distinct activity, so telemedicine variants are
mapped to their base activity to align the log vocabulary with the model
(Algorithm A1). The same conformance procedure is applied at two abstraction
levels — the collapsed *category* level (primary) and the finer *activity* level
— by selecting the corresponding activity attribute.

## A.3. Pseudocode

Throughout, a *trace* is the time-ordered sequence of events of one patient, and
an *optimal alignment* of a trace against the reference model is a minimum-cost
sequence of moves, where a *synchronous move* matches an observed event to a
model step, a *log-move* consumes an observed event with no model counterpart
(an insertion), and a *model-move* executes a required step with no observed
counterpart. Silent (routing) moves are cost-free.

```
Algorithm A1  Event-log preprocessing and vocabulary mapping
------------------------------------------------------------
Input : raw event log L_raw
Output: processed log L with model vocabulary and cohort labels

 1: for each event e in L_raw do
 2:     e.time  <- parse e.timestamp as UTC datetime
 3:     e.act   <- normalize case and accents of e.activity
 4:     e.icd   <- clean e.icd10  (remove separators; keep base code)
 5:     if e is a telemedicine variant then           # collapse modality
 6:         e.category_model <- base category of e
 7:         e.activity_model <- base activity of e
 8:     else
 9:         e.category_model <- e.category
10:         e.activity_model <- e.activity
11:     e.telemedicine <- [e was a telemedicine variant]
12: for each patient p do
13:     cohort(p)  <- rule over the dominant ICD-10 family of p
14:                   (C18->Colon, C20->Rectum, C50->Breast,
15:                    C34->Lung, C25->Pancreas, otherwise Others)
16:     quarter(p) <- calendar quarter of the first event of p
17: sort the events of every case by e.time
18: return L
```

```
Algorithm A2  Censoring-aware conformance and decomposition
-----------------------------------------------------------
Input : processed log L at abstraction level l; reference net N
        with initial marking m_i and final marking m_f
Output: per-trace conformance records; aggregate decomposition

 1: for each case c with trace sigma_c do
 2:     gamma <- OptimalAlignment(sigma_c, N, m_i, m_f)
 3:     S  <- number of synchronous moves in gamma
 4:     Lg <- number of log-moves        (insertions)
 5:     M  <- number of visible model-moves
 6:     tail     <- number of trailing visible model-moves in gamma
 7:                 (those after the last synchronous or log-move)
 8:     interior <- M - tail
 9:     genuine(c)   <- Lg + interior      # insertions + skipped interior steps
10:     censoring(c) <- tail               # required steps not yet reached
11:     prefix_fitness(c) <- S / (S + genuine(c))
12:     conformant(c)     <- [genuine(c) = 0]
13: G <- sum_c genuine(c);   C <- sum_c censoring(c)
14: pct_genuine   <- 100 * G / (G + C)
15: pct_censoring <- 100 * C / (G + C)
16: return {prefix_fitness, conformant}, (pct_genuine, pct_censoring)
```

Right-censoring is isolated at the level of the alignment: required steps that
remain only at the *end* of the optimal alignment are attributed to an
incomplete (not yet finished) trajectory, whereas insertions and required steps
skipped in the *interior* of the trace are genuine structural deviations.

```
Algorithm A3  Exposure-adjusted performance and modality analysis
-----------------------------------------------------------------
Input : processed log L; per-case output of Algorithm A2
Output: modality effect size; exposure-adjusted rates; patient-
        level tests; transition times; regression and survival fits

 1: for each case c do
 2:     n_events(c), tele_share(c) <- aggregate over the events of c
 3:     person_days(c) <- (last event time - first event time) of c, in days
 4:     n_oom(c)       <- number of out-of-model events of c
 5: # (i) Modality: are deviations specific to the remote modality?
 6: build contingency table [event repeated] x [telemedicine]
 7: V_modality <- Cramer's V of the table
 8: # (ii) Exposure: deviation rate per 100 patient-days
 9: rate(c) <- 100 * n_oom(c) / person_days(c)
10: fit Poisson GLM: n_oom ~ offset(log person_days) + log person_days
11: IRR_duration <- exp(coefficient of log person_days)
12: # (iii) Patient-level independence, corrected for multiplicity
13: for each cohort k do
14:     chi-square of [activity presence] x [cohort k vs rest]; report V_k
15: adjust p-values by the Benjamini-Hochberg (FDR) procedure
16: # (iv) Robust transition times
17: for each transition a -> b report mean, median and IQR of the delay
18: # (v) Regression and survival
19: OLS:        fitness      ~ log n_events + tele_share
20: Gamma GLM:  person_days  ~ tele_share + n_events        (log link)
21: Kaplan-Meier: time to pathway completion (event = terminal
22:               discharge), treating unfinished trajectories as censored
23: return V_modality, rate, IRR_duration, {V_k, FDR}, transition
24:        statistics, and the regression and survival fits
```

```
Algorithm A4  Sensitivity of the censoring / deviation decomposition
--------------------------------------------------------------------
Input : per-trace counts {S, Lg, interior, tail} from Algorithm A2;
        grid A = {0, 0.1, 0.2, ..., 1}; bootstrap size B
Output: decomposition and prefix-conformance as functions of the
        tail-generosity alpha, with bootstrap confidence intervals

 1: for alpha in A do
 2:     genuine_a(c)   <- Lg(c) + (1 - alpha) * interior(c)
 3:     censoring_a(c) <- tail(c) + alpha * interior(c)
 4:     pct_genuine(alpha) <- 100 * sum_c genuine_a(c)
 5:                            / sum_c (genuine_a(c) + censoring_a(c))
 6:     conformant_a(c)    <- [Lg(c) = 0 and (1 - alpha) * interior(c) = 0]
 7:     pct_conformant(alpha) <- 100 * mean_c conformant_a(c)
 8: for alpha in {0, 1} do
 9:     draw B case-resamples; report 95% CI of pct_genuine(alpha)
10: # alpha = 0: strict baseline; alpha = 1: worst case, in which every
11: # model-move is counted as censoring. Insertions (log-moves) are events
12: # that did occur and can never be censoring, so genuine deviation
13: # dominates whenever it dominates at alpha = 1.
14: return {pct_genuine(alpha), pct_conformant(alpha)} and bootstrap CIs
```

## A.4. Synthetic event log and reference implementation

The reference implementation consists of a generator that writes the synthetic
event log and a pipeline that executes Algorithms A1–A4 on it. The generator
produces fifty patients across the study cohorts, with the same event schema as
the processed study log, and reproduces the structural features that drive the
analysis: a telemedicine share close to that observed, a minority of out-of-model
activities, repeated activities within a case, entry at several different nodes
(not only the specialist consultation), and right-censoring — only a small
minority of trajectories reach a terminal discharge. Running the pipeline on this
log yields, in order, the conformance decomposition, the modality effect size,
the exposure-adjusted rates and incidence-rate ratio, the patient-level tests,
the robust transition times, the regression and survival fits, and the
sensitivity analysis, each written as a separate result table. All random
elements use fixed seeds, so a run is fully reproducible.
