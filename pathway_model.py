#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Guideline reference model and alignment classification.

The normative oncology pathway is built in code as an acyclic Petri net, so no
proprietary BPMN file is required to reproduce the analysis. The same procedure
applies at either abstraction level of the study by pointing the pipeline at the
corresponding column of the log (``CATEGORY_MODEL`` for the category level or
``ACTIVITY_MODEL`` for the finer activity level); this module uses the category
level, which is the primary analytical level of the manuscript.
"""
from __future__ import annotations

from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.petri_net.utils import petri_utils

# Ordered guideline pathway. Mandatory steps must be executed; optional steps
# may be skipped without penalty; the treatment step is an exclusive choice.
REFERENCE_PATHWAY = [
    ("mandatory", ["SPECIALIST MEDICAL CONSULTATION"]),
    ("mandatory", ["ONCOLOGY CONSULTATION"]),
    ("mandatory", ["ONCOLOGY COMMITTEE"]),
    ("optional",  ["PSYCHOLOGY"]),
    ("choice",    ["SYSTEMIC TREATMENT", "SURGERY", "RADIOTHERAPY COMMITTEE"]),
    ("mandatory", ["NURSING CARE"]),
    ("optional",  ["NON-MEDICAL CARE"]),
    ("mandatory", ["HOSPITAL DISCHARGE"]),
]

# Set of activities permitted by the reference model (its allowlist).
MODEL_ACTIVITIES = {a for _, labels in REFERENCE_PATHWAY for a in labels}

COHORT_EN = {"Colon": "Colon", "Mama": "Breast", "Otros": "Others",
             "Pulmon": "Lung", "Pancreas": "Pancreas", "Recto": "Rectum"}


def build_reference_net():
    """Construct the guideline pathway as an acyclic Petri net."""
    net = PetriNet("normative_oncology_pathway")
    counter = [0]

    def add_place(name):
        p = PetriNet.Place(name)
        net.places.add(p)
        return p

    def add_transition(label):
        counter[0] += 1
        t = PetriNet.Transition(f"t{counter[0]}", label)
        net.transitions.add(t)
        return t

    place = add_place("p0")
    im = Marking({place: 1})
    for i, (kind, labels) in enumerate(REFERENCE_PATHWAY):
        nxt = add_place(f"p{i + 1}")
        for label in labels:                     # exclusive alternatives
            t = add_transition(label)
            petri_utils.add_arc_from_to(place, t, net)
            petri_utils.add_arc_from_to(t, nxt, net)
        if kind == "optional":                   # silent skip transition
            tau = add_transition(None)
            petri_utils.add_arc_from_to(place, tau, net)
            petri_utils.add_arc_from_to(tau, nxt, net)
        place = nxt
    fm = Marking({place: 1})
    return net, im, fm


def classify(alignment):
    """Split alignment moves into synchronous / insertion / skipped steps.

    log-move        (a, '>>')    observed activity not permitted -> insertion
    model-move      ('>>', b)    required step not observed
        - trailing (after the last matched event)  -> right-censoring
        - interior                                  -> genuine skipped step
    silent move     ('>>', None) ignored (cost-free)

    prefix fitness      = sync / (sync + genuine deviation)
    prefix conformance  = 1 iff there is no genuine deviation (the observed
                          trace is a valid prefix of a conformant path).
    """
    seq = []  # 'S' sync, 'L' log-move, 'M' visible model-move
    sync = logmove = mm_visible = 0
    for move in alignment:
        a, b = move[0], move[1]
        if b is None:                            # silent transition
            continue
        if a != ">>" and b != ">>":
            sync += 1
            seq.append("S")
        elif a != ">>" and b == ">>":
            logmove += 1
            seq.append("L")
        elif a == ">>" and b != ">>":
            mm_visible += 1
            seq.append("M")

    tail = 0
    for s in reversed(seq):                      # trailing model-moves = censoring
        if s == "M":
            tail += 1
        else:
            break
    interior_mm = mm_visible - tail
    genuine = logmove + interior_mm
    prefix_fitness = sync / (sync + genuine) if (sync + genuine) > 0 else 1.0
    return dict(sync=sync, logmove=logmove, mm_visible=mm_visible,
                tail_mm=tail, interior_mm=interior_mm,
                genuine_dev=genuine, censoring=tail,
                prefix_fitness=prefix_fitness,
                prefix_conform=int(genuine == 0))
