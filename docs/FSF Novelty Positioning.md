# FSF Novelty Positioning (Version 1.0)

## Why FSF Exists

Most transcriptomic analysis frameworks focus on identifying important features.

Examples include:

* Differential expression analysis
* Feature selection
* Machine learning classification
* Biomarker prioritization
* Network analysis

These approaches primarily answer:

> Which features are important?

FSF addresses a different question:

> What type of signal does each feature represent under repeated perturbation?

This distinction forms the conceptual basis of the framework.

---

## Core Concept: Signal Identity

FSF introduces the concept of signal identity.

Signal identity describes the reproducible behavioral pattern exhibited by a feature across repeated perturbations of a dataset.

Rather than focusing solely on significance, effect size, or classification importance, FSF characterizes how consistently a signal behaves.

A feature may exhibit one of four identities:

* Stable Up
* Stable Down
* Stable Constant
* Unstable

These identities are independent of statistical significance.

---

## Signal Identity vs Differential Expression

Differential expression methods determine whether a feature differs between groups.

FSF evaluates whether the observed behavior remains reproducible under perturbation.

Two equally significant genes may possess different signal identities:

Gene A:
Stable Up

Gene B:
Unstable

Although both may be statistically significant, they represent fundamentally different biological behaviors.

---

## Signal Identity vs Feature Importance

Feature-selection methods identify features that contribute to prediction or discrimination.

FSF does not evaluate predictive importance.

FSF evaluates behavioral reproducibility.

A highly predictive feature may still be unstable.

A modestly predictive feature may exhibit highly stable behavior.

---

## Signal Identity vs Instability

Instability quantifies inconsistency.

FSF classifies behavioral identity.

Instability may contribute to interpretation of the Unstable class, but instability and signal identity are not equivalent concepts.

Signal identity is categorical.

Instability is quantitative.

---

## Signal Stratification

FSF proposes that biological signals can be stratified according to reproducible behavioral classes.

This process is termed signal stratification.

The framework therefore moves beyond signal detection and toward signal characterization.

---

## Future Direction: Signal Identity Transitions

A potential extension of FSF is the study of signal identity transitions.

Examples include:

Stable Up → Unstable

Stable Constant → Stable Up

Stable Down → Unstable

These transitions may reveal how biological signals respond to increasing heterogeneity, subgroup conflict, temporal divergence, or perturbation intensity.

Signal identity transitions are not part of FSF Version 1.0 but represent a future methodological direction.

---

## Central Claim

FSF does not identify whether a signal exists.

FSF identifies the behavioral identity of an already detected signal.

The framework therefore complements existing analytical methods rather than replacing them.

