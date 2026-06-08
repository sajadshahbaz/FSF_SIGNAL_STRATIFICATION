# FSF_CORE_DEFINITION.md

# Feature Signal Framework (FSF) Version 1.0

## Core Mathematical Definition

---

## 1. Purpose

Feature Signal Framework (FSF) is a signal stratification framework designed to characterize the reproducible behavioral identity of features under repeated perturbation.

FSF does not perform feature selection, differential expression analysis, classification, prediction, or ranking.

The purpose of FSF is to determine whether a feature exhibits a reproducible signal state across perturbations and to quantify its degree of stability.

---

## 2. Scope of Version 1.0

FSF Version 1.0 contains only:

1. Perturbation Framework
2. Signal State Assignment
3. Signal-State Probabilities
4. Signal Stratification Index (SSI)
5. Stability Deviation
6. Signal Class Assignment

FSF Version 1.0 does not include:

* Feature selection
* Signal-aware ranking
* Switching metrics
* Shifting metrics
* Transition-point detection
* Regime detection
* Adaptive-response detection

These topics are outside the scope of the present version.

---

## 3. Perturbation Framework

Let a feature (i) be evaluated across (N) perturbation experiments.

Each perturbation produces a directional signal assignment:

[
S_{ij} \in {-1,0,+1}
]

where:

[
+1 = Up
]

[
0 = Constant
]

[
-1 = Down
]

and

[
j = 1,\ldots,N
]

represents the perturbation index.

---

## 4. Signal-State Probabilities

For feature (i), the probability of each signal state is calculated as:

### Probability of Up

[
P_i(\text{up})
==============

\frac{#(S_{ij}=+1)}{N}
]

### Probability of Constant

[
P_i(\text{const})
=================

\frac{#(S_{ij}=0)}{N}
]

### Probability of Down

[
P_i(\text{down})
================

\frac{#(S_{ij}=-1)}{N}
]

subject to:

[
P_i(\text{up})
+
P_i(\text{const})
+
P_i(\text{down})
================

1
]

---

## 5. Dominant Signal Identity

The dominant signal identity is defined as the signal state with the highest probability.

[
DSI_i
=====

\operatorname*{argmax}
\left[
P_i(\text{up}),
P_i(\text{const}),
P_i(\text{down})
\right]
]

Possible dominant signal identities:

* Stable Up
* Stable Constant
* Stable Down

---

## 6. Signal Stratification Index (SSI)

The Signal Stratification Index (SSI) measures the strength of the dominant signal identity.

[
SSI_i
=====

\max
\left[
P_i(\text{up}),
P_i(\text{const}),
P_i(\text{down})
\right]
]

Properties:

[
\frac{1}{3}
\le
SSI_i
\le
1
]

Interpretation:

* Higher SSI indicates stronger dominance of a single signal state.
* Lower SSI indicates increased competition among signal states.

Examples:

### Perfect Stable Up

[
P(\text{up})=1
]

[
SSI=1
]

### Perfect Stable Down

[
P(\text{down})=1
]

[
SSI=1
]

### Perfect Stable Constant

[
P(\text{const})=1
]

[
SSI=1
]

### Maximum Signal Ambiguity

[
P(\text{up})
============

# P(\text{down})

# P(\text{const})

\frac13
]

[
SSI
===

\frac13
]

---

## 7. Stability Deviation

Stability Deviation quantifies departure from ideal signal stability.

For feature (i):

[
SD_i
====

## 1

SSI_i
]

Properties:

[
0
\le
SD_i
\le
\frac23
]

Interpretation:

### Perfect Stability

[
SSI=1
]

[
SD=0
]

### Increasing Instability

As SSI decreases,

[
SD
]

increases.

### Maximum Signal Ambiguity

[
SSI=\frac13
]

[
SD=\frac23
]

---

## 8. Signal Stability Classes

FSF Version 1.0 assigns features to stability classes using SSI.

### Highly Stable

[
SSI \ge 0.90
]

### Stable

[
0.75 \le SSI < 0.90
]

### Weakly Stable / Transitional

[
0.60 \le SSI < 0.75
]

### Instable

[
SSI < 0.60
]

Thresholds may be adjusted during benchmarking and sensitivity analysis.

---

## 9. Final FSF Output

For each feature, FSF returns:

1. Dominant Signal Identity (DSI)
2. Signal Stratification Index (SSI)
3. Stability Deviation (SD)
4. Stability Class

Example output:

| Feature | DSI             | SSI  | SD   | Class         |
| ------- | --------------- | ---- | ---- | ------------- |
| Gene_A  | Stable Up       | 0.95 | 0.05 | Highly Stable |
| Gene_B  | Stable Down     | 0.82 | 0.18 | Stable        |
| Gene_C  | Stable Constant | 0.67 | 0.33 | Weakly Stable |
| Gene_D  | Mixed           | 0.45 | 0.55 | Instable      |

---

## 10. Interpretation

FSF is a signal stratification framework.

FSF identifies the dominant reproducible signal identity of a feature and quantifies the degree of departure from ideal stability.

FSF does not infer biological function, causality, differential expression, predictive importance, or feature relevance.

Such downstream interpretations remain the responsibility of subsequent analytical frameworks.

---

## Version

FSF Core Definition Version 1.0

