# FSF Version 1.0 Locked Definition

## Core Identity

Feature Stability Framework (FSF) Version 1.0 is a signal stratification framework.

FSF does not perform differential expression analysis, feature selection, classification, prediction, clustering, causal inference, or biological mechanism discovery.

FSF operates downstream of external candidate-generation methods and characterizes the reproducible behavioral identity of candidate features under repeated perturbation.

## FSF v1.0 Components

FSF Version 1.0 contains only:

1. Perturbation framework
2. Directional effect estimate input
3. Signal state assignment
4. Signal-state probabilities
5. Signal Stratification Index (SSI)
6. Stability Deviation
7. Signal class assignment

FSF Version 1.0 does not include:

- feature selection
- signal-aware ranking
- switching metrics
- shifting metrics
- transition-point detection
- regime detection
- adaptive-response detection

## Effect Estimate

For each feature i and perturbation p, FSF requires a directional effect estimate:

E(i,p)

The default implementation uses log fold change:

E(i,p) = logFC(i,p)

Alternative directional estimates may be used if they are applied consistently across all perturbations.

## Signal State Assignment

Given a directional threshold tau:

UP:
E(i,p) > tau

DOWN:
E(i,p) < -tau

CONST:
-tau <= E(i,p) <= tau

Default:

tau = 0.5

## Signal-State Probabilities

For each feature i across N perturbations:

Pup(i) = number of UP assignments / N

Pdown(i) = number of DOWN assignments / N

Pconst(i) = number of CONST assignments / N

Constraint:

Pup(i) + Pdown(i) + Pconst(i) = 1

## Dominant Signal Identity

The dominant signal identity is the state with the highest probability:

DSI(i) = argmax[Pup(i), Pdown(i), Pconst(i)]

Possible dominant identities:

- Stable Up
- Stable Down
- Stable Constant

## Signal Stratification Index

SSI measures the strength of the dominant signal identity:

SSI(i) = max[Pup(i), Pdown(i), Pconst(i)]

Interpretation:

- SSI = 1 indicates perfect signal stability.
- Lower SSI indicates increasing competition among signal states.
- SSI is not a p-value.
- SSI is not a feature-importance score.
- SSI is not an instability score.

## Stability Deviation

Stability Deviation quantifies departure from ideal stability:

SD(i) = 1 - SSI(i)

Interpretation:

- SD = 0 indicates perfect stability.
- Higher SD indicates greater departure from a dominant stable signal identity.
- SD is a continuous stability-deviation parameter.

## Stability Classes

FSF v1.0 assigns features using SSI:

Highly Stable:
SSI >= 0.90

Stable:
0.75 <= SSI < 0.90

Weakly Stable / Transitional:
0.60 <= SSI < 0.75

Instable:
SSI < 0.60

## Final FSF Signal Classes

The final FSF v1.0 signal classes are:

- Highly Stable Up
- Stable Up
- Weakly Stable / Transitional Up
- Highly Stable Down
- Stable Down
- Weakly Stable / Transitional Down
- Highly Stable Constant
- Stable Constant
- Weakly Stable / Transitional Constant
- Instable

## Terminology Lock

In FSF Version 1.0, the term Instable refers to features that do not maintain a sufficiently dominant signal identity under perturbation.

The term Unstable is avoided in FSF v1.0 to prevent confusion with broader instability concepts used in separate methodological work.

## Core Output

For every feature, FSF returns:

- feature_id
- Pup
- Pdown
- Pconst
- dominant signal identity
- SSI
- Stability Deviation
- stability level
- final FSF signal class

## Core Claim

FSF identifies the behavioral identity of already detected candidate signals.

FSF does not claim to discover statistically significant features or replace existing analytical methods.
