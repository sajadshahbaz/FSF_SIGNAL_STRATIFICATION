# Feature Signal Stratification (FSF)

FSF is a small, base-R package for describing how consistently a feature shows
an Up, Down, or Constant response across repeated perturbations. It converts
signed effects into state probabilities, calculates the Signal Stratification
Index (SSI), assigns current FSF signal classes, and summarizes class
architecture across user-defined conditions.

FSF v1 performs signal stratification. It does **not** estimate effects, test
differential expression, select or rank features, run enrichment, or make
biological interpretations.

## Installation

FSF is not currently claimed to be available from CRAN. Build and install it
from a source checkout with standard R tooling:

```sh
R CMD build .
R CMD INSTALL FSF_0.1.0.9000.tar.gz
```

During repository development, `devtools::install()` may also be used if that
development tool is already available.

## Input

The minimum input is a long-form data frame with one row per observed
feature-perturbation pair:

| Column | Requirement |
|---|---|
| `feature_id` | Nonmissing, nonblank identifier |
| `perturbation_id` | Nonmissing, nonblank identifier |
| `effect` | Finite numeric signed effect |

Duplicate `feature_id + perturbation_id` keys are invalid unless independent
strata are explicitly supplied through `group_cols`.

## Minimal analysis

```r
effects <- data.frame(
  condition = rep(c("control", "treated"), each = 6),
  feature_id = rep(rep(c("f1", "f2"), each = 3), 2),
  perturbation_id = rep(paste0("p", 1:3), 4),
  effect = c(1, 1, 0, 0, 0, 0, 1, 1, 1, -1, -1, 0)
)

classified <- fsf_analyze(
  effects,
  tau = 0.5,
  group_cols = "condition"
)
classified

architecture <- fsf_architecture(classified, condition_col = "condition")
architecture
```

`tau` is the positive magnitude threshold used to assign directional states:

- `effect > tau`: Up
- `effect < -tau`: Down
- `-tau <= effect <= tau`: Constant

Exact `+tau` and `-tau` values are therefore Constant.

## Signal Identity and SSI

For each feature identity, FSF reports the proportions `p_up`, `p_down`, and
`p_const`. These sum to one. SSI is their exact maximum:

```text
SSI = max(p_up, p_down, p_const)
```

The dominant state is the state attaining SSI. Equal maxima produce the
`"tied"` dominant state.

The current stability regions are:

- **Low Stability:** `SSI <= 0.50`
- **Transitional:** `0.50 < SSI < 0.75`
- **Stable:** `0.75 <= SSI < 0.90`
- **Highly Stable:** `SSI >= 0.90`

Low Stability is a single non-directional signal class. A dominant state may
still be recorded for a Low-Stability feature, but it does not create classes
such as “Low Stability Up” or “Low Stability Down”. Directional signal classes
are assigned only when SSI is strictly greater than 0.50.

## Public API

- `fsf_assign_states()` assigns Up, Down, and Constant states.
- `fsf_signal_identity()` calculates state counts and probabilities.
- `fsf_classify()` classifies an already aggregated identity table.
- `fsf_analyze()` performs identity calculation and classification.
- `fsf_architecture()` calculates condition-level class composition.

See the installed help pages, for example `?fsf_analyze`, for full schemas,
validation behavior, and examples.

## Package versus research repository

The installed package contains the reusable, dataset-independent FSF
computational core and package tests. This repository additionally retains
research provenance, validation workflows, historical pipelines, scientific
authorities, and manuscript-generation material. Those repository resources
are not required for normal installed-package use and are excluded from the
package source tarball.

## Citation

No paper DOI or publication citation is asserted here. Until authoritative
publication metadata is available, use R's standard package citation:

```r
citation("FSF")
```

## Author

Sajad Shahbazi
Department of Animal Physiology and Development
Faculty of Biology, Adam Mickiewicz University, Poznań, Poland
