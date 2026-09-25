# Feature Signal Stratification (FSF)

## Current FSF v1 authority

The canonical current FSF v1 results are under `results/current_fsf_v1/`, with
current manuscript products under `results/current_fsf_v1/manuscript/`. See
[`docs/CURRENT_FSF_AUTHORITY.md`](docs/CURRENT_FSF_AUTHORITY.md) for the governing
authority contract. Historical and legacy material is described by
[`docs/FSF_LEGACY_PROVENANCE.md`](docs/FSF_LEGACY_PROVENANCE.md); it is
provenance only and is not current FSF authority.


## Reproducible computational release

The validated computational release is **FSF 0.1.0** at commit:

```text
247e038cc2c206886b03058b1845b42023076033
```

This commit is the frozen computational release underlying the current FSF v1
manuscript analyses and results. Later documentation-only commits do not define
a new computational or scientific release.


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
R CMD INSTALL FSF_0.1.0.tar.gz
```

During repository development, `devtools::install()` may also be used if that
development tool is already available.

### Installation from GitHub

Install the exact frozen commit with `remotes`:

```r
remotes::install_github(
  "sajadshahbaz/FSF_SIGNAL_STRATIFICATION",
  ref = "247e038cc2c206886b03058b1845b42023076033"
)
```

Or clone, check out, build, and install that commit with standard tools:

```sh
git clone https://github.com/sajadshahbaz/FSF_SIGNAL_STRATIFICATION.git
cd FSF_SIGNAL_STRATIFICATION
git checkout 247e038cc2c206886b03058b1845b42023076033
R CMD build .
R CMD INSTALL FSF_0.1.0.tar.gz
```

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
library(FSF)

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

`fsf_analyze()` returns one row per Signal Identity. With `group_cols`, those
columns appear first, followed by `feature_id`, `n_perturbations`, `n_up`,
`n_down`, `n_const`, `p_up`, `p_down`, `p_const`, `dominant_state`, `ssi`,
`stability_region`, `signal_class`, and `stability_deviation`. The three state
probabilities describe the observed Up, Down, and Constant proportions; SSI is
their maximum; `dominant_state` records the maximizing state (or `tied`);
`stability_region` records the SSI interval; and `signal_class` combines region
and direction only above the Low-Stability boundary.

`fsf_architecture()` returns all ten current signal classes for every condition,
including zero-count classes, as the condition column, `signal_class`,
`class_count`, and `class_proportion`.

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
research provenance, validation workflows, scientific authorities, and
manuscript-generation material. Historical payload identities are retained in
public checksum manifests. Those repository resources
are not required for normal installed-package use and are excluded from the
package source tarball.

## Reproducing the FSF v1 manuscript analysis

Package reproducibility and manuscript reproducibility have different scopes:

- **Package reproducibility:** the dataset-independent R package can be built,
  installed, and used without the manuscript data or workflows.
- **Manuscript reproducibility:** the governed workflows under
  `scripts/current/`, current authorities under `results/current_fsf_v1/`, and
  the external artifacts below are required. Material under archive or legacy
  locations is provenance only and is not a current scientific input.

Run current workflow entrypoints from the repository root. Set portable
locators rather than embedding a machine-specific path:

```sh
export FSF_EXTERNAL_ARTIFACT_ROOT=/path/to/fsf-external-artifacts
export FSF_ANNOTATION_MASTER=/path/to/FSF_v1_annotation_master.tsv

Rscript scripts/current/build_current_external_artifact_manifests.R
Rscript scripts/current/audit_current_manuscript_bundle.R \
  results/current_fsf_v1/manuscript
Rscript scripts/current/figures/build_current_manuscript_figures.R

FSF_REPO_ROOT="$PWD" \
FSF_ANNOTATION_MASTER="$FSF_ANNOTATION_MASTER" \
python3 scripts/current/build_current_biological_validation.py --build

FSF_REPO_ROOT="$PWD" \
FSF_ANNOTATION_MASTER="$FSF_ANNOTATION_MASTER" \
python3 scripts/current/build_current_biological_validation.py --verify
```

The manifest builder validates locked sizes, SHA-256 values, and schemas before
writing the governed manifests. The manuscript audit requires the manuscript
bundle root as its sole argument. The full figure build requires the frozen
annotation master; focused figure modes and their requirements are documented
by the figure builder's usage check. Biological validation uses the frozen FSF
feature authority plus the annotation master, and `--verify` rebuilds in an
isolated temporary location to detect output drift. These are the existing
entrypoints, not a replacement pipeline; see
[`docs/CURRENT_FSF_AUTHORITY.md`](docs/CURRENT_FSF_AUTHORITY.md) and the current
workflow scripts for their full contracts.

### External artifacts

The following current authorities are intentionally not stored as Git blobs.
The first three are relative to `FSF_EXTERNAL_ARTIFACT_ROOT`; the annotation
master is supplied directly through `FSF_ANNOTATION_MASTER`.

| Artifact | Purpose | Expected SHA-256 | How FSF locates it | Included in Git? |
|---|---|---|---|---|
| `current_fsf_v1/manuscript/enrichment/go/go_enrichment_all.tsv` | Complete GO enrichment result authority | `1b991db9a178bead03adf3fd96246d9b8e4b08e40828af300d29aca75b260411` | Relative to `FSF_EXTERNAL_ARTIFACT_ROOT` | No |
| `current_fsf_v1/manuscript/enrichment/kegg/kegg_enrichment_all.tsv` | Complete combined KO/PATHWAY enrichment result authority | `bcf428a508975b0ffb81324852243adcf8cc24f136de5ff63bb266a10150fa32` | Relative to `FSF_EXTERNAL_ARTIFACT_ROOT` | No |
| `current_fsf_v1/representative_features/current_fsf_stable_signal_annotated_catalog.tsv` | Stable/Highly Stable candidate catalog joined to biological annotation | `2002148dcaa38a1b4e84c173d78e8aa34495a2dbb98f0fe5738b1ef0445fabdc` | Relative to `FSF_EXTERNAL_ARTIFACT_ROOT` | No |
| `FSF_v1_annotation_master.tsv` | Frozen feature-level biological annotation authority | `25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584` | Exact path in `FSF_ANNOTATION_MASTER` | No |

The repository does not currently document a public download location for
these large artifacts. Their expected identities are governed by
`results/current_fsf_v1/manuscript/audit/external_artifacts.tsv`,
`results/current_fsf_v1/representative_features/LARGE_ARTIFACT_MANIFEST.tsv`,
and [`docs/FSF_V1_ANNOTATION_AUTHORITY.md`](docs/FSF_V1_ANNOTATION_AUTHORITY.md).
No DOI, release asset, or public archive locator is claimed until one is
assigned and verified.

### Manuscript output locations

- `results/current_fsf_v1/` contains the current numerical FSF authorities and
  governed derived products.
- `results/current_fsf_v1/manuscript/` is the current manuscript authority
  root.
- `results/current_fsf_v1/manuscript/source_data/` contains frozen figure source
  tables.
- `results/current_fsf_v1/manuscript/figures/` contains the figure manifest,
  audit, and authoritative main and supplementary figure outputs.
- `results/current_fsf_v1/manuscript/biological_validation/` contains the
  deterministic candidate registry, explicit literature/author decisions,
  validated outputs, summary, and manifest.

### Release validation

The exact frozen computational commit above passed these release gates:

| Check | Result |
|---|---|
| `R CMD build` | PASS |
| `R CMD check` | PASS: 0 errors, 0 warnings, 0 notes |
| Fresh isolated source installation | PASS |
| Fresh-session `library(FSF)` | PASS |
| Namespace loading and public API visibility | PASS |
| Installed examples | PASS |
| Installed test suite | PASS: 123 passed, 0 failed, 0 warnings, 0 skipped |
| `FSF_0_TO_100_RELEASE_GATE` | PASS |

This validation applies to
`247e038cc2c206886b03058b1845b42023076033`; a subsequent README-only commit
does not change that frozen computational identity.

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
