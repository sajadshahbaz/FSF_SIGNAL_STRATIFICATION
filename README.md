# FSF_SIGNAL_STRATIFICATION

## Feature Signal Framework (FSF)

FSF (Feature Signal Framework) is a signal-aware feature stratification framework that characterizes feature behavior under perturbation by estimating state probabilities and deriving interpretable signal classes.

Rather than focusing solely on differential magnitude or statistical significance, FSF evaluates how consistently a feature exhibits specific response patterns across perturbations.

---

## Conceptual Overview

FSF operates in four stages:

1. **Perturbation**

   * Resampling
   * Noise injection
   * Threshold variation

2. **State Probability Estimation**

   * P(up)
   * P(down)
   * P(const)

3. **Signal Quantification**

   * Signal Stratification Index (SSI)
   * Stability Deviation (SD)

4. **Signal Classification**

   * Stable Up
   * Stable Constant
   * Stable Down
   * Transitional
   * Instable

The resulting signal classes can be interpreted biologically and aggregated to reveal condition-level signal architectures.

---

## Repository Structure

```text
config/
docs/
manuscript/
results/
scripts/
```

### Key Components

#### Documentation

```text
docs/
```

Contains formal definitions, benchmark plans, claims, limitations, and project specifications.

#### Analysis Pipeline

```text
scripts/
```

Contains the complete FSF workflow, including:

* State assignment
* SSI calculation
* Stability deviation estimation
* Synthetic benchmarking
* Real-data analysis
* Annotation integration
* Enrichment analysis
* Figure generation

#### Manuscript Outputs

```text
results/manuscript/
```

Contains publication-ready figures and curated summary tables.

---

## Main Manuscript Figures

### Figure 1

Overview of the FSF framework and signal stratification workflow.

### Figure 2

FSF probability simplex and signal-state regions.

### Figure 3

Condition-level signal architectures and signal-class composition.

### Figure 4

Condition-specific SSI distributions.

### Figure 5

Biological theme distributions and condition profiles.

### Figure 6

Synthetic benchmark validation.

### Figure 7

Final conceptual model.

---

## Reproducing Benchmark Results

Generate synthetic benchmark outputs:

```bash
bash run_fsf_v1_synthetic_benchmark.sh
```

Generate benchmark figures:

```bash
Rscript scripts/36D_generate_benchmark_figures.R
```

Generate Figure 2:

```bash
Rscript scripts/36J_rebuild_Figure2_FSF_simplex.R
```

Generate Figure 3:

```bash
Rscript scripts/36K_rebuild_Figure3_publication.R
```

---

## Current Status

Repository frozen for manuscript preparation.

Version:
**FSF manuscript freeze v1.0**

---

## Citation

Citation information will be updated following manuscript submission and publication.

---

## Author

Dr. Sajad Shahbazi

Adam Mickiewicz University, Poznań, Poland

