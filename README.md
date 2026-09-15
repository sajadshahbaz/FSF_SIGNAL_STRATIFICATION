# FSF_SIGNAL_STRATIFICATION

# Feature Signal Framework (FSF)

FSF (Feature Signal Framework) is a signal-aware feature stratification framework designed to characterize feature behavior under repeated perturbations by estimating probabilistic response states and deriving interpretable stability classes.

Rather than focusing solely on statistical significance or effect magnitude, FSF quantifies how consistently individual features exhibit specific response patterns across perturbation experiments. The framework summarizes these behaviors using probabilistic signal states and the Signal Stratification Index (SSI), enabling biologically interpretable signal classification.

---

# Conceptual Overview

FSF operates in four major stages:

1. **Perturbation**

   - Resampling
   - Noise injection
   - Threshold variation

2. **State Probability Estimation**

   Estimation of

   - P(Up)
   - P(Down)
   - P(Constant)

3. **Signal Quantification**

   - Signal Stratification Index (SSI)

4. **Signal Classification**

   - Highly stable Up
   - Highly stable Down
   - Highly stable Constant
   - Transitional
   - Low Stability

The resulting feature-level signal classes can subsequently be aggregated to reveal condition-level signal architectures and biological response patterns.

---

# Repository Structure

```
config/
docs/
manuscript/
results/
scripts/
```

## Directory Description

### config/

Project configuration files.

### docs/

Project documentation, conceptual notes, manuscript planning, benchmark specifications, and repository documentation.

### manuscript/

Working manuscript files (not intended as a permanent publication archive).

### results/

Publication-ready figures, supplementary tables, benchmark outputs, and processed manuscript results.

### scripts/

Complete computational workflow including

- preprocessing
- state assignment
- SSI calculation
- synthetic benchmarking
- biological analyses
- annotation integration
- enrichment analysis
- manuscript figure generation

---

# Main Figures

## Figure 1

Overview of the FSF workflow.

## Figure 2

Probability simplex partitioned into FSF stability regions.

## Figure 3

Condition-level signal architecture and signal-class composition.

## Figure 4

SSI distributions across biological conditions.

## Figure 5

Biological interpretation of signal classes.

## Figure 6

Synthetic benchmark validation.

## Figure 7

Conceptual interpretation model.

---

# Software Requirements

The analysis pipeline was developed using

- R (≥ 4.5)
- ggplot2
- dplyr
- tidyr
- readr
- data.table
- patchwork
- cowplot

Publication-quality PNG figures are generated from vector PDF outputs using **Poppler**.

## Install Poppler

### Ubuntu

```bash
sudo apt install poppler-utils
```

### macOS

```bash
brew install poppler
```

### Conda

```bash
conda install -c conda-forge poppler
```

---

# Reproducing Figures

## Generate Figure 2

```bash
Rscript scripts/manuscript/make_figure2_simplex.R
```

The script produces the publication PDF

```
results/manuscript/figures/main/Figure2_FSF_probability_simplex_regions.pdf
```

Generate the publication PNG from the PDF

```bash
pdftoppm \
-png \
-r 300 \
results/manuscript/figures/main/Figure2_FSF_probability_simplex_regions.pdf \
results/manuscript/figures/main/Figure2_FSF_probability_simplex_regions

mv \
results/manuscript/figures/main/Figure2_FSF_probability_simplex_regions-1.png \
results/manuscript/figures/main/Figure2_FSF_probability_simplex_regions.png
```

---

# Figure Generation Policy

The canonical publication figures are generated as **vector PDF** files.

PNG figures are intentionally produced from the PDFs using **Poppler (`pdftoppm`)** rather than directly from R graphics devices. This ensures

- consistent font rendering
- identical appearance across operating systems
- publication-quality rasterization
- avoidance of Cairo/font rendering issues

The PDF files should always be regarded as the authoritative figure outputs.

---

# Synthetic Benchmark

Generate the synthetic benchmark

```bash
bash run_fsf_v1_synthetic_benchmark.sh
```

---

# Current Project Status

Current repository status

**FSF manuscript preparation**

This repository accompanies the development of the Feature Signal Framework (FSF) manuscript.

---

# Citation

Citation information will be updated following manuscript acceptance and publication.

---

# Author

**Dr. Sajad Shahbazi**

Department of Animal Physiology and Development  
Faculty of Biology  
Adam Mickiewicz University  
Poznań, Poland
