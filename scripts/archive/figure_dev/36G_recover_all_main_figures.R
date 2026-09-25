#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"
OUT  <- file.path(ROOT, "results/archive/pre_repair_manuscript/figures/main")
LOG  <- file.path(ROOT, "results/archive/pre_repair_manuscript/logs/36G_recover_all_main_figures.log")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sink(LOG, split = TRUE)

cat("Step 36G: recover full clean main figure set\n")
cat("Started:", as.character(Sys.time()), "\n\n")

savep <- function(p, name, w=8, h=5){
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width=w, height=h, device=cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width=w, height=h, dpi=300)
  cat("saved", name, "\n")
}

theme_pub <- function(){
  theme_bw(base_size=12) +
    theme(
      plot.title=element_text(face="bold", hjust=0),
      panel.grid.minor=element_blank(),
      axis.text.x=element_text(angle=35, hjust=1),
      legend.title=element_text(face="bold")
    )
}

bio <- read_tsv(
  file.path(ROOT, "results/real_data/final_manuscript_tables/Table_SignalClass_Biology_CURATED.tsv"),
  show_col_types=FALSE
)

arch <- read_tsv(
  file.path(ROOT, "results/real_data/final_manuscript_tables/Table_Condition_SignalArchitecture_CURATED.tsv"),
  show_col_types=FALSE
)

# ---------------- Figure 1 ----------------

workflow <- tibble(
  x=1:6, y=1,
  label=c(
    "Input matrix\nor candidate set",
    "Repeated\nperturbation",
    "State probabilities\nP(up), P(down), P(const)",
    "Signal Stratification\nIndex (SSI)",
    "Signal class\nassignment",
    "Biological\ninterpretation"
  )
)

p1 <- ggplot(workflow, aes(x,y)) +
  geom_segment(
    data=workflow |> filter(x<6),
    aes(x=x+0.22, xend=x+0.78, y=y, yend=y),
    arrow=arrow(length=unit(0.18,"cm")),
    linewidth=0.45
  ) +
  geom_label(aes(label=label), size=3.2, label.size=0.25, fill="white") +
  scale_x_continuous(limits=c(0.55,6.45), breaks=NULL) +
  scale_y_continuous(limits=c(0.92,1.08), breaks=NULL) +
  labs(title="Figure 1. FSF workflow") +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0))

savep(p1, "Figure1_Workflow", 12, 2.6)

# ---------------- Figure 2 ----------------

simplex <- expand.grid(
  p_up=seq(0,1,by=0.025),
  p_down=seq(0,1,by=0.025)
) |>
  mutate(p_const=1-p_up-p_down) |>
  filter(p_const >= -1e-9) |>
  mutate(
    max_p=pmax(p_up,p_down,p_const),
    region=case_when(
      p_const==max_p ~ "Stable-constant region",
      p_down==max_p ~ "Stable-down region",
      p_up==max_p ~ "Stable-up region",
      TRUE ~ "Boundary"
    ),
    x=p_down+0.5*p_const,
    y=p_const*sqrt(3)/2
  )

p2 <- ggplot(simplex, aes(x,y,fill=region)) +
  geom_point(shape=22, size=2.2, color="white", stroke=0.15) +
  annotate("segment", x=0, xend=1, y=0, yend=0, linewidth=0.4) +
  annotate("segment", x=0, xend=0.5, y=0, yend=sqrt(3)/2, linewidth=0.4) +
  annotate("segment", x=1, xend=0.5, y=0, yend=sqrt(3)/2, linewidth=0.4) +
  annotate("text", x=0.08, y=0.035, label="P(up)", size=4, hjust=0) +
  annotate("text", x=0.92, y=0.035, label="P(down)", size=4, hjust=1) +
  annotate("text", x=0.5, y=0.82, label="P(const)", size=4, vjust=1) +
  scale_x_continuous(limits=c(-0.05,1.05), expand=expansion(mult=0)) +
  scale_y_continuous(limits=c(-0.04,0.92), expand=expansion(mult=0)) +
  coord_equal() +
  labs(title="Figure 2. FSF probability space", x=NULL, y=NULL, fill="Dominant region") +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0), legend.position="right")

savep(p2, "Figure2_Probability_Simplex", 7.5, 6)

# ---------------- Figure 3 ----------------

p3a <- bio |>
  group_by(condition, fsf_signal_class) |>
  summarise(n_genes=sum(n_genes), .groups="drop") |>
  group_by(condition) |>
  mutate(prop=n_genes/sum(n_genes)) |>
  ungroup() |>
  ggplot(aes(condition, prop, fill=fsf_signal_class)) +
  geom_col(width=0.75) +
  labs(title="Figure 3A. Gene-weighted FSF signal-class composition",
       x="Condition", y="Proportion of genes", fill="FSF class") +
  theme_pub()

savep(p3a, "Figure3A_gene_weighted_signal_class_composition", 9, 5.5)

p3b <- bio |>
  group_by(condition, stability_level) |>
  summarise(n_genes=sum(n_genes), .groups="drop") |>
  group_by(condition) |>
  mutate(prop=n_genes/sum(n_genes)) |>
  ungroup() |>
  ggplot(aes(condition, prop, fill=stability_level)) +
  geom_col(width=0.75) +
  labs(title="Figure 3B. Gene-weighted stability-level composition",
       x="Condition", y="Proportion of genes", fill="Stability level") +
  theme_pub()

savep(p3b, "Figure3B_gene_weighted_stability_level_composition", 8, 5.2)

p3c <- arch |>
  ggplot(aes(condition, n_signal_classes, fill=signal_architecture_type)) +
  geom_col(width=0.7) +
  geom_text(aes(label=n_signal_classes), vjust=-0.3, size=4) +
  labs(title="Figure 3C. Condition-level FSF architecture",
       x="Condition", y="Number of FSF classes", fill="Architecture") +
  theme_pub()

savep(p3c, "Figure3C_condition_architecture", 8.5, 5)

# ---------------- Figure 4 ----------------

mfile <- file.path(ROOT, "results/fsf_metrics/fsf_ssi_stability_deviation.tsv")
if(file.exists(mfile)){
  m <- read_tsv(mfile, show_col_types=FALSE)
  names(m) <- make.names(names(m))
  ssi <- names(m)[grepl("ssi", names(m), ignore.case=TRUE)][1]
  dev <- names(m)[grepl("deviation", names(m), ignore.case=TRUE)][1]
  cls <- names(m)[grepl("class", names(m), ignore.case=TRUE)][1]

  if(!is.na(ssi)){
    p4a <- ggplot(m, aes(.data[[ssi]])) +
      geom_histogram(bins=50, fill="grey35", color="white") +
      labs(title="Figure 4A. SSI distribution", x="SSI", y="Number of features") +
      theme_pub()
    savep(p4a, "Figure4A_SSI_distribution", 7, 5)
  }

  if(!is.na(dev)){
    p4b <- ggplot(m, aes(.data[[dev]])) +
      geom_histogram(bins=50, fill="grey35", color="white") +
      labs(title="Figure 4B. Stability-deviation distribution",
           x="Stability deviation", y="Number of features") +
      theme_pub()
    savep(p4b, "Figure4B_stability_deviation_distribution", 7, 5)
  }

  if(!is.na(ssi) && !is.na(dev) && !is.na(cls)){
    p4c <- ggplot(m, aes(.data[[ssi]], .data[[dev]], color=.data[[cls]])) +
      geom_point(alpha=0.5, size=1) +
      labs(title="Figure 4C. SSI versus stability deviation",
           x="SSI", y="Stability deviation", color="Class") +
      theme_pub()
    savep(p4c, "Figure4C_SSI_vs_deviation", 8, 5.5)
  }
}

# ---------------- Figure 5 ----------------

p5a <- bio |>
  group_by(final_biological_theme) |>
  summarise(n_genes=sum(n_genes), .groups="drop") |>
  mutate(final_biological_theme=fct_reorder(final_biological_theme,n_genes)) |>
  ggplot(aes(final_biological_theme,n_genes)) +
  geom_col(fill="grey35") +
  coord_flip() +
  labs(title="Figure 5A. Gene-weighted curated biological themes",
       x="Theme", y="Number of genes") +
  theme_pub() +
  theme(axis.text.x=element_text(angle=0))

savep(p5a, "Figure5A_gene_weighted_theme_distribution", 8.5, 5.5)

p5b <- bio |>
  group_by(condition, final_biological_theme) |>
  summarise(n_genes=sum(n_genes), .groups="drop") |>
  group_by(condition) |>
  mutate(prop=n_genes/sum(n_genes)) |>
  ungroup() |>
  ggplot(aes(final_biological_theme, condition, fill=prop)) +
  geom_tile(color="white") +
  labs(title="Figure 5B. Gene-weighted biological theme landscape",
       x="Theme", y="Condition", fill="Gene proportion") +
  theme_pub()

savep(p5b, "Figure5B_gene_weighted_theme_heatmap", 11, 5.5)

p5c <- bio |>
  filter(condition %in% c("des","uv","gam")) |>
  ggplot(aes(fsf_signal_class, n_genes, fill=final_biological_theme)) +
  geom_col(width=0.75) +
  facet_wrap(~condition, scales="free_x") +
  labs(title="Figure 5C. Gene-weighted DES, UV and GAM biological profiles",
       x="FSF class", y="Number of genes", fill="Theme") +
  theme_pub()

savep(p5c, "Figure5C_gene_weighted_DES_UV_GAM_profiles", 12, 6)

# ---------------- Figure 6A ----------------

cfile <- file.path(ROOT, "results/synthetic/synthetic_recovery_confusion_matrix.tsv")
if(file.exists(cfile)){
  cm <- read_tsv(cfile, show_col_types=FALSE)
  names(cm) <- make.names(names(cm))
  p6a <- ggplot(cm, aes(.data[[names(cm)[2]]], .data[[names(cm)[1]]], fill=.data[[names(cm)[3]]])) +
    geom_tile(color="white") +
    geom_text(aes(label=.data[[names(cm)[3]]]), size=3.5) +
    labs(title="Figure 6A. Synthetic truth versus FSF class",
         x="Predicted FSF class", y="Synthetic truth", fill="Count") +
    theme_pub()
  savep(p6a, "Figure6A_synthetic_truth_confusion", 8, 6)
}

# ---------------- Figure 6B / 6C flexible recovery ----------------

ngdir <- file.path(ROOT, "results/synthetic/noise_gradient")
ngfiles <- list.files(ngdir, pattern="\\.tsv$", full.names=TRUE)
cat("Noise files found:\n"); print(basename(ngfiles))

if(length(ngfiles)>0){
  for(f in ngfiles){
    x <- read_tsv(f, show_col_types=FALSE)
    cat("\n", basename(f), "\n")
    print(names(x))
  }
}

f <- file.path(ngdir, "publication_noise_gradient_summary.tsv")
if(file.exists(f)){
  x <- read_tsv(f, show_col_types=FALSE)
  names(x) <- make.names(names(x))
  ncolx <- names(x)[grepl("noise", names(x), ignore.case=TRUE)][1]
  nums <- names(x)[sapply(x,is.numeric)]
  nums <- setdiff(nums, ncolx)
  if(!is.na(ncolx) && length(nums)>0){
    long <- x |> pivot_longer(cols=all_of(nums), names_to="metric", values_to="value")
    p6b <- ggplot(long, aes(.data[[ncolx]], value, linetype=metric)) +
      geom_line(linewidth=0.8) + geom_point(size=2) +
      labs(title="Figure 6B. Noise-gradient robustness",
           x="Noise level", y="Metric value", linetype="Metric") +
      theme_pub()
    savep(p6b, "Figure6B_noise_gradient_robustness", 8.5, 5)
  }
}

f <- file.path(ngdir, "publication_noise_gradient_dominant_recovery.tsv")
if(file.exists(f)){
  x <- read_tsv(f, show_col_types=FALSE)
  names(x) <- make.names(names(x))
  ncolx <- names(x)[grepl("noise", names(x), ignore.case=TRUE)][1]
  cls <- names(x)[grepl("class|truth|state", names(x), ignore.case=TRUE)][1]
  val <- names(x)[grepl("recover|accuracy|rate|proportion|mean|value", names(x), ignore.case=TRUE)][1]
  if(!is.na(ncolx) && !is.na(cls) && !is.na(val)){
    p6c <- ggplot(x, aes(.data[[ncolx]], .data[[val]], linetype=.data[[cls]])) +
      geom_line(linewidth=0.8) + geom_point(size=2) +
      labs(title="Figure 6C. Class-specific recovery across noise",
           x="Noise level", y="Recovery", linetype="Class") +
      theme_pub()
    savep(p6c, "Figure6C_class_specific_noise_recovery", 8.5, 5)
  } else {
    cat("Figure6C not generated. Detected columns:\n")
    print(c(noise=ncolx, class=cls, value=val))
  }
}

# ---------------- Figure 6D ----------------

tdir <- file.path(ROOT, "results/synthetic/tau_sensitivity")
tf <- list.files(tdir, pattern="\\.tsv$", full.names=TRUE)
if(length(tf)>0){
  tau <- bind_rows(lapply(tf, function(z) read_tsv(z, show_col_types=FALSE)))
  names(tau) <- make.names(names(tau))
  tcol <- names(tau)[grepl("tau|threshold", names(tau), ignore.case=TRUE)][1]
  cls <- names(tau)[grepl("class", names(tau), ignore.case=TRUE)][1]
  val <- names(tau)[sapply(tau,is.numeric)][1]
  if(!is.na(tcol) && !is.na(val)){
    p6d <- ggplot(tau, aes(.data[[tcol]], .data[[val]], linetype=.data[[cls]])) +
      geom_line(linewidth=0.8) + geom_point(size=2) +
      labs(title="Figure 6D. Tau sensitivity",
           x="Tau / threshold", y="Value", linetype="Class") +
      theme_pub()
    savep(p6d, "Figure6D_tau_sensitivity", 8, 5)
  }
}

# ---------------- Figure 7 ----------------

model <- tibble(
  x=c(1,2,3), y=1,
  label=c(
    "Stable-dominated\nLT, HT, OSM\nMostly high-stability signal classes",
    "Mixed stable-instability\nDES\nStable and instable programs coexist",
    "Multi-layer architecture\nUV, GAM\nStable, transitional and instable layers coexist"
  )
)

p7 <- ggplot(model, aes(x,y)) +
  geom_segment(
    data=model |> filter(x<3),
    aes(x=x+0.25, xend=x+0.75, y=y, yend=y),
    arrow=arrow(length=unit(0.18,"cm")),
    linewidth=0.5
  ) +
  geom_label(aes(label=label), size=3.35, label.size=0.25, fill="white") +
  scale_x_continuous(limits=c(0.55,3.45), breaks=NULL) +
  scale_y_continuous(limits=c(0.92,1.08), breaks=NULL) +
  labs(title="Figure 7. Final FSF interpretation model") +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=0))

savep(p7, "Figure7_Final_Model", 10.5, 3)

cat("\nFinal inventory:\n")
print(sort(list.files(OUT, pattern="\\.pdf$")))

cat("\nFinished:", as.character(Sys.time()), "\n")
sink()
