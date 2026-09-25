#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

ROOT <- "/media/saji/5E06441D0643F5152/FSF_SIGNAL_STRATIFICATION"

MAIN <- file.path(ROOT,"results/archive/pre_repair_manuscript/figures/main")
SUPP <- file.path(ROOT,"results/archive/pre_repair_manuscript/figures/supplementary")

theme_pub <- function() {
  theme_bw(11) +
    theme(
      plot.title = element_text(face="bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle=45,hjust=1)
    )
}

# ====================================================
# FIGURE 6B
# Noise Gradient
# ====================================================

f <- file.path(
  ROOT,
  "results/synthetic/noise_gradient/publication_noise_gradient_summary.tsv"
)

if(file.exists(f)) {

  x <- read_tsv(f, show_col_types = FALSE)

  names(x) <- make.names(names(x))

  noise_col <- names(x)[str_detect(names(x),"noise")][1]

  num_cols <- names(x)[sapply(x,is.numeric)]
  num_cols <- setdiff(num_cols, noise_col)

  if(length(num_cols)>0){

    long <- x %>%
      pivot_longer(
        cols = all_of(num_cols),
        names_to="metric",
        values_to="value"
      )

    p <- ggplot(
      long,
      aes(
        x=.data[[noise_col]],
        y=value,
        linetype=metric
      )
    )+
      geom_line(linewidth=0.9)+
      geom_point(size=2)+
      labs(
        title="Figure 6B. Noise-gradient robustness",
        x="Noise level",
        y="Metric value",
        linetype="Metric"
      )+
      theme_pub()

    ggsave(
      file.path(MAIN,"Figure6B_noise_gradient.pdf"),
      p,width=8,height=5
    )

    ggsave(
      file.path(MAIN,"Figure6B_noise_gradient.png"),
      p,width=8,height=5,dpi=300
    )
  }
}

# ====================================================
# FIGURE 6C
# Recovery by Signal Class
# ====================================================

f <- file.path(
 ROOT,
 "results/synthetic/noise_gradient/publication_noise_gradient_dominant_recovery.tsv"
)

if(file.exists(f)) {

  x <- read_tsv(f, show_col_types = FALSE)

  names(x) <- make.names(names(x))

  noise_col <- names(x)[str_detect(names(x),"noise")][1]

  value_col <- names(x)[
    str_detect(
      names(x),
      "recover|accuracy|rate|mean|proportion"
    )
  ][1]

  class_col <- names(x)[
    str_detect(
      names(x),
      "class|truth|state"
    )
  ][1]

  if(
    !is.na(noise_col) &
    !is.na(value_col)
  ){

    p <- ggplot(
      x,
      aes(
        x=.data[[noise_col]],
        y=.data[[value_col]],
        linetype=.data[[class_col]]
      )
    )+
      geom_line(linewidth=1)+
      geom_point(size=2)+
      labs(
        title="Figure 6C. Recovery across signal classes",
        x="Noise level",
        y="Recovery",
        linetype="Signal class"
      )+
      theme_pub()

    ggsave(
      file.path(MAIN,"Figure6C_class_recovery.pdf"),
      p,width=8,height=5
    )

    ggsave(
      file.path(MAIN,"Figure6C_class_recovery.png"),
      p,width=8,height=5,dpi=300
    )
  }
}

# ====================================================
# SUPPLEMENTARY S2
# Tau sensitivity
# ====================================================

tau_dir <- file.path(
 ROOT,
 "results/synthetic/tau_sensitivity"
)

if(dir.exists(tau_dir)){

  ff <- list.files(
    tau_dir,
    pattern="\\.tsv$",
    full.names=TRUE
  )

  if(length(ff)>0){

    dat <- bind_rows(
      lapply(ff,function(z){

        read_tsv(
          z,
          show_col_types = FALSE
        ) %>%
        mutate(source_file=basename(z))

      })
    )

    names(dat) <- make.names(names(dat))

    tau_col <- names(dat)[
      str_detect(
        names(dat),
        "tau|threshold"
      )
    ][1]

    val_col <- names(dat)[
      sapply(dat,is.numeric)
    ][1]

    p <- ggplot(
      dat,
      aes(
        x=.data[[tau_col]],
        y=.data[[val_col]]
      )
    )+
      geom_line()+
      geom_point()+
      facet_wrap(~source_file,scales="free_y")+
      labs(
        title="Supplementary Figure S2. Tau sensitivity",
        x="Tau",
        y="Metric"
      )+
      theme_pub()

    ggsave(
      file.path(SUPP,"Supplementary_FigureS2_tau_sensitivity.pdf"),
      p,width=10,height=6
    )

    ggsave(
      file.path(SUPP,"Supplementary_FigureS2_tau_sensitivity.png"),
      p,width=10,height=6,dpi=300
    )
  }
}

cat("Benchmark figures generated.\n")
