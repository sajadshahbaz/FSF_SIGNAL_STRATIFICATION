#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "source_data")
figure_dir <- file.path(root, "results", "current_fsf_v1", "manuscript", "figures")
builder <- file.path(root, "scripts", "current", "figures", "build_current_manuscript_figures.R")

figure_ids <- c(paste0("Figure", 1:7), paste0("FigureS", 1:4))
source_files <- c(
  "figure2_simplex_source.tsv", "figure3_architecture_source.tsv",
  "figure4_ssi_ecdf_source.tsv", "figure5_go_source.tsv",
  "figure5_kegg_source.tsv", "figure6_baseline_source.tsv",
  "figure6_noise_source.tsv", "figureS2_tau_source.tsv",
  "figureS3_architecture_source.tsv", "figureS4_ssi_distribution_source.tsv"
)
stopifnot(all(file.exists(file.path(source_dir, source_files))))
stopifnot(file.exists(builder))

outputs <- unlist(lapply(figure_ids, function(id) file.path(figure_dir, paste0(id, c(".png", ".pdf")))))
stopifnot(all(file.exists(outputs)))
stopifnot(all(file.info(outputs)$size > 1000L))

manifest_path <- file.path(figure_dir, "figure_manifest.tsv")
audit_path <- file.path(figure_dir, "figure_audit.md")
stopifnot(file.exists(manifest_path), file.exists(audit_path))
manifest <- read.delim(manifest_path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
expected_manifest_columns <- c(
  "figure_id", "panel", "source_authority", "generation_script",
  "output_png", "output_pdf", "numerical_status",
  "interpretation_status", "notes"
)
stopifnot(identical(names(manifest), expected_manifest_columns))
stopifnot(identical(manifest$figure_id, figure_ids))
stopifnot(anyDuplicated(manifest$figure_id) == 0L)

architecture <- read.delim(
  file.path(source_dir, "figure3_architecture_source.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
expected <- data.frame(
  condition = rep(c("DES", "GAM", "HT", "LT", "OSM", "UV"), each = 4),
  stability_region = rep(c("Low Stability", "Transitional", "Stable", "Highly Stable"), 6),
  region_count = c(
    5837, 0, 0, 9350,
    1085, 2106, 6127, 5869,
    0, 0, 0, 15187,
    1211, 0, 0, 13976,
    0, 0, 0, 15187,
    2948, 6697, 3151, 2391
  ),
  stringsAsFactors = FALSE
)
observed <- architecture[, c("condition", "stability_region", "region_count")]
observed <- observed[order(match(observed$condition, unique(expected$condition)),
                           match(observed$stability_region, unique(expected$stability_region))), ]
rownames(observed) <- NULL
stopifnot(identical(observed$condition, expected$condition))
stopifnot(identical(observed$stability_region, expected$stability_region))
stopifnot(identical(as.integer(observed$region_count), as.integer(expected$region_count)))
stopifnot(all(rowsum(architecture$region_count, architecture$condition)[, 1] == 15187L))
uv <- architecture[architecture$condition == "UV", ]
stopifnot(uv$stability_region[which.max(uv$region_count)] == "Transitional")

tau <- read.delim(
  file.path(source_dir, "figureS2_tau_source.tsv"),
  check.names = FALSE, quote = "", stringsAsFactors = FALSE
)
stopifnot(identical(names(tau), c("tau", "truth_scenario", "metric", "category", "value")))
builder_text <- readLines(builder, warn = FALSE)
stopifnot(any(grepl("ggplot(s2, aes(tau, value", builder_text, fixed = TRUE)))
stopifnot(!any(grepl("ggplot(s2, aes(tau, tau", builder_text, fixed = TRUE)))

pdf_text <- character()
for (pdf in file.path(figure_dir, paste0(figure_ids, ".pdf"))) {
  txt <- tempfile(fileext = ".txt")
  status <- system2("pdftotext", c(pdf, txt))
  stopifnot(status == 0L, file.exists(txt))
  pdf_text <- c(pdf_text, readLines(txt, warn = FALSE))
  unlink(txt)
}
forbidden_terms <- c("Instability", "instable", "weakly stable", "weakly_stable")
pdf_text_lower <- tolower(pdf_text)
stopifnot(!any(vapply(tolower(forbidden_terms), function(x) any(grepl(x, pdf_text_lower, fixed = TRUE)), logical(1))))

figure5_text_file <- tempfile(fileext = ".txt")
stopifnot(system2("pdftotext", c(file.path(figure_dir, "Figure5.pdf"), figure5_text_file)) == 0L)
figure5_text <- readLines(figure5_text_file, warn = FALSE)
unlink(figure5_text_file)
dropped <- c(
  "PROFILE_002", "PROFILE_003", "PROFILE_008",
  "PROFILE_009", "PROFILE_013", "PROFILE_026",
  "GO_THEME_001", "KEGG_THEME_001", "KEGG_THEME_002",
  "Other KEGG", "manual review"
)
stopifnot(!any(vapply(dropped, function(x) any(grepl(x, figure5_text, fixed = TRUE)), logical(1))))

audit <- readLines(audit_path, warn = FALSE)
stopifnot(any(grepl("tau-against-tau", audit, fixed = TRUE)))
stopifnot(any(grepl("x = `tau` and y = `value`", audit, fixed = TRUE)))
stopifnot(any(grepl("UV shown as Transitional-dominant", audit, fixed = TRUE)))

protected <- c(
  "R", "results/current_fsf_v1/manuscript/source_data",
  "results/current_fsf_v1/manuscript/review",
  "scripts/current/audit_current_manuscript_bundle.R",
  "scripts/current/build_current_manuscript_tables.R"
)
status <- system2("git", c("diff", "--quiet", "--", protected))
stopifnot(status == 0L)

cat("current manuscript figure tests: PASS\n")
