#!/usr/bin/env Rscript

# Build evidence-only biological interpretation review tables. Historical
# mappings are inventory records, never automatically accepted authorities.

.review_root <- file.path("results", "current_fsf_v1", "manuscript")
.review_out <- file.path(.review_root, "review")

.write_review <- function(x, name, root = .review_out) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, file.path(root, name), na = "NA")
}

.historical_rules <- function() {
  add <- function(prefix, source, type, labels, rules, targets, mode = "manual") {
    data.frame(
      historical_rule_id = sprintf("%s_%02d", prefix, seq_along(labels)),
      source_script = source, rule_type = type, historical_label = labels,
      mapping_or_rule = rules, historical_target = targets,
      manual_or_deterministic = mode,
      current_review_status = "REQUIRES_AUTHOR_REVIEW",
      stringsAsFactors = FALSE
    )
  }
  go_labels <- c(
    "Translation / ribosome", "Protein targeting / secretion",
    "Genome maintenance / DNA repair", "Chromosome / nuclear organization",
    "Cell cycle / cell division", "Development / morphogenesis",
    "Neural / projection organization", "Movement / behavior / taxis",
    "Cell adhesion / junction", "Metabolism / oxidation-reduction",
    "Broad regulatory / cellular process", "Unclassified / manual review"
  )
  go_rules <- c(
    "GO:0022626;GO:0022625;GO:0022627;GO:0005840;GO:0003735;GO:0002181;GO:0006412;GO:0043043",
    "GO:0006612;GO:0006613;GO:0006614;GO:0045047;GO:0072599;GO:0070972",
    "GO:0006259;GO:0006260;GO:0006271;GO:0006273;GO:0006281;GO:0006974;GO:0000723;GO:0000724;GO:0000725;GO:0000726",
    "GO:0005694;GO:0051276;GO:0000228;GO:0044427;GO:0044454;GO:0098813;GO:0005634;GO:0005654",
    "GO:0007049;GO:0022402;GO:0000278;GO:0000280;GO:0048285;GO:0140013;GO:0051321;GO:1903046;GO:1903047;GO:0000070;GO:0000075;GO:0000077",
    "GO:0035239;GO:0048729;GO:0060562;GO:0009887;GO:0048598;GO:0060541;GO:0001501;GO:0002009",
    "GO:0007411;GO:0097485;GO:0031175;GO:0048812;GO:0061564",
    "GO:0040011;GO:0006935;GO:0050920;GO:0007610;GO:0008045;GO:0008038",
    "GO:0007155;GO:0007156;GO:0098609;GO:0005912;GO:0098742",
    "GO:0016491;GO:0055114;GO:0044281;GO:0006082;GO:0019752;GO:0043436;GO:0071704;GO:0044238",
    "regex:^GO:00", "fallback:TRUE"
  )
  kegg_labels <- c(
    "Translation / ribosome", "Genome maintenance / DNA repair",
    "Cell cycle regulation", "Protein turnover / proteostasis",
    "Energy / core metabolism", "Stress signaling / signal transduction",
    "Cell adhesion / extracellular interaction", "Cytoskeleton / cellular remodeling",
    "Neural signaling", "Endocrine / physiological regulation",
    "Broad proliferation / regulatory pathway", "Other KEGG / manual review"
  )
  kegg_rules <- c(
    "ko03010|map03010|Ribosome|K029|K028",
    "ko03030|map03030|ko034|map034|K107|K108|K109|DNA|repair|replication|recombination",
    "ko04110|map04110|ko04111|map04111|ko04113|map04113|cell cycle",
    "ko03050|map03050|ko04120|map04120|proteasome|ubiquitin|K056",
    "ko00190|map00190|ko01100|map01100|ko01110|map01110|ko01120|map01120|ko01200|map01200|oxidative|metabolism",
    "ko04010|map04010|ko04013|map04013|ko04022|map04022|ko04151|map04151|ko046|map046|signaling|MAPK",
    "ko045|map045|adhesion|junction|ECM", "ko048|map048|cytoskeleton|actin",
    "ko047|map047|synapse|neuro|neuron", "ko049|map049|endocrine|hormone|insulin",
    "ko052|map052|cancer", "fallback:TRUE"
  )
  collapse_labels <- c(
    "Translation / protein targeting", "Genome maintenance / nuclear regulation",
    "Cell cycle / division", "Proteostasis / protein turnover",
    "Energy metabolism / redox regulation", "Stress signaling / signal transduction",
    "Cell adhesion / extracellular interaction", "Development / morphogenesis",
    "Neural / behavioral organization", "No clear enrichment",
    "Weak / nonspecific enrichment", "Broad cellular regulation", "Mixed functional program"
  )
  collapse_rules <- c(
    "translation|ribosome|ribosomal|protein targeting|secretion",
    "dna repair|genome maintenance|replication|recombination|chromosome|nuclear organization|nucleotide excision|homologous recombination|mismatch repair",
    "cell cycle|cell division|mitotic|meiosis|spindle|chromosome segregation",
    "proteostasis|protein turnover|proteasome|ubiquitin|folding|chaperone",
    "energy|metabolism|oxidation|oxidative|oxidoreduction|mitochond|carbon metabolism|oxidative phosphorylation",
    "stress signaling|signal transduction|mapk|calcium|phosphatidylinositol|signaling",
    "adhesion|junction|extracellular|ecm|cell adhesion",
    "development|morphogenesis|tissue|structural remodeling|organ development",
    "neural|neuron|synapse|axon|projection|behavior|movement|taxis|locomotion",
    "both GO and KEGG report no significant enrichment",
    "no significant plus broad/other/manual/unclassified", "broad/other/manual/unclassified",
    "fallback:TRUE"
  )
  rbind(
    add("GO", "scripts/pipeline/32_GO_enrichment_interpretation_tables.R",
        "GO_THEME_MAPPING", go_labels, go_rules, go_labels),
    add("KEGG", "scripts/pipeline/34_build_biological_profiles.R",
        "KEGG_THEME_MAPPING", kegg_labels, kegg_rules, kegg_labels),
    add("PROFILE", "scripts/pipeline/35_curate_final_biological_themes.R",
        "COMBINED_PROFILE_MAPPING", collapse_labels, collapse_rules, collapse_labels),
    add("ARCH", "scripts/pipeline/35_curate_final_biological_themes.R",
        "CONDITION_ARCHITECTURE", c("Highly stable-dominated architecture",
          "Mixed stable-instabile architecture", "Multi-layer signal architecture",
          "Unassigned architecture"),
        c("condition in LT,OSM,HT", "condition == DES", "condition in UV,GAM", "fallback"),
        c("LT;OSM;HT", "DES", "UV;GAM", "other")),
    add("PRIORITY", "scripts/pipeline/35_curate_final_biological_themes.R",
        "EVIDENCE_PRIORITY", c("low_no_clear", "high", "moderate", "limited", "low_fallback"),
        c("theme == No clear enrichment", "GO >100 and KEGG >10", "GO >20 or KEGG >5",
          "GO >0 or KEGG >0", "fallback"), c("low", "high", "moderate", "limited", "low")),
    add("SUPPORT", "scripts/pipeline/34_build_biological_profiles.R",
        "EVIDENCE_PRESENCE", c("both", "GO_only", "KEGG_only", "neither"),
        c("GO>0 and KEGG>0", "GO>0 and KEGG=0", "GO=0 and KEGG>0", "GO=0 and KEGG=0"),
        c("GO_and_KEGG_supported", "GO_supported", "KEGG_supported", "no_significant_enrichment"),
        "deterministic"),
    add("NARRATIVE34", "scripts/pipeline/34_build_biological_profiles.R",
        "GENERATED_NARRATIVE", "signal-class profile sentence",
        "condition/class shows profile; GO theme; KEGG theme", "manuscript_interpretation"),
    add("NARRATIVE35", "scripts/pipeline/35_curate_final_biological_themes.R",
        "GENERATED_NARRATIVE", "curated support sentence",
        "condition/class is assigned theme based on evidence; support level", "final_interpretation")
  )
}

.rule_evidence <- function(rules, terms, ids, type) {
  out <- rules[rules$rule_type == type, , drop = FALSE]
  out$current_matching_term_count <- 0L
  out$current_matching_gene_set_count <- 0L
  out$current_matching_terms <- NA_character_
  out$evidence_classification <- "AMBIGUOUS_REQUIRES_AUTHOR_REVIEW"
  for (i in seq_len(nrow(out))) {
    rule <- out$mapping_or_rule[[i]]
    hit <- if (startsWith(rule, "fallback:")) rep(FALSE, length(terms)) else if (type == "GO_THEME_MAPPING" && !startsWith(rule, "regex:")) {
      terms %in% strsplit(rule, ";", fixed = TRUE)[[1L]]
    } else {
      pattern <- sub("^regex:", "", rule)
      grepl(pattern, terms, ignore.case = TRUE)
    }
    matched <- sort(unique(terms[hit]), method = "radix")
    out$current_matching_term_count[[i]] <- length(matched)
    out$current_matching_gene_set_count[[i]] <- length(unique(ids[hit]))
    out$current_matching_terms[[i]] <- if (length(matched)) paste(matched, collapse = ";") else NA_character_
    if (!startsWith(rule, "fallback:") && length(matched)) {
      out$evidence_classification[[i]] <- "SUPPORTED_UNCHANGED"
    } else if (!startsWith(rule, "fallback:")) {
      out$evidence_classification[[i]] <- "NO_CURRENT_EVIDENCE"
    }
  }
  out$author_decision <- NA_character_
  out$author_notes <- NA_character_
  out
}

.compact_evidence <- function(x, term, n = 10L) {
  if (!nrow(x)) return(x)
  groups <- split(seq_len(nrow(x)), x$gene_set_id)
  out <- do.call(rbind, lapply(groups, function(i) {
    y <- x[i, , drop = FALSE]
    y <- y[order(y$adjusted_p_value, y$p_value, -y$enrichment_ratio,
                 y[[term]], method = "radix"), , drop = FALSE]
    y <- y[seq_len(min(n, nrow(y))), , drop = FALSE]
    y$review_rank <- seq_len(nrow(y))
    y$contributing_feature_count <- lengths(strsplit(y$contributing_feature_ids, ";", fixed = TRUE))
    y$term_description <- NA_character_
    y$review_status <- "REQUIRES_AUTHOR_REVIEW"
    y
  }))
  rownames(out) <- NULL
  out
}

.complete_evidence <- function(x, main) {
  missing <- main[!main$gene_set_id %in% x$gene_set_id, , drop = FALSE]
  if (!nrow(missing)) return(x)
  blank <- x[rep(NA_integer_, nrow(missing)), , drop = FALSE]
  common <- intersect(names(blank), names(missing))
  blank[common] <- missing[common]
  blank$review_status <- "REQUIRES_AUTHOR_REVIEW"
  blank$term_description <- NA_character_
  out <- rbind(x, blank)
  out[order(out$gene_set_id, out$review_rank, method = "radix", na.last = TRUE), , drop = FALSE]
}

.evidence_summary <- function(x, term, prefix) {
  sets <- unique(x$gene_set_id)
  rows <- lapply(sets, function(id) {
    y <- x[x$gene_set_id == id, , drop = FALSE]
    data.frame(gene_set_id = id,
      summary = paste0(prefix, " significant terms=", nrow(y), "; top IDs=",
                       paste(head(y[[term]], 5L), collapse = ";"),
                       "; minimum adjusted p=", format(min(y$adjusted_p_value), scientific = TRUE)),
      stringsAsFactors = FALSE)
  })
  if (!length(rows)) return(data.frame(gene_set_id = character(), summary = character()))
  do.call(rbind, rows)
}

.architecture_review <- function(region, signal) {
  conditions <- c("DES", "GAM", "HT", "LT", "OSM", "UV")
  historical <- c(DES = "Mixed stable-instabile architecture", GAM = "Multi-layer signal architecture",
                  HT = "Highly stable-dominated architecture", LT = "Highly stable-dominated architecture",
                  OSM = "Highly stable-dominated architecture", UV = "Multi-layer signal architecture")
  do.call(rbind, lapply(conditions, function(cond) {
    x <- region[region$condition == cond, ]
    dominant <- x$stability_region[x$region_count == max(x$region_count)]
    minor <- x$stability_region[x$region_count > 0 & x$region_count < max(x$region_count)]
    s <- signal[signal$condition == cond, ]
    dir <- aggregate(class_count ~ dominant_state, s, sum)
    data.frame(condition = cond,
      historical_architecture_label = historical[[cond]],
      dominant_regions = paste(dominant, collapse = ";"),
      minor_regions = if (length(minor)) paste(minor, collapse = ";") else "none",
      transitional_count = x$region_count[x$stability_region == "Transitional"],
      directional_composition = paste0(dir$dominant_state, "=", dir$class_count, collapse = ";"),
      directly_supported = paste0("observed region counts: ", paste0(x$stability_region, "=", x$region_count, collapse = ";")),
      not_supported = "causal mechanism, biological activation, adaptation, or fitness effect",
      historical_profile_status = if (cond == "DES") "SUPPORTED_BUT_LABEL_UPDATE" else if (cond == "UV") "PARTIALLY_SUPPORTED" else "AMBIGUOUS_REQUIRES_AUTHOR_REVIEW",
      review_status = "REQUIRES_AUTHOR_REVIEW", author_decision = NA_character_, author_notes = NA_character_,
      stringsAsFactors = FALSE)
  }))
}

.representative_review <- function(go, kegg) {
  path <- file.path("results", "current_fsf_v1", "representative_features",
                    "current_fsf_top50_annotated_stable_directional_signals_compact.tsv")
  reps <- readr::read_tsv(path, show_col_types = FALSE)
  membership_count <- function(ids, evidence) {
    vapply(ids, function(id) sum(vapply(strsplit(evidence$contributing_feature_ids, ";", fixed = TRUE),
                                        function(z) id %in% z, logical(1L))), integer(1L))
  }
  reps$current_go_supporting_term_count <- membership_count(reps$feature_id, go)
  reps$current_kegg_supporting_term_count <- membership_count(reps$feature_id, kegg)
  reps$historical_selection_reason <- "frozen annotated Stable/Highly Stable directional representative ranking"
  reps$current_evidence_status <- ifelse(
    reps$current_go_supporting_term_count + reps$current_kegg_supporting_term_count > 0,
    "RETAIN", "RETAIN_WITH_REVISED_RATIONALE")
  reps$review_status <- "REQUIRES_AUTHOR_REVIEW"
  reps$author_decision <- NA_character_
  reps$author_notes <- NA_character_
  reps
}

build_current_biological_review_package <- function(root = .review_out) {
  if (!requireNamespace("readr", quietly = TRUE)) stop("readr is required.")
  manifest <- readr::read_tsv(file.path(.review_root, "tables/current_enrichment_gene_set_summary.tsv"), show_col_types = FALSE)
  main <- manifest[manifest$gene_set_family == "main_class", ]
  go_all <- readr::read_tsv(file.path(.review_root, "enrichment/go/go_enrichment_significant.tsv"), show_col_types = FALSE)
  kegg_all <- readr::read_tsv(file.path(.review_root, "enrichment/kegg/kegg_enrichment_significant.tsv"), show_col_types = FALSE)
  go <- go_all[go_all$gene_set_id %in% main$gene_set_id, ]
  kegg <- kegg_all[kegg_all$gene_set_id %in% main$gene_set_id, ]
  go_review <- .complete_evidence(.compact_evidence(go, "go_term", 10L), main)
  kegg_review <- .complete_evidence(.compact_evidence(kegg, "kegg_term", 10L), main)
  rules <- .historical_rules()
  go_rules <- .rule_evidence(rules, go$go_term, go$gene_set_id, "GO_THEME_MAPPING")
  kegg_rules <- .rule_evidence(rules, kegg$kegg_term, kegg$gene_set_id, "KEGG_THEME_MAPPING")

  go_sum <- .evidence_summary(go, "go_term", "GO")
  kg_sum <- .evidence_summary(kegg, "kegg_term", "KEGG")
  names(go_sum)[2L] <- "current_GO_evidence_summary"
  names(kg_sum)[2L] <- "current_KEGG_evidence_summary"
  historical <- readr::read_tsv(file.path("results", "manuscript", "tables", "main",
                                           "Table_SignalClass_Biology_CURATED.tsv"), show_col_types = FALSE)
  historical$key <- paste(toupper(historical$condition),
    ifelse(historical$stability_level == "instable", "Low Stability",
      ifelse(historical$stability_level == "weakly_stable_transitional", "Transitional",
        ifelse(historical$stability_level == "highly_stable", "Highly Stable", "Stable"))),
    ifelse(historical$dominant_state == "mixed", "tied", historical$dominant_state), sep = "\r")
  main$key <- paste(main$condition, main$stability_region, main$dominant_state, sep = "\r")
  hi <- match(main$key, historical$key)
  profile <- main[c("condition", "stability_region", "dominant_state", "signal_class", "gene_set_id")]
  profile$region_or_class <- ifelse(is.na(profile$signal_class), profile$stability_region, profile$signal_class)
  profile$historical_profile <- historical$final_biological_theme[hi]
  profile$historical_GO_theme <- historical$GO_theme[hi]
  profile$historical_KEGG_theme <- historical$KEGG_theme[hi]
  profile$historical_GO_significant_count <- historical$significant_GO_terms_FDR005[hi]
  profile$historical_KEGG_significant_count <- historical$significant_KEGG_terms_FDR005[hi]
  profile <- merge(profile, go_sum, by = "gene_set_id", all.x = TRUE, sort = FALSE)
  profile <- merge(profile, kg_sum, by = "gene_set_id", all.x = TRUE, sort = FALSE)
  profile$current_GO_evidence_summary[is.na(profile$current_GO_evidence_summary)] <- "GO significant terms=0"
  profile$current_KEGG_evidence_summary[is.na(profile$current_KEGG_evidence_summary)] <- "KEGG significant terms=0"
  profile$architecture_context <- paste(profile$condition, profile$stability_region, profile$dominant_state, sep = "/")
  go_n <- table(go$gene_set_id)[profile$gene_set_id]; go_n[is.na(go_n)] <- 0
  kg_n <- table(kegg$gene_set_id)[profile$gene_set_id]; kg_n[is.na(kg_n)] <- 0
  profile$current_GO_significant_count <- as.integer(go_n)
  profile$current_KEGG_significant_count <- as.integer(kg_n)
  compare_count <- function(old, current) {
    ifelse(is.na(old), "NO_HISTORICAL_COMPARATOR",
      ifelse(old == 0 & current == 0, "ABSENT_UNCHANGED",
        ifelse(old == 0 & current > 0, "NEW_EVIDENCE",
          ifelse(old > 0 & current == 0, "LOST_EVIDENCE",
            ifelse(current > old, "STRENGTHENED_TERM_COUNT",
              ifelse(current < old, "WEAKENED_TERM_COUNT", "RETAINED_TERM_COUNT"))))))
  }
  profile$GO_evidence_change <- compare_count(profile$historical_GO_significant_count, profile$current_GO_significant_count)
  profile$KEGG_evidence_change <- compare_count(profile$historical_KEGG_significant_count, profile$current_KEGG_significant_count)
  profile$theme_evidence_change <- paste0("GO_", profile$GO_evidence_change, ";KEGG_", profile$KEGG_evidence_change)
  profile$historical_profile_status <- ifelse(
    profile$condition == "UV", "AMBIGUOUS_REQUIRES_AUTHOR_REVIEW",
    ifelse(go_n + kg_n == 0 & profile$historical_GO_significant_count + profile$historical_KEGG_significant_count > 0,
      "NO_LONGER_SUPPORTED",
      ifelse(go_n + kg_n == 0, "NO_CURRENT_EVIDENCE",
        ifelse(profile$stability_region == "Low Stability", "SUPPORTED_BUT_LABEL_UPDATE",
          ifelse(profile$GO_evidence_change == "RETAINED_TERM_COUNT" & profile$KEGG_evidence_change == "RETAINED_TERM_COUNT",
            "SUPPORTED_UNCHANGED", "PARTIALLY_SUPPORTED")))))
  profile$review_status <- "REQUIRES_AUTHOR_REVIEW"
  profile$author_decision <- NA_character_
  profile$author_notes <- NA_character_

  region <- readr::read_tsv(file.path(.review_root, "tables/current_region_architecture.tsv"), show_col_types = FALSE)
  signal <- readr::read_tsv(file.path(.review_root, "tables/current_signal_class_architecture.tsv"), show_col_types = FALSE)
  architecture <- .architecture_review(region, signal)
  reps <- .representative_review(go_all, kegg_all)

  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  .write_review(rules, "historical_interpretation_inventory.tsv", root)
  .write_review(go_review, "current_go_evidence.tsv", root)
  .write_review(kegg_review, "current_kegg_evidence.tsv", root)
  .write_review(go_rules, "go_theme_review.tsv", root)
  .write_review(kegg_rules, "kegg_theme_review.tsv", root)
  .write_review(profile, "combined_profile_review.tsv", root)
  .write_review(profile[profile$condition == "UV", ], "uv_specific_review.tsv", root)
  .write_review(architecture, "condition_architecture_review.tsv", root)
  .write_review(reps, "representative_feature_review.tsv", root)

  statuses <- table(profile$historical_profile_status)
  summary <- c(
    "# FSF biological interpretation review package", "",
    "This directory contains evidence views and historical-rule candidates for author review. It contains no accepted biological theme authority and no final manuscript prose.", "",
    paste0("- Current main-class gene sets reviewed: ", nrow(main)),
    paste0("- Historical interpretation rules inventoried: ", nrow(rules)),
    paste0("- GO historical mappings with current matching evidence: ", sum(go_rules$current_matching_term_count > 0)),
    paste0("- KEGG historical mappings with current matching evidence: ", sum(kegg_rules$current_matching_term_count > 0)),
    paste0("- Profile status counts: ", paste(names(statuses), as.integer(statuses), sep = "=", collapse = "; ")), "",
    "## Historical mapping evidence", "",
    paste0("- GO mappings with matching current terms: ", paste(go_rules$historical_label[go_rules$evidence_classification == "SUPPORTED_UNCHANGED"], collapse = "; ")),
    paste0("- KEGG mappings with matching current terms: ", paste(kegg_rules$historical_label[kegg_rules$evidence_classification == "SUPPORTED_UNCHANGED"], collapse = "; ")),
    paste0("- Historical mappings without current matches: ", paste(kegg_rules$historical_label[kegg_rules$evidence_classification == "NO_CURRENT_EVIDENCE"], collapse = "; ")),
    "- Broad/default and combined-profile mappings remain manual biological judgments even where their source terms are present.", "",
    "## Condition review", "",
    "- DES, GAM, and LT contain Low Stability rows requiring the locked label migration; their non-UV memberships are unchanged.",
    "- HT and OSM are entirely Highly Stable in the region authority; biological wording still requires author review.",
    "- Statistical count differences versus historical curated tables are recorded per gene set in combined_profile_review.tsv.", "",
    "## UV review gate", "",
    "UV is Transitional-dominated (6,697 of 15,187); its historical profiles remain REQUIRES_AUTHOR_REVIEW because current class membership differs materially.", "",
    "## Figure 5 candidates", "",
    "Candidate mappings are the GO/KEGG rules marked SUPPORTED_UNCHANGED in the review tables. Presence of matching enriched terms does not constitute author acceptance of a theme.", "",
    "## Figure 7 candidate structure", "",
    "Current condition architecture -> current signal class -> GO/KEGG statistical evidence -> author-reviewed theme/profile. Causal or adaptive language is outside the automatic evidence layer.", "",
    "Author decisions and notes are intentionally blank."
  )
  writeLines(summary, file.path(root, "review_summary.md"), useBytes = TRUE)
  invisible(list(rules = rules, go = go_review, kegg = kegg_review,
                 profiles = profile, architecture = architecture, representatives = reps))
}

audit_current_biological_review_package <- function(root = .review_out) {
  required <- c("historical_interpretation_inventory.tsv", "current_go_evidence.tsv",
    "current_kegg_evidence.tsv", "go_theme_review.tsv", "kegg_theme_review.tsv",
    "combined_profile_review.tsv", "uv_specific_review.tsv", "condition_architecture_review.tsv",
    "representative_feature_review.tsv", "review_summary.md")
  if (!all(file.exists(file.path(root, required)))) stop("Review package is incomplete.")
  inventory <- readr::read_tsv(file.path(root, required[[1L]]), show_col_types = FALSE)
  profile <- readr::read_tsv(file.path(root, "combined_profile_review.tsv"), show_col_types = FALSE)
  uv <- readr::read_tsv(file.path(root, "uv_specific_review.tsv"), show_col_types = FALSE)
  architecture <- readr::read_tsv(file.path(root, "condition_architecture_review.tsv"), show_col_types = FALSE)
  reps <- readr::read_tsv(file.path(root, "representative_feature_review.tsv"), show_col_types = FALSE)
  if (anyDuplicated(inventory$historical_rule_id) || nrow(profile) != 40L ||
      anyDuplicated(profile$gene_set_id) || nrow(uv) != 13L || nrow(architecture) != 6L ||
      anyDuplicated(architecture$condition)) stop("Review key coverage failed.")
  for (x in list(profile, architecture, reps)) {
    if (!all(c("author_decision", "author_notes") %in% names(x)) ||
        any(!is.na(x$author_decision)) || any(!is.na(x$author_notes))) {
      stop("Author fields must remain blank.")
    }
  }
  predicted <- unlist(profile[c("stability_region", "signal_class")])
  if (any(grepl("Instability|instable|weakly_stable|weakly stable", predicted, ignore.case = TRUE))) {
    stop("Obsolete current terminology detected.")
  }
  script <- readLines("scripts/current/build_current_biological_review_package.R", warn = FALSE)
  if (any(grepl("SSI *(>=|>|<) *0\\.60", script))) stop("Obsolete executable boundary detected.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  build_current_biological_review_package()
  audit_current_biological_review_package()
  cat("current biological review package: PASS\n")
}
