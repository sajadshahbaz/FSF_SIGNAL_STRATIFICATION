#' Assign FSF directional states
#'
#' Contract stub for assigning Up, Down, or Constant to finite signed effects.
#' The current scientific lock defines Up as `effect > tau`, Down as
#' `effect < -tau`, and Constant as `-tau <= effect <= tau`.
#'
#' @param data A data frame containing `feature_id`, `perturbation_id`, and
#'   `effect` columns.
#' @param tau One positive finite numeric scalar. The default is `0.5`.
#'
#' @return When implemented, the validated input with a `state` column whose
#'   values are `up`, `down`, or `constant`.
#' @export
fsf_assign_states <- function(data, tau = 0.5) {
  stop("Not implemented: Phase 2A contract stub", call. = FALSE)
}

#' Calculate FSF Signal Identity
#'
#' Contract stub for aggregating feature-perturbation states into state counts
#' and the Signal Identity `(p_up, p_down, p_const)` for each feature.
#'
#' @param data A validated data frame of feature-perturbation states or signed
#'   effects accepted by the future implementation.
#' @param tau One positive finite numeric scalar used only when `data` contains
#'   signed effects. The default is `0.5`.
#' @param group_cols `NULL`, `character(0)`, or a character vector naming
#'   columns that define independent grouping strata. `character(0)` is
#'   normalized to `NULL`.
#'
#' @return When implemented, one row per Signal Identity. When grouping is
#'   requested, grouping columns occur first in user-supplied order, followed
#'   by `feature_id`, counts, and probabilities.
#' @export
fsf_signal_identity <- function(data, tau = 0.5, group_cols = NULL) {
  stop("Not implemented: Phase 2A contract stub", call. = FALSE)
}

#' Classify FSF Signal Identities
#'
#' Contract stub for deriving dominant state, SSI, current stability region,
#' current signal class, and auxiliary Stability Deviation from Signal Identity.
#'
#' @param identity A data frame containing `feature_id`, `p_up`, `p_down`, and
#'   `p_const`, with optional state-count columns.
#'
#' @return When implemented, the identity columns plus `dominant_state`, `ssi`,
#'   `stability_region`, `signal_class`, and `stability_deviation`.
#' @export
fsf_classify <- function(identity) {
  stop("Not implemented: Phase 2A contract stub", call. = FALSE)
}

#' Run the complete FSF feature-level analysis
#'
#' Contract stub for the generic composition of state assignment, Signal
#' Identity calculation, and current-lock classification. It does not perform
#' differential-expression analysis, feature selection, or feature ranking.
#'
#' @param data A long-form data frame containing `feature_id`,
#'   `perturbation_id`, and finite numeric `effect` columns.
#' @param tau One positive finite numeric scalar. The default is `0.5`.
#' @param group_cols `NULL`, `character(0)`, or a character vector naming
#'   columns that define independent grouping strata. Signal Identity is
#'   calculated within `group_cols` plus `feature_id`.
#'
#' @return When implemented, one row per Signal Identity following the stable
#'   feature-level output contract, with grouping columns first when supplied.
#' @export
fsf_analyze <- function(data, tau = 0.5, group_cols = NULL) {
  stop("Not implemented: Phase 2A contract stub", call. = FALSE)
}

#' Calculate condition-level FSF signal architecture
#'
#' Contract stub for counting all ten current signal classes within each
#' condition and calculating directly comparable class proportions.
#'
#' @param data A classified feature data frame containing a condition column
#'   and `signal_class`.
#' @param condition_col A single character string naming the condition column.
#'
#' @return When implemented, all ten signal classes for every condition in
#'   deterministic order, including zero-count classes. The first column keeps
#'   the name supplied by `condition_col`, followed by `signal_class`,
#'   `class_count`, and `class_proportion`.
#' @export
fsf_architecture <- function(data, condition_col = "condition") {
  stop("Not implemented: Phase 2A contract stub", call. = FALSE)
}
