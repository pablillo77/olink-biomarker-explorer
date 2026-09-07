# R/qc.R
# Pure QC functions: no Shiny here. Keeping analysis logic separate from UI
# means every function can be unit-tested with testthat without spinning up
# an app, and reused identically in the Quarto notebook.

#' Per-sample QC summary
#'
#' One row per SampleID with warning rate, below-LOD rate, and a control flag.
#' `pct_warning` is what we threshold on to flag failing samples.
#'
#' @param df Long-format NPX tibble.
#' @return Tibble, one row per SampleID.
sample_qc_summary <- function(df) {
  df |>
    dplyr::mutate(.control = is_control_sample(df)) |>
    dplyr::group_by(.data$SampleID, .data$.control) |>
    dplyr::summarise(
      n_assays      = dplyr::n(),
      n_warning     = sum(.data$QC_Warning == "Warning", na.rm = TRUE),
      pct_warning   = .data$n_warning / .data$n_assays,
      n_below_lod   = sum(.data$NPX < .data$LOD, na.rm = TRUE),
      pct_below_lod = .data$n_below_lod / .data$n_assays,
      .groups = "drop"
    ) |>
    dplyr::rename(is_control = .data$.control) |>
    dplyr::arrange(dplyr::desc(.data$pct_warning))
}

#' Per-assay (per-protein) QC summary
#'
#' @param df Long-format NPX tibble.
#' @return Tibble, one row per OlinkID/Assay.
assay_qc_summary <- function(df) {
  df |>
    dplyr::group_by(.data$OlinkID, .data$Assay, .data$Panel) |>
    dplyr::summarise(
      n_samples     = dplyr::n(),
      n_warning     = sum(.data$QC_Warning == "Warning", na.rm = TRUE),
      pct_warning   = .data$n_warning / .data$n_samples,
      n_below_lod   = sum(.data$NPX < .data$LOD, na.rm = TRUE),
      pct_below_lod = .data$n_below_lod / .data$n_samples,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$pct_below_lod))
}

#' Flag samples whose QC-warning rate exceeds a threshold
#'
#' @param df Long-format NPX tibble.
#' @param threshold Fraction of assays with "Warning" above which a sample
#'   is flagged (default 0.10 = 10%).
#' @return Character vector of SampleIDs to flag.
flag_failing_samples <- function(df, threshold = 0.10) {
  qc <- sample_qc_summary(df)
  qc$SampleID[qc$pct_warning > threshold]
}

#' Remove excluded samples from the dataset
#'
#' @param df Long-format NPX tibble.
#' @param excluded Character vector of SampleIDs to drop.
#' @return Filtered tibble.
apply_sample_exclusions <- function(df, excluded = character(0)) {
  dplyr::filter(df, !.data$SampleID %in% excluded)
}

#' Overall missingness overview
#'
#' Uses OlinkAnalyze's own `MissingFreq` column (fraction of samples where
#' that assay was missing, precomputed by Olink) plus our own below-LOD rate.
#'
#' @param df Long-format NPX tibble.
#' @return A one-row tibble of overall stats, for a value-box style summary.
missingness_overview <- function(df) {
  tibble::tibble(
    n_na_npx        = sum(is.na(df$NPX)),
    pct_na_npx      = mean(is.na(df$NPX)),
    pct_below_lod   = mean(df$NPX < df$LOD, na.rm = TRUE),
    mean_missfreq   = mean(df$MissingFreq, na.rm = TRUE),
    max_missfreq    = max(df$MissingFreq, na.rm = TRUE)
  )
}
