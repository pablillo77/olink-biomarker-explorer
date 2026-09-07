# R/data_access.R
# Data-access layer. Everything the app knows about *where data comes from*
# lives here, so the Shiny code never touches OlinkAnalyze datasets directly.
# Roxygen-style comments (#') document each function; they also feed
# devtools::document() if this ever becomes a package.

#' Load one of Olink's synthetic demo NPX datasets
#'
#' `npx_data1` and `npx_data2` ship inside OlinkAnalyze. They are randomly
#' generated example data, not real patient measurements.
#'
#' @param which "npx_data1" or "npx_data2".
#' @return A long-format tibble: one row per (SampleID, OlinkID).
load_npx_demo <- function(which = c("npx_data1", "npx_data2")) {
  which <- match.arg(which)
  # getExportedValue() fetches a lazy-loaded dataset by name without
  # attaching the whole package (keeps the app's namespace clean).
  df <- getExportedValue("OlinkAnalyze", which)
  tibble::as_tibble(df)
}

#' Read a user-uploaded NPX export (csv/xlsx/parquet/zip)
#'
#' Thin wrapper so the app has one place to change if the OlinkAnalyze
#' reader API changes again (it was rewritten in 5.0.0).
#'
#' @param path Path to the uploaded file.
#' @return Long-format tibble.
read_npx_upload <- function(path) {
  OlinkAnalyze::read_npx(filename = path)
}

#' Identify Olink control samples
#'
#' Olink external controls are named CONTROL_SAMPLE_* in the demo data.
#' Real exports may also carry a SampleType column; we honour it if present.
#'
#' @param df Long-format NPX tibble.
#' @return Logical vector, TRUE for control-sample rows.
is_control_sample <- function(df) {
  by_name <- grepl("^CONTROL", df$SampleID, ignore.case = TRUE)
  if ("SampleType" %in% names(df)) {
    by_type <- !is.na(df$SampleType) & df$SampleType != "SAMPLE"
    return(by_name | by_type)
  }
  by_name
}

#' Summarise the shape of an NPX dataset
#'
#' @param df Long-format NPX tibble.
#' @return Named list of counts used by the app's overview panel.
npx_summary <- function(df) {
  ctrl <- is_control_sample(df)
  list(
    rows            = nrow(df),
    cols            = ncol(df),
    n_samples       = dplyr::n_distinct(df$SampleID),
    n_study_samples = dplyr::n_distinct(df$SampleID[!ctrl]),
    n_control       = dplyr::n_distinct(df$SampleID[ctrl]),
    n_assays        = dplyr::n_distinct(df$OlinkID),
    n_panels        = dplyr::n_distinct(df$Panel),
    panels          = sort(unique(df$Panel)),
    qc_warning_rows = sum(df$QC_Warning == "Warning", na.rm = TRUE),
    pct_below_lod   = mean(df$NPX < df$LOD, na.rm = TRUE),
    n_na_npx        = sum(is.na(df$NPX))
  )
}
