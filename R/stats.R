# R/stats.R
# Analysis functions: PCA and differential abundance. Pure, no Shiny.
# Both take the raw long-format NPX tibble and an optional vector of
# SampleIDs to exclude (wired from the QC module), so every downstream
# tab respects whatever the user flagged in QC.

#' Run PCA on the samples x proteins matrix
#'
#' Pivots long -> wide, scales each protein (mean 0, sd 1), and runs
#' base R's prcomp(). Scaling matters here: NPX values already sit on
#' a comparable log2 scale across proteins, but some proteins are far
#' more variable than others, and PCA is variance-driven — without
#' scaling, high-variance proteins would dominate the components for
#' reasons that have nothing to do with biology.
#'
#' @param df Long-format NPX tibble.
#' @param excluded SampleIDs to drop before computing PCA.
#' @param n_components How many PCs to keep in the returned scores.
#' @return list(scores = tibble with metadata + PC1..PCn,
#'              variance_explained = numeric vector,
#'              loadings = tibble of protein loadings on PC1/PC2)
run_pca <- function(df, excluded = character(0), n_components = 5) {
  df <- apply_sample_exclusions(df, excluded)
  df <- df[!is_control_sample(df), ]

  meta <- df |>
    dplyr::distinct(.data$SampleID, .data$Subject, .data$Treatment,
                     .data$Site, .data$Time)

  wide <- df |>
    dplyr::select(.data$SampleID, .data$Assay, .data$NPX) |>
    tidyr::pivot_wider(names_from = .data$Assay, values_from = .data$NPX)

  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$SampleID
  # Drop any protein with zero variance (would blow up scaling).
  mat <- mat[, apply(mat, 2, stats::sd) > 0, drop = FALSE]

  pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  n_components <- min(n_components, ncol(pca$x))

  scores <- tibble::as_tibble(pca$x[, seq_len(n_components), drop = FALSE])
  scores$SampleID <- rownames(pca$x)
  scores <- dplyr::left_join(meta, scores, by = "SampleID")

  variance_explained <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(n_components)]

  loadings <- tibble::as_tibble(pca$rotation[, 1:2], rownames = "Assay") |>
    dplyr::rename(PC1_loading = "PC1", PC2_loading = "PC2")

  list(scores = scores, variance_explained = variance_explained, loadings = loadings)
}

#' Differential abundance: per-protein Treatment effect, accounting for
#' repeated measures within Subject via a linear mixed model.
#'
#' Each Subject contributes multiple samples (one per Time), so treating
#' all samples as independent (a plain t-test) understates the true
#' variance and inflates significance. We fit, per protein:
#'   NPX ~ Treatment + (1 | Subject)
#' via OlinkAnalyze::olink_lmer(), which loops this across every protein
#' and applies Benjamini-Hochberg correction. It returns an ANOVA-style
#' table (F-statistic via lmerTest) rather than a signed coefficient, so
#' it tells us HOW significant the Treatment effect is but not which
#' direction it points. We compute that direction separately as a plain
#' mean difference (Treated - Untreated) and join it in, for the volcano
#' plot's x-axis.
#'
#' Requires the `lmerTest` package (loaded implicitly by olink_lmer) in
#' addition to `lme4` — install both, or olink_lmer errors or produces
#' an unusable result.
#'
#' @param df Long-format NPX tibble.
#' @param excluded SampleIDs to drop before testing.
#' @return Tibble with one row per protein: Assay, estimate (mean diff),
#'   p.value, Adjusted_pval.
run_differential_abundance <- function(df, excluded = character(0)) {
  df <- apply_sample_exclusions(df, excluded)
  df <- df[!is_control_sample(df), ]

  lmer_res <- OlinkAnalyze::olink_lmer(
    df       = df,
    variable = "Treatment",
    random   = "Subject"
  )

  # Plain mean difference per protein — the effect's direction/magnitude,
  # independent of the mixed model's significance test.
  effect_size <- df |>
    dplyr::group_by(.data$Assay, .data$Treatment) |>
    dplyr::summarise(mean_npx = mean(.data$NPX, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$Treatment, values_from = .data$mean_npx) |>
    dplyr::mutate(estimate = .data$Treated - .data$Untreated) |>
    dplyr::select(.data$Assay, .data$estimate)

  lmer_res |>
    dplyr::left_join(effect_size, by = "Assay") |>
    dplyr::select(.data$Assay, .data$OlinkID, .data$Panel, .data$estimate,
                  .data$p.value, .data$Adjusted_pval) |>
    dplyr::arrange(.data$Adjusted_pval)
}
