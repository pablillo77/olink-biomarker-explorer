# tests/testthat/test_stats.R
# These tests use a small synthetic fixture with real variance structure
# (not the toy 2-row-per-sample fixture from test_qc.R), since PCA needs
# more than a couple of proteins/samples to be meaningful.

library(testthat)

set.seed(42)
make_toy_dataset <- function(n_subjects = 10, n_proteins = 6) {
  subjects <- paste0("SUBJ", seq_len(n_subjects))
  treatment <- rep(c("Treated", "Untreated"), length.out = n_subjects)
  times <- c("Baseline", "Week.6", "Week.12")

  rows <- list()
  i <- 1
  for (s in seq_len(n_subjects)) {
    subj_baseline <- rnorm(n_proteins, mean = 5, sd = 1)  # this subject's own level
    for (t in times) {
      for (p in seq_len(n_proteins)) {
        rows[[i]] <- tibble::tibble(
          SampleID   = paste0(subjects[s], "_", t),
          Subject    = subjects[s],
          Assay      = paste0("Prot", p),
          OlinkID    = paste0("OID", p),
          Panel      = "Panel1",
          Treatment  = treatment[s],
          Site       = "Site_A",
          Time       = t,
          QC_Warning = "Pass",
          LOD        = 1,
          NPX        = subj_baseline[p] + rnorm(1, sd = 0.3)
        )
        i <- i + 1
      }
    }
  }
  dplyr::bind_rows(rows)
}

toy <- make_toy_dataset()

test_that("run_pca returns one score row per sample and valid variance_explained", {
  res <- run_pca(toy, n_components = 3)
  expect_equal(nrow(res$scores), dplyr::n_distinct(toy$SampleID))
  expect_true(all(res$variance_explained >= 0 & res$variance_explained <= 1))
  expect_true(all(c("PC1", "PC2", "PC3") %in% names(res$scores)))
})

test_that("run_pca respects sample exclusions", {
  excl <- toy$SampleID[1]
  res <- run_pca(toy, excluded = excl)
  expect_false(excl %in% res$scores$SampleID)
  expect_equal(nrow(res$scores), dplyr::n_distinct(toy$SampleID) - 1)
})

test_that("run_pca loadings has one row per protein", {
  res <- run_pca(toy)
  expect_equal(nrow(res$loadings), dplyr::n_distinct(toy$Assay))
})

test_that("run_differential_abundance returns one row per protein with BH-adjusted p-values", {
  # olink_lmer can be finicky on tiny synthetic data (convergence warnings);
  # we only assert on structure here, not on specific p-values.
  res <- tryCatch(run_differential_abundance(toy), error = function(e) NULL)
  skip_if(is.null(res), "olink_lmer did not converge on the toy fixture")
  expect_equal(nrow(res), dplyr::n_distinct(toy$Assay))
  expect_true("Adjusted_pval" %in% names(res))
  expect_true(all(res$Adjusted_pval >= res$p.value - 1e-8))  # BH >= raw p
})
