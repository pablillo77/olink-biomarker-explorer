# tests/testthat/test_qc.R
# Test the pure QC functions against a tiny hand-built long-format tibble,
# not npx_data1 — small, deterministic fixtures make failures easy to read.

library(testthat)

toy_npx <- tibble::tibble(
  SampleID   = c("S1", "S1", "S2", "S2", "CONTROL_SAMPLE_AS 1", "CONTROL_SAMPLE_AS 1"),
  OlinkID    = rep(c("OID001", "OID002"), 3),
  Assay      = rep(c("ProtA", "ProtB"), 3),
  Panel      = "Panel1",
  QC_Warning = c("Pass", "Warning", "Pass", "Pass", "Pass", "Pass"),
  NPX        = c(5, 1, 6, 7, 4, 4),
  LOD        = c(2, 2, 2, 2, 2, 2),
  MissingFreq = c(0.1, 0.1, 0.2, 0.2, 0, 0)
)

test_that("is_control_sample flags Olink control rows only", {
  flags <- is_control_sample(toy_npx)
  expect_equal(flags, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
})

test_that("sample_qc_summary computes correct warning and below-LOD rates", {
  s <- sample_qc_summary(toy_npx)
  s1 <- s[s$SampleID == "S1", ]
  expect_equal(s1$n_assays, 2)
  expect_equal(s1$n_warning, 1)
  expect_equal(s1$pct_warning, 0.5)
  expect_equal(s1$n_below_lod, 1)   # NPX values 5 and 1 vs LOD 2: only the 1 is below LOD
})

test_that("flag_failing_samples respects the threshold", {
  flagged <- flag_failing_samples(toy_npx, threshold = 0.4)
  expect_true("S1" %in% flagged)   # 50% warning rate > 40%
  expect_false("S2" %in% flagged)  # 0% warning rate
})

test_that("apply_sample_exclusions removes only the requested samples", {
  filtered <- apply_sample_exclusions(toy_npx, excluded = "S1")
  expect_false("S1" %in% filtered$SampleID)
  expect_true(all(c("S2", "CONTROL_SAMPLE_AS 1") %in% filtered$SampleID))
  expect_equal(nrow(filtered), 4)
})

test_that("missingness_overview returns sane summary stats", {
  m <- missingness_overview(toy_npx)
  expect_equal(m$n_na_npx, 0)
  expect_true(m$pct_below_lod >= 0 && m$pct_below_lod <= 1)
  expect_equal(round(m$mean_missfreq, 4), round(mean(toy_npx$MissingFreq), 4))
})

test_that("assay_qc_summary aggregates per assay, not per sample", {
  a <- assay_qc_summary(toy_npx)
  expect_equal(nrow(a), 2)  # OID001, OID002
  oid002 <- a[a$OlinkID == "OID002", ]
  expect_equal(oid002$n_warning, 1)  # only S1's OID002 row is "Warning"
})
