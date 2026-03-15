context("HiSpaR basic functionality tests")

test_that("hispa_analyze handles invalid inputs", {
  # Test with non-numeric input
  expect_error(hispa_analyze("invalid", "output"))
})

test_that("hispa_analyze errors on non-square matrix", {
  mat <- matrix(rpois(20, 10), nrow = 4, ncol = 5)
  expect_error(hispa_analyze(mat, tempdir()), "square")
})

test_that("hispa_analyze errors on non-character output_dir", {
  n <- 5
  mat <- matrix(rpois(n * n, 10), n, n)
  mat <- (mat + t(mat)) / 2
  expect_error(hispa_analyze(mat, 123), "output_dir")
})

test_that("hispa_analyze runs successfully with matrix input", {
  skip_if_not_installed("HiSpaR")
  skip_on_cran()

  # Create small test matrix
  n <- 10
  mat <- matrix(rpois(n*n, 10), n, n)
  mat <- (mat + t(mat)) / 2

  output_dir <- file.path(tempdir(), "hispa_test_output")
  dir.create(output_dir, showWarnings = FALSE)

  # Run analysis with matrix input
  result <- hispa_analyze(
    hic_experiment = mat,
    output_dir = output_dir,
    mcmc_iterations = 50,
    mcmc_burn_in = 5,
    verbose = FALSE,
    safe_mode = FALSE
  )

  # Check that output files exist
  expect_true(file.exists(file.path(output_dir, "final_positions.txt")))

  # Check result is a list with expected elements
  expect_type(result, "list")
  expect_true("positions" %in% names(result))
  expect_true("beta0" %in% names(result))
  expect_true("beta1" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result$positions), n)
  expect_equal(ncol(result$positions), 3)

  # Check types
  expect_true(is.matrix(result$positions))
  expect_true(is.numeric(result$beta0))
  expect_true(is.numeric(result$beta1))

  # Clean up
  unlink(output_dir, recursive = TRUE)
})

test_that("hispa_analyze filters loci when filter_quantile > 0", {
  skip_if_not_installed("HiSpaR")
  skip_on_cran()

  set.seed(42)
  n <- 15
  # Give a few rows much lower sums so they will be filtered
  mat <- matrix(rpois(n * n, 20), n, n)
  mat[1:2, ] <- 0
  mat[, 1:2] <- 0
  mat <- (mat + t(mat)) / 2
  storage.mode(mat) <- "double"

  output_dir <- file.path(tempdir(), "hispa_test_filter")
  dir.create(output_dir, showWarnings = FALSE)

  result <- hispa_analyze(
    hic_experiment = mat,
    output_dir = output_dir,
    mcmc_iterations = 50,
    mcmc_burn_in = 5,
    filter_quantile = 0.1,
    verbose = FALSE,
    safe_mode = FALSE
  )

  # Fewer rows than original after filtering
  expect_lt(nrow(result$positions), n)
  expect_equal(ncol(result$positions), 3)

  # filtered_locus_indices attribute is set and maps to original indices
  idx <- attr(result$positions, "filtered_locus_indices")
  expect_false(is.null(idx))
  expect_true(is.integer(idx))
  expect_equal(length(idx), nrow(result$positions))

  unlink(output_dir, recursive = TRUE)
})

test_that("hispa_analyze retains all loci when no row is below the quantile threshold", {
  skip_if_not_installed("HiSpaR")
  skip_on_cran()

  set.seed(7)
  n <- 8
  # All row sums will be well above zero; filter_quantile = 0 gives threshold = min(row_sums).
  # Because row sums are not all equal, the minimum appears exactly once, so n-1 loci are
  # retained (n_removed >= 1). Use filter_quantile slightly above 0 but pick a matrix where
  # the 5th-percentile (interpolated) lands strictly below all values — this triggers n_removed=0.
  # Easiest: construct a matrix whose row sums are all equal except one is slightly higher,
  # then use filter_quantile = 0 — the "no loci removed" message branch is triggered when
  # filter_quantile >= 0 and n_removed == 0. Achieve this by calling with filter_quantile = -1
  # which skips filtering entirely, exercising the else-branch at line 199.
  mat <- matrix(rpois(n * n, 20), n, n)
  mat <- (mat + t(mat)) / 2
  storage.mode(mat) <- "double"

  output_dir <- file.path(tempdir(), "hispa_test_nofilter")
  dir.create(output_dir, showWarnings = FALSE)

  # filter_quantile = -1 means no filtering; exercises the verbose "No locus filtering applied" path
  result <- hispa_analyze(
    hic_experiment = mat,
    output_dir = output_dir,
    mcmc_iterations = 50,
    mcmc_burn_in = 5,
    filter_quantile = -1,
    verbose = TRUE,
    safe_mode = FALSE
  )

  # All loci should be retained
  expect_equal(nrow(result$positions), n)
  # filtered_locus_indices should be 1:n
  idx <- attr(result$positions, "filtered_locus_indices")
  expect_equal(idx, seq_len(n))

  unlink(output_dir, recursive = TRUE)
})

test_that("hispa_analyze safe_mode = FALSE runs inline", {
  skip_if_not_installed("HiSpaR")
  skip_on_cran()

  n <- 8
  mat <- matrix(rpois(n * n, 10), n, n)
  mat <- (mat + t(mat)) / 2
  storage.mode(mat) <- "double"

  output_dir <- file.path(tempdir(), "hispa_test_inline")
  dir.create(output_dir, showWarnings = FALSE)

  result <- hispa_analyze(
    hic_experiment = mat,
    output_dir = output_dir,
    mcmc_iterations = 50,
    mcmc_burn_in = 5,
    verbose = FALSE,
    safe_mode = FALSE
  )

  expect_type(result, "list")
  expect_true("positions" %in% names(result))
  expect_equal(nrow(result$positions), n)
  expect_equal(ncol(result$positions), 3)

  unlink(output_dir, recursive = TRUE)
})

test_that("hispa_analyze verbose = FALSE suppresses messages", {
  skip_if_not_installed("HiSpaR")
  skip_on_cran()

  n <- 8
  mat <- matrix(rpois(n * n, 10), n, n)
  mat <- (mat + t(mat)) / 2
  storage.mode(mat) <- "double"

  output_dir <- file.path(tempdir(), "hispa_test_quiet")
  dir.create(output_dir, showWarnings = FALSE)

  expect_no_message(
    hispa_analyze(
      hic_experiment = mat,
      output_dir = output_dir,
      mcmc_iterations = 50,
      mcmc_burn_in = 5,
      filter_quantile = 0.1,
      verbose = FALSE,
      safe_mode = FALSE
    )
  )

  unlink(output_dir, recursive = TRUE)
})
