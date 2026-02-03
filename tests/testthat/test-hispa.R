context("HiSpaR basic functionality tests")

test_that("hispa_analyze handles invalid inputs", {
  # Test with non-numeric input
  expect_error(hispa_analyze("invalid", "output"))
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
    verbose = FALSE
  )
  
  # Check that output files exist
  expect_true(file.exists(file.path(output_dir, "final_positions.txt")))
  expect_true(file.exists(file.path(output_dir, "log_likelihood_trace.txt")))
  
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
