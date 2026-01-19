context("HiSpaR basic functionality tests")

test_that("hispa_analyze handles invalid inputs", {
  # Test with non-existent file
  expect_error(hispa_analyze("nonexistent.txt", "output"))
})

test_that("hispa_analyze runs successfully", {
  skip_if_not_installed("HiSpaR")
  skip_on_cran()
  
  # Create small test matrix and save to file
  n <- 20
  mat <- matrix(rpois(n*n, 10), n, n)
  mat <- (mat + t(mat)) / 2
  
  # Save to temporary file
  input_file <- file.path(tempdir(), "test_matrix.txt")
  write.table(mat, input_file, row.names = FALSE, col.names = FALSE)
  
  output_dir <- file.path(tempdir(), "hispa_test_output")
  dir.create(output_dir, showWarnings = FALSE)
  
  # Run analysis - returns output directory path
  result <- hispa_analyze(
    input_file = input_file,
    output_dir = output_dir,
    mcmc_iterations = 100,
    mcmc_burn_in = 10,
    verbose = FALSE
  )
  
  # Check that output files exist
  expect_true(file.exists(file.path(output_dir, "final_positions.txt")))
  expect_true(file.exists(file.path(output_dir, "log_likelihood_trace.txt")))
  
  # Clean up
  unlink(input_file)
  unlink(output_dir, recursive = TRUE)
  
  # Check dimensions
  expect_equal(nrow(result$position_matrix), n)
  expect_equal(ncol(result$position_matrix), 3)
  
  # Check class
  expect_s3_class(result, "hispa_result")
})
