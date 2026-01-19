#' Run HiSpa MCMC Analysis
#'
#' @description
#' Performs hierarchical Bayesian inference of 3D chromatin structure from
#' Hi-C contact matrix using MCMC sampling. This function follows the exact
#' workflow from the HiSpa C++ implementation.
#'
#' @param contact_matrix Numeric matrix, Hi-C contact frequencies (symmetric, n x n).
#'   Alternatively, a character string path to a contact matrix file.
#' @param output_dir Character string specifying the output directory path.
#' @param mcmc_iterations Integer, number of MCMC iterations (default: 6000).
#' @param num_clusters Integer, number of clusters for hierarchical analysis. 
#'   If 0 (default), automatically determined as sqrt(n).
#' @param mcmc_burn_in Integer, number of burn-in iterations to discard (default: 0).
#' @param mcmc_initial_sd Numeric, initial standard deviation for MCMC proposals (default: 0.1).
#' @param mcmc_sd_floor Numeric, minimum allowed standard deviation (default: 0.0001).
#' @param mcmc_sd_ceil Numeric, maximum allowed standard deviation (default: 0.3).
#' @param use_cluster_init Logical, use cluster-based initialization instead of 
#'   random initialization (default: FALSE).
#' @param cluster_init_iterations Integer, number of iterations for cluster 
#'   initialization MCMC (default: 1000).
#' @param cluster_initial_sd Numeric, initial standard deviation for cluster 
#'   initialization MCMC (default: 0.1).
#' @param save_samples Logical, whether to save MCMC trace samples (default: FALSE).
#' @param sample_interval Integer, save samples every n iterations (default: 50).
#' @param verbose Logical, enable verbose output (default: TRUE).
#'
#' @return Invisibly returns TRUE if successful. All analysis results
#'   are saved as text files in the output directory (see Details).
#'
#' @details
#' This function implements the complete HiSpa workflow:
#' \enumerate{
#'   \item Load contact matrix from file
#'   \item Assign loci to clusters (k-means)
#'   \item Build cluster relationships
#'   \item Initialize structure (random or cluster-based)
#'   \item Assemble global structure
#'   \item Run main MCMC algorithm
#'   \item Save results to output directory
#' }
#' 
#' All results are automatically saved to the output directory:
#' \itemize{
#'   \item \strong{final_positions.txt} - Final inferred 3D coordinates (n x 3 matrix)
#'   \item \strong{initial_positions.txt} - Initial positions before MCMC (n x 3 matrix)
#'   \item \strong{log_likelihood_trace.txt} - MCMC log-likelihood values (convergence diagnostic)
#'   \item \strong{block_timings.txt} - Computation time for each MCMC block
#'   \item \strong{mcmc_log.txt} - Detailed analysis log with parameter values and settings
#' }
#' 
#' Read results using standard R functions:
#' \code{final_pos <- read.table("output_dir/final_positions.txt")}
#' \code{ll_trace <- scan("output_dir/log_likelihood_trace.txt")}
#'
#' @examples
#' # Load example contact matrix
#' data(su1_contact_mat)
#' 
#' # Check dimensions
#' dim(su1_contact_mat)
#' 
#' \donttest{
#' # Example 1: Run analysis with matrix input
#' hispa_analyze(
#'   contact_matrix = su1_contact_mat,
#'   output_dir = tempdir(),
#'   mcmc_iterations = 100,
#'   mcmc_burn_in = 10,
#'   use_cluster_init = FALSE,
#'   verbose = TRUE
#' )
#' 
#' # Read results from output directory
#' final_pos <- as.matrix(read.table(file.path(tempdir(), "final_positions.txt")))
#' head(final_pos)
#' }
#'
#' @importFrom utils write.table
#' @export
hispa_analyze <- function(
    contact_matrix,
    output_dir,
    mcmc_iterations = 6000L,
    num_clusters = 0L,
    mcmc_burn_in = 0L,
    mcmc_initial_sd = 0.1,
    mcmc_sd_floor = 0.0001,
    mcmc_sd_ceil = 0.3,
    use_cluster_init = FALSE,
    cluster_init_iterations = 1000L,
    cluster_initial_sd = 0.1,
    save_samples = FALSE,
    sample_interval = 50L,
    verbose = TRUE
) {
  # Handle input: either matrix or file path
  temp_file_created <- FALSE
  input_file <- NULL
  
  if (is.character(contact_matrix)) {
    # Input is a file path
    input_file <- contact_matrix
    if (!file.exists(input_file)) {
      stop("Input file does not exist: ", input_file)
    }
  } else if (is.matrix(contact_matrix) || is.data.frame(contact_matrix)) {
    # Input is a matrix - save to temporary file
    contact_matrix <- as.matrix(contact_matrix)
    
    # Validate matrix
    if (nrow(contact_matrix) != ncol(contact_matrix)) {
      stop("contact_matrix must be square (n x n)")
    }
    
    # Create temporary file
    input_file <- tempfile(pattern = "hispa_input_", fileext = ".txt")
    write.table(contact_matrix, input_file, 
                row.names = FALSE, col.names = FALSE)
    temp_file_created <- TRUE
    
    if (verbose) {
      message("Saved contact matrix (", nrow(contact_matrix), " x ", 
              ncol(contact_matrix), ") to temporary file")
    }
  } else {
    stop("contact_matrix must be either a matrix or a file path (character string)")
  }
  
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Ensure cleanup happens even if there's an error
  tryCatch({
    # Call C++ function (returns TRUE on success)
    result <- .hispa_analyze_cpp(
                    input_file,
                    output_dir,
                    as.integer(mcmc_iterations),
                    as.integer(num_clusters),
                    as.integer(mcmc_burn_in),
                    as.numeric(mcmc_initial_sd),
                    as.numeric(mcmc_sd_floor),
                    as.numeric(mcmc_sd_ceil),
                    as.logical(use_cluster_init),
                    as.integer(cluster_init_iterations),
                    as.numeric(cluster_initial_sd),
                    as.logical(save_samples),
                    as.integer(sample_interval),
                    as.logical(verbose))
    
    # Return success status invisibly
    invisible(result)
  }, finally = {
    # Clean up temporary file if we created one
    if (temp_file_created && file.exists(input_file)) {
      unlink(input_file)
    }
  })
}
