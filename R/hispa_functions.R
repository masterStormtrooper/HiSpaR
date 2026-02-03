#' Run HiSpa MCMC Analysis
#'
#' @description
#' Performs hierarchical Bayesian inference of 3D chromatin structure from
#' Hi-C contact matrix using MCMC sampling. This function follows the exact
#' workflow from the HiSpa C++ implementation.
#'
#' @param hic_experiment A Bioconductor `HiCExperiment` object or a numeric matrix
#'   containing Hi-C contact data. If a `HiCExperiment` object, the function extracts 
#'   the contact matrix using `gi2cm()` from the interactions and converts it to a 
#'   numeric matrix for analysis. If a matrix, it is used directly.
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
#' @param filter_quantile Numeric, quantile threshold for filtering loci by contact counts 
#'   (default: -1, no filtering). If >= 0, loci with row sums below this quantile are removed. 
#'   For example, 0.1 removes loci in the bottom 10\% of contact counts.
#' @param safe_mode Logical, whether to run analysis in a safe subprocess using callr 
#'   (default: TRUE). If TRUE and callr is available, analysis runs in an isolated subprocess 
#'   to prevent R crashes from affecting the parent session. Set to FALSE to run inline.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \strong{positions} - A numeric matrix of dimensions n x 3 containing the 
#'       final inferred 3D coordinates of loci. Rows correspond to loci, columns are (x, y, z).
#'     \item \strong{beta0} - The final intercept parameter of the log-distance relationship.
#'     \item \strong{beta1} - The final slope parameter of the log-distance relationship.
#'   }
#'   If filtering was applied, the positions matrix has an attribute `filtered_locus_indices` 
#'   containing the original indices of the retained loci (before filtering).
#'   Additionally, all analysis results are saved as text files in the output directory.
#'
#' @details
#' This function implements the complete HiSpa workflow:
#' \enumerate{
#'   \item Filter loci by contact count (optional)
#'   \item Load contact matrix from file
#'   \item Assign loci to clusters (k-means)
#'   \item Build cluster relationships
#'   \item Initialize structure (random or cluster-based)
#'   \item Assemble global structure
#'   \item Run main MCMC algorithm
#'   \item Save results to output directory
#' }
#' 
#' \bold{Locus Filtering:} By default, no filtering is applied (filter_quantile = -1). 
#' If filter_quantile >= 0, loci with contact counts below the specified quantile 
#' are removed before analysis. For example, filter_quantile = 0.1 removes loci in 
#' the bottom 10\% of contact counts. This improves computational efficiency and 
#' focuses analysis on loci with sufficient data.
#' 
#' All results are automatically saved to the output directory:
#' \itemize{
#'   \item \strong{final_positions.txt} - Final inferred 3D coordinates (n x 3 matrix)
#'   \item \strong{initial_positions.txt} - Initial positions before MCMC (n x 3 matrix)
#' }
#' 
#' Read results using standard R functions:
#' \code{final_pos <- read.table("output_dir/final_positions.txt")}
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
#' su_res <- hispa_analyze(
#'   hic_experiment = su1_contact_mat,
#'   output_dir = tempdir(),
#'   mcmc_iterations = 100,
#'   mcmc_burn_in = 10,
#'   use_cluster_init = FALSE,
#'   verbose = TRUE
#' )
#' 
#' # check results
#' dim(su_res$positions)  # Should be n x 3
#' }
#'
#' @importFrom utils write.table
#' @importFrom stats quantile
#' @importFrom HiCExperiment interactions gi2cm
#' @importFrom Matrix as.matrix
#' @export
hispa_analyze <- function(
    hic_experiment,
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
    verbose = TRUE,
    filter_quantile = -1,
    safe_mode = TRUE
) {
  # Handle both matrix and HiCExperiment inputs
  if (is.matrix(hic_experiment)) {
    # Input is already a matrix, use it directly
    contact_matrix <- hic_experiment
    
    # Ensure it's numeric and double precision
    if (!is.numeric(contact_matrix)) {
      stop("Input matrix must be numeric")
    }
    storage.mode(contact_matrix) <- "double"
    
  } else if (inherits(hic_experiment, "HiCExperiment")) {
    # Input is a HiCExperiment object, extract the contact matrix
    if (!requireNamespace("HiContacts", quietly = TRUE)) {
      stop("Package 'HiContacts' is required to extract contact matrices. Please install it.")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required for matrix conversion. Please install it.")
    }

    # Extract contact matrix using gi2cm from HiCExperiment
    # gi2cm returns a sparse matrix, so we need Matrix::as.matrix first, then base::as.matrix
    # cm <- Matrix::as.matrix(gi2cm(interactions(hic_experiment), 'count'))
    contact_matrix <- as.matrix(hic_experiment, use.scores = "count")
    
    storage.mode(contact_matrix) <- "double"
    
  } else {
    stop("`hic_experiment` must be either a numeric matrix or a HiCExperiment object")
  }
  
  if (nrow(contact_matrix) != ncol(contact_matrix)) {
    stop("contact matrix must be square")
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  # Filter loci by row sum if filter_quantile >= 0
  if (filter_quantile >= 0) {
    row_sums <- rowSums(contact_matrix)
    
    # Calculate threshold based on quantile
    threshold <- quantile(row_sums, probs = filter_quantile, names = FALSE)
    # Include rows with rowsum > threshold
    keep_idx <- row_sums > threshold
    
    n_before <- nrow(contact_matrix)
    n_after <- sum(keep_idx)
    n_removed <- n_before - n_after
    
    if (n_removed > 0) {
      if (verbose) {
        cat("Filtering loci by contact count:\n")
        cat("  Original number of loci: ", n_before, "\n")
        cat("  Threshold (", filter_quantile, " quantile): ", threshold, "\n")
        cat("  Loci removed: ", n_removed, "\n")
        cat("  Final number of loci: ", n_after, "\n")
      }
      
      # Filter the contact matrix
      contact_matrix <- contact_matrix[keep_idx, keep_idx]
      
      # Store original locus indices for later mapping
      original_indices <- which(keep_idx)
    } else {
      if (verbose) {
        cat("No loci filtered (all loci above threshold)\n")
      }
      original_indices <- seq_len(nrow(contact_matrix))
    }
  } else {
    # No filtering
    if (verbose) {
      cat("No locus filtering applied (filter_quantile = ", filter_quantile, ")\n")
    }
    original_indices <- seq_len(nrow(contact_matrix))
  }

  # Decide whether to run in subprocess
  use_subprocess <- isTRUE(safe_mode) && requireNamespace("callr", quietly = TRUE)
  if (isTRUE(safe_mode) && !use_subprocess) {
    message("safe_mode requested but package 'callr' not available; running inline. Install callr for crash isolation.")
  }

  # Prepare args for C++ call (now passing contact_matrix directly, not file path)
  if (use_subprocess) {
    result <- callr::r_safe(
      function(cm, outdir, iter, ncl, burn, init_sd, sd_fl, sd_ce, use_cl, cl_iter, cl_sd, samp, samp_int, verb) {
        pkg <- asNamespace("HiSpaR")
        fun <- get('.hispa_analyze_cpp', envir = pkg, inherits = FALSE)
        fun(cm, outdir, iter, ncl, burn, init_sd, sd_fl, sd_ce, use_cl, cl_iter, cl_sd, samp, samp_int, verb)
      },
      args = list(
        contact_matrix,
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
        as.logical(verbose)
      ),
      show = isTRUE(verbose)
    )
  } else {
    result <- .hispa_analyze_cpp(
      contact_matrix,
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
      as.logical(verbose)
    )
  }

  # Add original indices as attribute to the positions matrix in the result list
  if (!is.null(result$positions)) {
    attr(result$positions, "filtered_locus_indices") <- original_indices
  }
  
  result
}
