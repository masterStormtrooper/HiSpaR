#' HiSpaR: Hierarchical Inference of Spatial Positions from Hi-C Data
#'
#' @description
#' HiSpaR provides R bindings for the HiSpa C++ library, enabling hierarchical
#' Bayesian inference of 3D chromatin structures from Hi-C contact matrices.
#'
#' @details
#' The package provides functions to:
#' \itemize{
#'   \item Analyze Hi-C contact matrices using MCMC sampling
#'   \item Infer 3D chromatin structures
#'   \item Use prior information from other datasets
#'   \item Perform cluster-based hierarchical analysis
#' }
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{hispa_analyze}}: Run complete HiSpa analysis
#' }
#'
#' @docType package
#' @name HiSpaR-package
#' @aliases HiSpaR
#' @useDynLib HiSpaR, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
