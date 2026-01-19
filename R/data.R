#' Example Hi-C Contact Matrix
#'
#' @description
#' Hi-C contact matrix from Drosophila melanogaster chromosome SU1
#' containing normalized contact frequencies between genomic loci.
#' This dataset is provided as an example for testing and demonstrating
#' the HiSpaR package functionality.
#'
#' @format A symmetric numeric matrix with 649 rows and 649 columns, where 
#'   entry (i,j) represents the normalized contact frequency between 
#'   genomic loci i and j. The matrix is:
#'   \describe{
#'     \item{Dimensions}{649 x 649}
#'     \item{Type}{Numeric matrix}
#'     \item{Symmetry}{Symmetric (contact_matrix[i,j] = contact_matrix[j,i])}
#'     \item{Diagonal}{Zero or near-zero (self-contacts)}
#'     \item{Range}{Non-negative contact frequencies}
#'   }
#'
#' @details
#' This contact matrix represents chromatin interaction frequencies derived
#' from Hi-C experiments on Drosophila melanogaster chromosome SU1. Higher
#' values indicate more frequent spatial proximity between genomic loci in
#' the 3D nuclear space.
#' 
#' The data can be used directly with \code{\link{hispa_analyze}} by first
#' saving it to a text file.
#'
#' @source Drosophila melanogaster Hi-C data, chromosome SU1
#'
#' @examples
#' # Load the example data
#' data(su1_contact_mat)
#' 
#' # Check dimensions
#' dim(su1_contact_mat)
#' 
#' # Check if matrix is symmetric
#' isSymmetric(su1_contact_mat)
#' 
#' # Summary statistics
#' summary(as.vector(su1_contact_mat))
#' 
#' \donttest{
#' # Visualize contact matrix
#' image(su1_contact_mat, 
#'       main = "Hi-C Contact Matrix (SU1)",
#'       xlab = "Genomic Locus", 
#'       ylab = "Genomic Locus")
#' }
#'
#' @seealso \code{\link{hispa_analyze}} for running the analysis
#' @keywords datasets
"su1_contact_mat"