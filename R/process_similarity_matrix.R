#' @title Pre-process Similarity Matrix
#'
#' @description
#' Prepares a similarity matrix for downstream analysis by ensuring non-negative values
#' and zero diagonal. Negative values are replaced with small positive values (0.001 * minimum absolute value).
#'
#' @param similarity_matrix A square numeric matrix representing pairwise similaritys.
#'        Can contain negative values which will be automatically handled.
#'
#' @return A processed probability matrix with:
#' \itemize{
#'   \item Diagonal elements set to 0
#'   \item Negative values replaced with 0.001 * minimum absolute value
#' }
#'
#' @export
#'
#' @examples
#' # Create a sample similarity matrix
#' connect_mat <- matrix(c(0, 1.2, -0.5,
#'                     1.2, 0, 0.8,
#'                     -0.5, 0.8, 0), nrow = 3)
#'
#' # Process the matrix
#' processed_mat <- process_similarity_matrix(connect_mat)
process_similarity_matrix <- function(similarity_matrix) {
  # Input validation
  if (!is.matrix(similarity_matrix)) {
    stop("Input must be a matrix")
  }

  if (nrow(similarity_matrix) != ncol(similarity_matrix)) {
    stop("Similarity matrix must be square")
  }

  if (!is.numeric(similarity_matrix)) {
    stop("Similarity matrix must contain numeric values")
  }


  # Set diagonal to 0
  diag(similarity_matrix) <- 0

  # Handle negative values
  if (any(similarity_matrix < 0, na.rm = TRUE)) {
    min_abs_val <- min(abs(similarity_matrix[similarity_matrix != 0]), na.rm = TRUE)
    replacement_val <- 0.001 * min_abs_val

    warning(sprintf(
      "%d negative values found in similarity matrix. Replacing with %.2e",
      sum(similarity_matrix < 0, na.rm = TRUE),
      replacement_val
    ))

    similarity_matrix[similarity_matrix < 0] <- replacement_val
  }

  return(similarity_matrix)
}
