#' @title Create Meta-Cell Expression Matrix
#'
#' @description Generates a meta-cell expression matrix by aggregating individual cells' expression profiles
#' according to the specified meta-cell groupings.
#'
#' @param exp_matrix A gene expression matrix (genes x cells) where columns are cells and rows are genes
#' @param meta_cells_list A list where each element contains indices of cells belonging to one meta-cell
#'
#' @return A data frame representing the meta-cell expression matrix (genes x meta-cells),
#'         where each meta-cell's expression is the sum of its constituent cells' expression
#'
#' @export
#'
#' @examples
#' # Create example data
#' exp_mat <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
#' meta_list <- list(1:3, 4:6, 7:10)
#'
#' # Generate meta-cell matrix
#' meta_matrix <- create_meta_cell_matrix(exp_mat, meta_list)
create_meta_cell_matrix <- function(exp_matrix, meta_cells_list){
  # Input validation
  if (!is.matrix(exp_matrix)) {
    stop("exp_matrix must be a matrix")
  }
  if (!is.list(meta_cells_list)) {
    stop("meta_cells_list must be a list")
  }
  if (any(!unlist(lapply(meta_cells_list, is.numeric)))) {
    stop("All elements of meta_cells_list must be numeric vectors")
  }

  # Initialize result matrix
  n_meta <- length(meta_cells_list)
  meta_matrix <- matrix(
    nrow = nrow(exp_matrix),
    ncol = n_meta,
    dimnames = list(
      rownames(exp_matrix),
      paste0("MetaCell_", seq_len(n_meta))
    )
  )

    # Create each meta-cell profile
    for (i in seq_len(n_meta)) {
      cell_indices <- meta_cells_list[[i]]

      # Validate cell indices
      if (any(cell_indices < 1 | cell_indices > ncol(exp_matrix))) {
        stop(sprintf("Invalid cell indices in meta-cell %d", i))
      }

      # Extract and sum expression profiles
      cell_subset <- exp_matrix[, cell_indices, drop = FALSE]
      meta_matrix[, i] <- if (ncol(cell_subset) == 1) {
        cell_subset
      } else {
        rowSums(cell_subset)
      }
    }

    return(as.data.frame(meta_matrix))
}



