#' @title  Calculate affinity probability matrix
#' @description
#' Computes the affinity probabilities from current metacell (group of cells) to other nodes based on
#' their average similarity.
#' @param total_sets Character or numeric vector of all cell identifiers in the network.
#' @param current_sets Character or numeric vector of cell identifiers in the current metacell.
#' @param similarity_matrix Square numeric matrix representing pairwise similaritys between cells.
#' Row and column names should match node identifiers in `total_sets`.
#' @return A numeric vector of transition probabilities where:
#' \itemize{
#'  \item Nodes in the `current_sets` have affinity 0
#'  \item Other nodes have affinity probabilites proportional to their average similarity to current nodes
#' }
#' @export
#'
#' @examples
#' # Create example data
#' all_cells <- paste0("Cell", 1:100)
#' current_metacell <- c("Cell1", "Cell5", "Cell10")
#' sim_mat <- matrix(runif(10000), nrow = 100, dimnames = list(all_cells, all_cells))
#'
#' # Calculate affinity probabilities
#' aff_probs <- calculate_affinity_probabilities(
#'   total_sets = all_cells,
#'   current_sets = current_metacell,
#'   similarity_matrix = sim_mat
#' )
calculate_affinity_probabilities <- function(total_sets, current_sets, similarity_matrix) {
  # Input validation
  if (!all(current_sets %in% total_sets)) {
    stop("All current_sets must be present in total_sets")
  }

  if (!identical(rownames(similarity_matrix), colnames(similarity_matrix))) {
    stop("Similarity matrix must be symmetric (identical row and column names)")
  }

  if (any(is.na(similarity_matrix))) {
    stop("Similarity matrix contains NA values")
  }

  # Initialize result vector with 0 for current metacell members
  affinities <- setNames(rep(0, length(total_sets)), total_sets)

  # Identify target cells (not in current metacell)
  target_cells <- setdiff(total_sets, current_sets)

  if (length(target_cells) == 0) {
    return(affinities)  # All cells are in current metacell
  }

  # Calculate mean similarities for target cells
  mean_similarities <- vapply(target_cells, function(cell) {
    mean(similarity_matrix[cell, current_sets, drop = FALSE])
  }, numeric(1))

  # Normalize to get probabilities
  affinities[target_cells] <- mean_similarities / sum(mean_similarities)

  return(affinities)
}
