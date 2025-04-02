#' @title Identifiy Cell Indices for Meta-Cell Construction
#' @description This function groups individual cells into meta-cells based on sequencing depth and
#' cell-cell similarity matrix, aiming to create meta-cells with balanced library sizes.
#' @param seq_dep Numeric vector of sequencing depths (library sizes) for each cell.
#' @param similarity_matrix Square similarity probability between cells.
#' @param meta_cells_num Integer specifying the target number meta-cells to create.
#' @param seed Integer for random seed. Set for reproducible sampling when creating meta-cells.
#' @return A list where each element contains the indices of cells belonging to one meta-cell.
#' @export
#'
#' @examples
#' #Example usage:
#' seq_depth <- c(1000,1501,800,1200,2000,900,100,700)
#' similarity <- matrix(runif(64), nrow=8)
#' meta_indices <- meta_cells_index(seq_depth, similarity, 3)
meta_cells_index <- function(seq_dep, similarity_matrix, meta_cells_num, seed=123){
  set.seed(seed)
  # Input validation
  if (!is.numeric(seq_dep)) stop("seq_dep must be a numeric vector")
  if (!is.matrix(similarity_matrix) || nrow(similarity_matrix) != ncol(similarity_matrix)) {
    stop("similarity_matrix must be a square matrix")
  }
  if (length(seq_dep) != nrow(similarity_matrix)) {
    stop("Length of seq_dep must match dimension of similarity_matrix")
  }
  if (!is.numeric(meta_cells_num) || meta_cells_num <= 0) {
    stop("meta_cells_num must be a positive integer")
  }

  total_sets <- seq_along(seq_dep)
  target_lib_size <- sum(seq_dep)/meta_cells_num
  meta_cells <- list()
  iteration=0

  while (any(seq_dep<Inf)){
    iteration = iteration + 1

    # Start with cell having smallest remaining library size
    start_node <- which.min(seq_dep)

    if (seq_dep[start_node] < target_lib_size){
      last_depth <- 0
      current_depth <- seq_dep[start_node]
      current_meta_cell <- c(start_node)

       while (current_depth < target_lib_size){
        next_probs <- calculate_affinity_probabilities(total_sets, current_meta_cell, similarity_matrix)

        if (sum(is.nan(next_probs))>0) {
          break
        }

        # Sample next cell weighted by affinity probabilities
        #set.seed(seed)
        new_node <- sample.int(length(next_probs), size=1, prob=next_probs)

        last_depth <- current_depth
        current_depth <- current_depth + seq_dep[new_node]
        current_meta_cell <- c(current_meta_cell, new_node)
      }
      if (abs(current_depth - target_lib_size) > abs(last_depth - target_lib_size)){
        current_meta_cell <- current_meta_cell[-length(current_meta_cell)]
      }
    }else{
      current_meta_cell <- start_node
    }

    # Store results and update tracking variables
    meta_cells[[iteration]] <- current_meta_cell

    # Remove assigned cells from further consideration
    similarity_matrix[current_meta_cell,] <- 0
    similarity_matrix[,current_meta_cell] <- 0
    seq_dep[current_meta_cell] <- Inf
  }

  return(meta_cells)
}



