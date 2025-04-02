#' @title Create Library Size-Stabilized Meta-Cells Within Specified Groups
#' @description Generates meta-cells with stabilized library sizes separately for each group of cells,
#' then combines the results. This preserves biological heterogeneity while reducing
#' technical noise within each group.
#' @param count_matrix A gene expression count matrix (genes x cells).
#' @param similarity_matrix Optional pre-computed cell-cell similarity matrix (default: NULL).
#' @param meta_cells_num Target total number of meta-cells (default: 1000).
#' @param group_id Vector specifying group membership for each cell (length must match ncol(count_matrix))
#' @param seed Random seed for reproducibility (default: 123)
#' @return A combined meta-cell expression matrix where:
#' \itemize{
#'   \item Rows represent genes (same as input)
#'   \item Columns represent meta-cells with group prefixes
#'   \item Values are aggregated counts of constituent cells
#' }
#'
#' @export
#'
#' @examples
#' # Example usage:
#' data(sample_counts) # genes x cells matrix
#' groups <- rep(c("A","B"), each=100) # 200 cells
#' metacells <- lib_stabilized_metacells_by_group(
#'   count_matrix = as.matrix(sample_counts),
#'   meta_cells_num = 20,
#'   group_id = groups
#' )
lib_stabilized_metacells_by_group <- function(count_matrix, meta_cells_num=1000,  group_id, similarity_matrix=NULL, seed=123){

  # Input validation
  if (!is.matrix(count_matrix)) {
    stop("count_matrix must be a matrix")
  }
  if (length(group_id) != ncol(count_matrix)) {
    stop("group_id length must match number of cells in count_matrix")
  }
  if (meta_cells_num <= 0) {
    stop("meta_cells_num must be a positive integer")
  }
  if (!is.null(similarity_matrix) &&
      (nrow(similarity_matrix) != ncol(count_matrix) ||
       ncol(similarity_matrix) != ncol(count_matrix))) {
    stop("similarity_matrix dimensions must match count_matrix columns")
  }

  cells_seq_dep <- apply(count_matrix, 2, sum)
  seq_dep_sum <- sum(cells_seq_dep)
  group_seq_dep <- list()
  count_matrices <- list()
  group_id <- as.character(group_id)
  similarity_matrices <- list()
  unique_groups <- unique(group_id)

  # Split data by groups
  for (id in unique_groups){
    sub_data <- count_matrix[, group_id==id]
    count_matrices[[as.character(id)]] <- sub_data
    group_seq_dep[[as.character(id)]] <- sum(apply(sub_data, 2, sum))

    if (!is.null(similarity_matrix)){
      sub_similarity_matrix <- similarity_matrix[,group_id==id]
      similarity_matrices[[as.character(id)]] <- sub_similarity_matrix
    }
    else{
      similarity_matrices[[as.character(id)]] <- NULL
      }
  }



  # Calculate meta-cells per group proportional to sequencing depth
  sub_cells_num <- round((unlist(group_seq_dep)/seq_dep_sum) * meta_cells_num)

  # Process each group separately
  sub_output_list <-list()
  for (i in names(sub_cells_num)){
    message("Processing group: ", i," (", ncol(count_matrices[[i]]), " cells)")
    #if (sub_cells_num[i]==0){sub_cells_num[i] <- 1}
    tryCatch({
      sub_output <- lib_stabilized_metacells(
        count_matrices[[i]],
        similarity_matrices[[i]],
        meta_cells_num = sub_cells_num[i],
        seed = seed
      )
      colnames(sub_output) <- paste(i, colnames(sub_output), sep="_")
      sub_output_list[[i]] <- sub_output
    }, error = function(e) {
      warning("Failed to process group ", i, ": ", e$message)
    })
  }

  # Combine all group results
  output <- do.call(cbind, sub_output_list[!sapply(sub_output_list, is.null)])
  if (is.null(output)) {
    stop("No valid meta-cells were created from any group")
  }
  message("Successfully created ", ncol(output), " meta-cells across ",
          length(unique_groups), " groups")
  return(output)
}
