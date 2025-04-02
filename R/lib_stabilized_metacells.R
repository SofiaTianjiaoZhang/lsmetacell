#' @title Create Library Size-Stabilized Meta-Cells
#' @description Generate meta-cells with stabilized library sizes by aggregating similar cells while
#' accounting for sequencing depth variation. This helps reduce technical noise in single-cell RNA-seq data.
#' @param count_matrix The gene expression count matrix (genes  x cells)
#' @param similarity_matrix Optional pre-computed cell-cell similarity matrix (default: NULL)
#' @param meta_cells_num Target number of meta-cells to create (default: 1000)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return A meta-cell expression matrix where:
#' \itemize{
#'   \item Rows represent genes (same as input)
#'   \item Columns represent meta-cells
#'   \item Values are aggregated counts of constituent cells
#' }
#' @export
#'
#'@importFrom WGCNA `cor`
#'
#' @examples
#' # Example usage:
#' data(sample_counts) # A genes x cells count matrix
#' metacells <- lib_stabilized_metacells(as.matrix(sample_counts), meta_cells_num = 5)
lib_stabilized_metacells <- function(count_matrix, similarity_matrix=NULL, meta_cells_num=1000, seed=123){
  # Input validation
  if (!is.matrix(count_matrix)) {
    stop("count_matrix must be a matrix")
  }
  if (any(count_matrix < 0)) {
    stop("count_matrix contains negative values")
  }
  if (meta_cells_num <= 0 || meta_cells_num > ncol(count_matrix)) {
    stop("meta_cells_num must be between 1 and number of cells")
  }

  # Set random seed if provided
  if (!is.numeric(seed)) {
    stop("seed must be an number")
  }


  # Create similarity matrix if not provided
   if (is.null(similarity_matrix)){
     message("Computing cell-cell similaritys by pearson correlation...")
    similarity_matrix <- process_similarity_matrix(cor(log2(count_matrix+0.00001)))
   }

  if (!is.null(similarity_matrix) &&
      (nrow(similarity_matrix) != ncol(count_matrix) ||
       ncol(similarity_matrix) != ncol(count_matrix))) {
    stop("similarity_matrix dimensions must match count_matrix columns")
  }

  # Calculate transition probabilities
  message("Calculating cell affinities...")
  similarity_matrix <- process_similarity_matrix(similarity_matrix)

  seq_dep <- apply(count_matrix, 2, sum)
  mean_cells <- ncol(count_matrix)/meta_cells_num

  #Filter outlier cells with extremely high library sizes
  filtered <- seq_dep < (median(seq_dep) * mean_cells * 1.5)
  if (sum(filtered) < ncol(count_matrix)) {
    message("Filtered out ", ncol(count_matrix) - sum(filtered), " outlier cells")
    count_matrix <- count_matrix[, filtered, drop = FALSE]
    similarity_matrix <- similarity_matrix[filtered, filtered, drop = FALSE]
    seq_dep <- seq_dep[filtered]
  }

  #Construct meta-cells
  message("Building meta-cells...")
  meta_index <- meta_cells_index(seq_dep, similarity_matrix, meta_cells_num=meta_cells_num, seed=seed)

  # Create final meta-cell matrix
  meta_cells <- create_meta_cell_matrix(count_matrix, meta_index)
  rownames(meta_cells) <- rownames(count_matrix)
  message(paste("Successfully created", ncol(meta_cells), "meta-cells"))

  return(meta_cells)
}
