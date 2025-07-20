#' Compute minimum distances between query and reference cell types
#'
#' @param seurat_obj A Seurat object containing spatial coordinates in metadata.
#' @param x_col Name of column for x-coordinates in metadata.
#' @param y_col Name of column for y-coordinates in metadata.
#' @param group_col Column indicating cell type or grouping.
#' @param query_groups Character vector of cell groups to use as queries.
#' @param reference_groups Character vector of cell groups to use as references.
#'
#' @return A data frame with nearest reference cell info for each query cell.
#' @export
compute_min_celltype_distances <- function(seurat_obj,
                                           x_col = "x",
                                           y_col = "y",
                                           group_col = "subset",
                                           query_groups,
                                           reference_groups) {
  meta <- seurat_obj@meta.data
  meta$cell_id <- rownames(meta)

  if (!all(c(x_col, y_col) %in% colnames(meta))) {
    stop("Specified x_col or y_col does not exist in metadata.")
  }

  query_cells <- meta[meta[[group_col]] %in% query_groups, ]
  ref_cells   <- meta[meta[[group_col]] %in% reference_groups, ]

  if (nrow(query_cells) == 0) stop("No cells found for query_groups.")
  if (nrow(ref_cells) == 0) stop("No cells found for reference_groups.")

  query_coords <- as.matrix(query_cells[, c(x_col, y_col)])
  ref_coords   <- as.matrix(ref_cells[, c(x_col, y_col)])

  result <- data.frame()
  query_label <- paste(sort(unique(query_cells[[group_col]])), collapse = ",")
  reference_label <- paste(sort(unique(ref_cells[[group_col]])), collapse = ",")

  for (i in seq_len(nrow(query_coords))) {
    q <- query_coords[i, ]
    dists <- sqrt((ref_coords[,1] - q[1])^2 + (ref_coords[,2] - q[2])^2)
    min_idx <- which.min(dists)

    result <- rbind(result, data.frame(
      query_cell_id     = query_cells$cell_id[i],
      query_group       = query_label,
      query_x           = q[1],
      query_y           = q[2],
      reference_cell_id = ref_cells$cell_id[min_idx],
      reference_group   = reference_label,
      reference_x       = ref_coords[min_idx, 1],
      reference_y       = ref_coords[min_idx, 2],
      min_distance      = dists[min_idx]
    ))
  }

  return(result)
}
