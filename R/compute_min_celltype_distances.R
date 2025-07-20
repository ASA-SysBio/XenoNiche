#' Compute minimum distances between query and reference cell types
#'
#' For each query cell, computes the Euclidean distance to the nearest reference cell,
#' using spatial coordinates (e.g., x and y centroid positions). Supports combining
#' or separating query and reference groups.
#'
#' @param seurat_obj A Seurat object containing spatial coordinates in its metadata.
#' @param x_col Column name for x-coordinates (default: "x").
#' @param y_col Column name for y-coordinates (default: "y").
#' @param group_col Column in metadata defining cell group labels (e.g., cell type).
#' @param query_groups Vector of group names to be used as query cells.
#' @param reference_groups Vector of group names to be used as reference cells.
#' @param combine_query_groups Logical. Whether to combine all query groups into one (default: TRUE).
#' @param combine_reference_groups Logical. Whether to combine all reference groups into one (default: TRUE).
#'
#' @return A data frame with nearest reference cell info for each query cell.
#' @export
compute_min_celltype_distances <- function(seurat_obj,
                                           x_col = "x",
                                           y_col = "y",
                                           group_col = "subset",
                                           query_groups,
                                           reference_groups,
                                           combine_query_groups = TRUE,
                                           combine_reference_groups = TRUE) {
  # Extract metadata from the Seurat object
  meta <- seurat_obj@meta.data
  meta$cell_id <- rownames(meta)  # Add cell IDs as a new column

  # Check that the specified coordinate columns exist
  if (!all(c(x_col, y_col) %in% colnames(meta))) {
    stop("Specified x_col or y_col does not exist in metadata.")
  }

  # Initialize an empty dataframe to store results
  result <- data.frame()

  # Case 1: Combine all query and reference groups into single groups
  if (combine_query_groups && combine_reference_groups) {
    # Subset cells
    query_cells <- meta[meta[[group_col]] %in% query_groups, ]
    ref_cells   <- meta[meta[[group_col]] %in% reference_groups, ]

    if (nrow(query_cells) == 0) stop("No cells found for query_groups.")
    if (nrow(ref_cells) == 0) stop("No cells found for reference_groups.")

    # Convert coordinates to matrix format
    query_coords <- as.matrix(query_cells[, c(x_col, y_col)])
    ref_coords   <- as.matrix(ref_cells[, c(x_col, y_col)])

    # Create group labels by concatenating group names
    query_label <- paste(sort(unique(query_cells[[group_col]])), collapse = ",")
    reference_label <- paste(sort(unique(ref_cells[[group_col]])), collapse = ",")

    # For each query cell, compute the distance to all reference cells and select the minimum
    for (i in seq_len(nrow(query_coords))) {
      q <- query_coords[i, ]
      dists <- sqrt((ref_coords[,1] - q[1])^2 + (ref_coords[,2] - q[2])^2)
      min_idx <- which.min(dists)

      # Store result
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

  } else {
    # Case 2: Separate groups (or partially combined depending on the flags)

    # Prepare lists of query and reference groups based on whether they should be combined
    query_list <- if (combine_query_groups) list(query_groups) else as.list(query_groups)
    reference_list <- if (combine_reference_groups) list(reference_groups) else as.list(reference_groups)

    # Loop through each query group (or combined group)
    for (q_group in query_list) {
      q_cells <- meta[meta[[group_col]] %in% q_group, ]
      if (nrow(q_cells) == 0) next
      q_coords <- as.matrix(q_cells[, c(x_col, y_col)])
      q_label <- paste(sort(unique(q_cells[[group_col]])), collapse = ",")

      # Loop through each reference group (or combined group)
      for (r_group in reference_list) {
        r_cells <- meta[meta[[group_col]] %in% r_group, ]
        if (nrow(r_cells) == 0) next
        r_coords <- as.matrix(r_cells[, c(x_col, y_col)])
        r_label <- paste(sort(unique(r_cells[[group_col]])), collapse = ",")

        # For each query cell in this group, compute distances to current reference group
        for (i in seq_len(nrow(q_coords))) {
          q <- q_coords[i, ]
          dists <- sqrt((r_coords[,1] - q[1])^2 + (r_coords[,2] - q[2])^2)
          min_idx <- which.min(dists)

          # Store result
          result <- rbind(result, data.frame(
            query_cell_id     = q_cells$cell_id[i],
            query_group       = q_label,
            query_x           = q[1],
            query_y           = q[2],
            reference_cell_id = r_cells$cell_id[min_idx],
            reference_group   = r_label,
            reference_x       = r_coords[min_idx, 1],
            reference_y       = r_coords[min_idx, 2],
            min_distance      = dists[min_idx]
          ))
        }
      }
    }
  }

  return(result)
}
