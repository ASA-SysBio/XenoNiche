#' Label cells by proximity to target groups within specified FOVs
#'
#' For each field-of-view (FOV), this function labels cells as:
#' \itemize{
#'   \item \strong{"target"}: cells whose \code{column_of_idents} value is in \code{cell_ids} and fall in the FOV
#'   \item \strong{"surrounding"}: non-target cells within \code{radius} (coordinate units) of any target cell in that FOV
#'   \item \strong{"others"}: remaining cells within the FOV
#'   \item \strong{"not selected for analysis"}: cells not in any of the provided FOVs
#' }
#'
#' Centroids are extracted from \code{object@images[[fov]][["centroids"]]} via
#' \code{Seurat::GetTissueCoordinates()}, then de-duplicated and intersected with the cells in the Seurat object.
#'
#' @param seurat_obj A \code{Seurat} object containing spatial images with centroid coordinates.
#' @param cell_ids Character vector of **group labels** (e.g., cell types) that define the targets.
#' @param column_of_idents Character scalar; name of the metadata column that holds the group labels to match against \code{cell_ids}.
#' @param radius Numeric radius (same units as image coordinates) used to label \emph{surrounding} cells around targets.
#' @param fov_list Character vector of image/FOV names to evaluate (must match \code{names(seurat_obj@images)}).
#' @param column_name Character scalar; the metadata column name to write the resulting labels into
#'   (default \code{"proximity_label"}). If the column already exists, it will be overwritten (with a message).
#'
#' @return The input \code{Seurat} object with a new metadata column \code{column_name}
#' containing a factor with levels \code{c("target","surrounding","others","not selected for analysis")}.
#'
#' @examples
#' \donttest{
#' # seurat_obj <- label_cells_by_proximity(
#' #   seurat_obj,
#' #   cell_ids = c("Fibroblast","Perivascular"),
#' #   column_of_idents = "subset",
#' #   radius = 50,
#' #   fov_list = c("img1","img2"),
#' #   column_name = "prox_labels"
#' # )
#' # table(seurat_obj[["prox_labels"]][,1])
#' }
#'
#' @export
#' @importFrom Seurat Cells GetTissueCoordinates
label_cells_by_proximity <- function(seurat_obj,
                                     cell_ids,
                                     column_of_idents,
                                     radius,
                                     fov_list,
                                     column_name = "proximity_label") {
  
  # --- validations -------------------------------------------------------------
  if (!is.character(column_of_idents) || length(column_of_idents) != 1) {
    stop("`column_of_idents` must be a single character string naming a metadata column.")
  }
  if (!column_of_idents %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("`column_of_idents` ('%s') not found in seurat_obj@meta.data.", column_of_idents))
  }
  if (!is.character(column_name) || length(column_name) != 1) {
    stop("`column_name` must be a single character string.")
  }
  if (!is.numeric(radius) || length(radius) != 1 || is.na(radius) || radius <= 0) {
    stop("`radius` must be a single positive numeric value.")
  }
  if (length(fov_list) == 0) {
    warning("[label_cells_by_proximity] `fov_list` is empty; all cells will be 'not selected for analysis'.")
  }
  if (!all(fov_list %in% names(seurat_obj@images))) {
    missing_fov <- setdiff(fov_list, names(seurat_obj@images))
    stop(sprintf("FOV(s) not found in object@images: %s", paste(missing_fov, collapse = ", ")))
  }
  
  # --- helper: extract centroid table (x,y; rownames = cell ids) from all images
  get_centroids <- function(object, quiet = FALSE) {
    imgs <- names(object@images)
    if (length(imgs) == 0) stop("[get_centroids] object@images is empty.")
    
    centroids_list <- lapply(imgs, function(i) {
      coords <- tryCatch(GetTissueCoordinates(object[[i]][["centroids"]]),
                         error = function(e) NULL)
      if (is.null(coords)) return(NULL)
      df <- as.data.frame(coords)
      if (!"cell" %in% names(df)) df$cell <- rownames(df)
      if (!all(c("x","y","cell") %in% names(df))) {
        if (!quiet) message(sprintf("[get_centroids] Skip '%s': missing x/y/cell.", i))
        return(NULL)
      }
      df[, c("cell","x","y")]
    })
    
    centroids_df <- do.call(rbind, Filter(Negate(is.null), centroids_list))
    if (is.null(centroids_df) || nrow(centroids_df) == 0) {
      stop("[get_centroids] Could not extract centroids from any image.")
    }
    
    # de-duplicate, keep only cells present in the object
    centroids_df <- centroids_df[!duplicated(centroids_df$cell), ]
    rownames(centroids_df) <- centroids_df$cell
    centroids_df$cell <- NULL
    
    centroids_df[intersect(rownames(centroids_df), colnames(object)), , drop = FALSE]
  }
  
  # --- compute centroids once --------------------------------------------------
  centroids_df <- get_centroids(seurat_obj, quiet = TRUE)
  
  # --- initialize global vector with "not selected for analysis" ---------------
  all_cells <- colnames(seurat_obj)
  proximity_labels <- rep("not selected for analysis", length(all_cells))
  names(proximity_labels) <- all_cells
  
  # --- select target cells by metadata column ----------------------------------
  # targets = cells whose meta[column_of_idents] is in cell_ids
  meta <- seurat_obj@meta.data
  target_cells <- rownames(meta)[meta[[column_of_idents]] %in% cell_ids]
  
  if (length(target_cells) == 0) {
    message(sprintf(
      "[label_cells_by_proximity] No targets found where %s %%in%% cell_ids; all FOV cells will be 'others'.",
      column_of_idents
    ))
  }
  
  # --- iterate per FOV, label target/surrounding/others ------------------------
  for (fov in fov_list) {
    cells_in_fov <- Cells(seurat_obj[[fov]])
    if (length(cells_in_fov) == 0) next
    
    # Everyone in the FOV starts as "others"
    fov_labels <- rep("others", length(cells_in_fov))
    names(fov_labels) <- cells_in_fov
    
    # Targets that are in this FOV
    fov_targets <- intersect(target_cells, cells_in_fov)
    if (length(fov_targets) > 0) {
      fov_labels[fov_targets] <- "target"
      
      # Coordinates for FOV cells and for FOV targets
      fov_coords    <- centroids_df[cells_in_fov, , drop = FALSE]
      target_coords <- centroids_df[fov_targets, , drop = FALSE]
      
      # Vectorized min distance to any target (use squared distances; compare to radius^2)
      dx <- outer(fov_coords$x, target_coords$x, "-")
      dy <- outer(fov_coords$y, target_coords$y, "-")
      min_sqdist <- apply(dx^2 + dy^2, 1, min)
      
      is_surrounding <- (min_sqdist <= radius^2) & !(names(fov_labels) %in% fov_targets)
      fov_labels[is_surrounding] <- "surrounding"
    }
    
    # Write back labels for this FOV into the global vector
    proximity_labels[names(fov_labels)] <- fov_labels
  }
  
  # Factor with stable level order (helps consistent legends/plots)
  proximity_labels <- factor(
    proximity_labels,
    levels = c("target","surrounding","others","not selected for analysis")
  )
  
  # --- store in metadata under column_name -------------------------------------
  if (column_name %in% colnames(seurat_obj@meta.data)) {
    message(sprintf("Overwriting existing metadata column '%s'.", column_name))
  }
  seurat_obj[[column_name]] <- proximity_labels
  
  # Return the modified object so users can continue piping
  return(seurat_obj)
}
