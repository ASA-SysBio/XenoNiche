#' Label cells by proximity to target cells within specified FOVs
#'
#' For each field-of-view (FOV), this function labels cells as:
#' \itemize{
#'   \item \strong{"target"}: cells explicitly provided by \code{cell_ids} that fall in the FOV
#'   \item \strong{"surrounding"}: non-target cells within \code{radius} (in pixel/coordinate units) of any target cell in that FOV
#'   \item \strong{"others"}: remaining cells within the FOV
#'   \item \strong{"not selected for analysis"}: cells not in any of the provided FOVs
#' }
#'
#' Centroids are extracted from \code{object@images[[fov]][["centroids"]]} via
#' \code{Seurat::GetTissueCoordinates()}, then de-duplicated and intersected with the cells in the Seurat object.
#'
#' @param seurat_obj A \code{Seurat} object containing spatial images with centroid coordinates.
#' @param cell_ids Character vector of target cell barcodes/IDs to define the proximity query.
#' @param radius Numeric radius (same units as image coordinates) used to label \emph{surrounding} cells around targets.
#' @param fov_list Character vector of image/FOV names to evaluate (must match \code{names(seurat_obj@images)}).
#' @param column_name (Ignored in return; kept for API parity.) Intended metadata column name if user later stores results.
#'
#' @return A factor named by cell IDs, with levels:
#'   \code{c("target","surrounding","others","not selected for analysis")}.
#'
#' @examples
#' \donttest{
#' # labels <- label_cells_by_proximity(seurat_obj, cell_ids = c("cellA","cellB"),
#' #                                    radius = 50, fov_list = c("img1","img2"))
#' }
#'
#' @export
#' @importFrom Seurat Cells GetTissueCoordinates
label_cells_by_proximity <- function(seurat_obj,
                                     cell_ids,
                                     radius,
                                     fov_list,
                                     column_name = "proximity_label") {

  # --- helper: extract centroid table (x,y, rownames = cell ids) from all images ---
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

  # --- compute centroids table once ---
  centroids_df <- get_centroids(seurat_obj, quiet = TRUE)

  # --- initialize full vector with "not selected for analysis" ---
  all_cells <- colnames(seurat_obj)
  proximity_labels <- rep("not selected for analysis", length(all_cells))
  names(proximity_labels) <- all_cells

  # --- iterate per FOV, label target/surrounding/others ---
  if (length(fov_list) == 0) {
    warning("[label_cells_by_proximity] fov_list is empty; returning all 'not selected for analysis'.")
    return(factor(proximity_labels,
                  levels = c("target","surrounding","others","not selected for analysis")))
  }

  if (!all(fov_list %in% names(seurat_obj@images))) {
    missing_fov <- setdiff(fov_list, names(seurat_obj@images))
    stop(sprintf("FOV(s) not found in object@images: %s", paste(missing_fov, collapse = ", ")))
  }

  # Basic validation
  if (!is.numeric(radius) || length(radius) != 1 || is.na(radius) || radius <= 0) {
    stop("`radius` must be a single positive numeric value.")
  }
  if (length(cell_ids) == 0) {
    # still proceed; all cells in FOVs become "others"
    message("[label_cells_by_proximity] `cell_ids` is empty; no 'target' or 'surrounding' labels will be assigned.")
  }

  for (fov in fov_list) {
    # Cells in this FOV/image
    cells_in_fov <- Cells(seurat_obj[[fov]])
    if (length(cells_in_fov) == 0) next

    # Start everyone in the FOV as "others"
    fov_labels <- rep("others", length(cells_in_fov))
    names(fov_labels) <- cells_in_fov

    # Targets in this FOV
    fov_targets <- intersect(cell_ids, cells_in_fov)
    if (length(fov_targets) > 0) {
      fov_labels[fov_targets] <- "target"

      # coordinates for just the FOV cells
      fov_coords <- centroids_df[cells_in_fov, , drop = FALSE]

      # Vectorized distance labeling to avoid re-labelling targets
      target_coords <- centroids_df[fov_targets, , drop = FALSE]
      # Compute pairwise squared distances (broadcast) and take min over targets
      # (small optimization: use squared distances to avoid many sqrt calls)
      dx <- outer(fov_coords$x, target_coords$x, "-")
      dy <- outer(fov_coords$y, target_coords$y, "-")
      min_sqdist <- apply(dx^2 + dy^2, 1, min)

      # surrounding if within radius and not a target itself
      is_surrounding <- (min_sqdist <= radius^2) & !(names(fov_labels) %in% fov_targets)
      fov_labels[is_surrounding] <- "surrounding"
    }

    # Write back labels for this FOV into the global vector
    proximity_labels[names(fov_labels)] <- fov_labels
  }

  factor(proximity_labels,
         levels = c("target","surrounding","others","not selected for analysis"))
}

