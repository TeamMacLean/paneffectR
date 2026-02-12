#' UpSet Plot Visualization
#'
#' @name plot-upset
#' @keywords internal
NULL

#' Plot UpSet diagram of orthogroup sharing
#'
#' Creates an UpSet plot showing intersections of orthogroups across assemblies.
#' Useful for visualizing core, accessory, and unique orthogroups.
#'
#' @param pa A `pa_matrix` object.
#' @param min_size Integer. Minimum intersection size to show (default: 1).
#' @param max_sets Integer. Maximum number of sets to display. If NULL (default),
#'   all assemblies are shown.
#' @param order_by Character. How to order intersections: "freq" (by size,
#'   default) or "degree" (by number of sets in intersection).
#' @param ... Additional arguments passed to UpSetR::upset().
#'
#' @return Draws an UpSet plot (base graphics). Returns invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_upset(pa)
#' plot_upset(pa, min_size = 2)
#' plot_upset(pa, order_by = "degree")
#' }
plot_upset <- function(pa,
                       min_size = 1,
                       max_sets = NULL,
                       order_by = "freq",
                       ...) {
  # Validate input
  if (!inherits(pa, "pa_matrix")) {
    cli::cli_abort("{.arg pa} must be a {.cls pa_matrix} object")
  }

  # Validate order_by parameter
  valid_order <- c("freq", "degree")
  if (!order_by %in% valid_order) {
    cli::cli_abort(
      "{.arg order_by} must be one of {.val {valid_order}}, not {.val {order_by}}"
    )
  }

  # Check for empty matrix
  if (nrow(pa$matrix) == 0) {
    cli::cli_abort("Cannot create UpSet plot: no orthogroups in the matrix")
  }

  # Check for single assembly
  if (ncol(pa$matrix) < 2) {
    cli::cli_abort(
      "UpSet plots require at least 2 assemblies, but only {ncol(pa$matrix)} found"
    )
  }

  # Convert to binary if needed (UpSetR expects 0/1 integers)
  mat <- pa$matrix
  mat[is.na(mat)] <- 0
  mat[mat > 0] <- 1
  storage.mode(mat) <- "integer"

  # UpSetR expects data frame with columns as sets (assemblies)
  # Our matrix has rows=orthogroups, cols=assemblies, which is correct
  upset_df <- as.data.frame(mat)

  # Apply max_sets filter if specified
  if (!is.null(max_sets) && max_sets < ncol(upset_df)) {
    # Select assemblies with most orthogroups
    set_sizes <- colSums(upset_df)
    top_sets <- names(sort(set_sizes, decreasing = TRUE))[seq_len(max_sets)]
    upset_df <- upset_df[, top_sets, drop = FALSE]
  }

  # Filter by min_size will be handled by UpSetR parameter

  # Create UpSet plot
  UpSetR::upset(
    upset_df,
    nsets = ncol(upset_df),
    nintersects = NA,  # Show all intersections
    order.by = order_by,
    decreasing = TRUE,
    cutoff = min_size,
    mb.ratio = c(0.6, 0.4),
    sets.bar.color = "#56B4E9",
    main.bar.color = "#009E73",
    matrix.color = "#009E73",
    keep.order = TRUE,
    set_size.show = TRUE,
    text.scale = 1.2,
    ...
  )

  invisible(NULL)
}
