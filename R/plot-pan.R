#' Pan-Genome Structure Visualization
#'
#' Functions for visualizing the structure of the pan-genome (or pan-effectorome),
#' showing the distribution of proteins across core, accessory, rare, unique,
#' and singleton categories.
#'
#' @name plot-pan
#' @keywords internal
NULL

#' Categorize orthogroups by presence pattern
#'
#' Internal helper that assigns each orthogroup to a category based on
#' how many assemblies contain it.
#'
#' @param clusters An `orthogroup_result` object.
#' @param pa Optional pre-built `pa_matrix`. If NULL, built internally.
#' @param rare_threshold Proportion threshold for "Rare" category (default 0.10).
#'   Orthogroups present in <= this proportion of assemblies (but > 1) are "Rare".
#' @param accessory_threshold Proportion threshold for "Accessory" category (default 0.60).
#'   Orthogroups present in <= this proportion (but > rare_threshold) are "Accessory".
#'
#' @return A tibble with columns: orthogroup_id, n_present, n_assemblies, category
#' @keywords internal
categorize_orthogroups <- function(clusters,
                                   pa = NULL,
                                   rare_threshold = 0.10,
                                   accessory_threshold = 0.60) {
  # Build PA matrix if not provided

if (is.null(pa)) {
    pa <- build_pa_matrix(clusters, type = "binary")
  }

  n_assemblies <- ncol(pa$matrix)
  presence_count <- rowSums(pa$matrix)

 # Calculate thresholds as counts
  rare_count <- max(1, floor(n_assemblies * rare_threshold))
  accessory_count <- floor(n_assemblies * accessory_threshold)

  # Categorize
  tibble::tibble(
    orthogroup_id = names(presence_count),
    n_present = as.integer(presence_count),
    n_assemblies = n_assemblies,
    category = dplyr::case_when(
      n_present == n_assemblies ~ "Core",
      n_present == 1 ~ "Unique",
      n_present <= rare_count ~ "Rare",
      n_present <= accessory_count ~ "Accessory",
      TRUE ~ "Accessory"
    )
  )
}

#' Plot pan-genome structure
#'
#' Creates a bar chart showing the composition of the pan-genome by category:
#' Core (present in all assemblies), Accessory, Rare, Unique (present in one
#' assembly), and optionally Singletons (proteins not in any orthogroup).
#'
#' @param clusters An `orthogroup_result` object.
#' @param pa Optional pre-built `pa_matrix`. If NULL, built internally.
#' @param rare_threshold Proportion threshold for "Rare" category (default 0.10).
#'   Orthogroups present in <= this proportion of assemblies (but > 1) are "Rare".
#' @param accessory_threshold Proportion threshold for "Accessory" category (default 0.60).
#'   Orthogroups present in > rare_threshold and <= this proportion are "Accessory".
#' @param exclude_singletons Logical. If TRUE, singletons are not shown (default FALSE).
#' @param colors Named vector of colors for each category. If NULL, uses default palette.
#' @param show_counts Logical. If TRUE, display count labels on bars (default TRUE).
#' @param title Plot title. If NULL, uses default.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_pan_structure(clusters)
#'
#' # Custom thresholds
#' plot_pan_structure(clusters, rare_threshold = 0.15, accessory_threshold = 0.50)
#'
#' # Without singletons
#' plot_pan_structure(clusters, exclude_singletons = TRUE)
#' }
plot_pan_structure <- function(clusters,
                               pa = NULL,
                               rare_threshold = 0.10,
                               accessory_threshold = 0.60,
                               exclude_singletons = FALSE,
                               colors = NULL,
                               show_counts = TRUE,
                               title = NULL) {
  # Validate input
  if (!inherits(clusters, "orthogroup_result")) {
    cli::cli_abort("{.arg clusters} must be an {.cls orthogroup_result} object")
 }

  # Categorize orthogroups
  categories <- categorize_orthogroups(
    clusters, pa,
    rare_threshold = rare_threshold,
    accessory_threshold = accessory_threshold
  )

  # Count by category
  category_counts <- categories |>
    dplyr::count(.data$category, name = "n")

  # Add singletons if not excluded
  if (!exclude_singletons) {
    n_sing <- n_singletons(clusters)
    category_counts <- dplyr::bind_rows(
      category_counts,
      tibble::tibble(category = "Singleton", n = n_sing)
    )
  }

  # Set category order
  cat_levels <- c("Core", "Accessory", "Rare", "Unique")
  if (!exclude_singletons) {
    cat_levels <- c(cat_levels, "Singleton")
  }
  category_counts <- category_counts |>
    dplyr::mutate(category = factor(.data$category, levels = cat_levels)) |>
    dplyr::filter(!is.na(.data$category))

  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "Core" = "#2166AC",
      "Accessory" = "#67A9CF",
      "Rare" = "#D1E5F0",
      "Unique" = "#FDDBC7",
      "Singleton" = "#B2182B"
    )
  }

  # Default title
  if (is.null(title)) {
    title <- "Pan-Genome Structure"
  }

  # Build plot
  p <- ggplot2::ggplot(category_counts, ggplot2::aes(x = .data$category, y = .data$n, fill = .data$category)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = title,
      x = "Category",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Add count labels
  if (show_counts) {
    max_n <- max(category_counts$n, na.rm = TRUE)
    p <- p +
      ggplot2::geom_text(ggplot2::aes(label = .data$n), vjust = -0.5) +
      ggplot2::ylim(0, max_n * 1.15)
  }

  p
}

#' Plot assembly composition
#'
#' Creates a stacked bar chart showing the composition of each assembly's
#' proteins by category: Core, Accessory, Rare, Unique, and optionally Singletons.
#'
#' @inheritParams plot_pan_structure
#' @param position Bar position: "stack" for counts (default), "fill" for proportions.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_assembly_composition(clusters)
#'
#' # Show proportions instead of counts
#' plot_assembly_composition(clusters, position = "fill")
#'
#' # Without singletons
#' plot_assembly_composition(clusters, exclude_singletons = TRUE)
#' }
plot_assembly_composition <- function(clusters,
                                      pa = NULL,
                                      rare_threshold = 0.10,
                                      accessory_threshold = 0.60,
                                      exclude_singletons = FALSE,
                                      colors = NULL,
                                      position = "stack",
                                      title = NULL) {
  # Validate input
  if (!inherits(clusters, "orthogroup_result")) {
    cli::cli_abort("{.arg clusters} must be an {.cls orthogroup_result} object")
  }

  # Categorize orthogroups
  categories <- categorize_orthogroups(
    clusters, pa,
    rare_threshold = rare_threshold,
    accessory_threshold = accessory_threshold
  )

  # Count proteins per category per assembly
  og_by_assembly <- clusters$orthogroups |>
    dplyr::left_join(categories, by = "orthogroup_id") |>
    dplyr::count(.data$assembly, .data$category)

  # Add singletons if not excluded
  if (!exclude_singletons) {
    singletons <- get_singletons(clusters)
    if (nrow(singletons) > 0) {
      singleton_by_asm <- singletons |>
        dplyr::count(.data$assembly) |>
        dplyr::mutate(category = "Singleton")
      og_by_assembly <- dplyr::bind_rows(og_by_assembly, singleton_by_asm)
    }
  }

  # Set category order
  cat_levels <- c("Core", "Accessory", "Rare", "Unique")
  if (!exclude_singletons) {
    cat_levels <- c(cat_levels, "Singleton")
  }
  og_by_assembly <- og_by_assembly |>
    dplyr::mutate(category = factor(.data$category, levels = cat_levels)) |>
    dplyr::filter(!is.na(.data$category))

  # Default colors
  if (is.null(colors)) {
    colors <- c(
      "Core" = "#2166AC",
      "Accessory" = "#67A9CF",
      "Rare" = "#D1E5F0",
      "Unique" = "#FDDBC7",
      "Singleton" = "#B2182B"
    )
  }

  # Default title
  if (is.null(title)) {
    title <- "Assembly Composition"
  }

  # Y-axis label depends on position
  y_label <- if (position == "fill") "Proportion" else "Number of Proteins"

  # Build plot
  ggplot2::ggplot(og_by_assembly, ggplot2::aes(x = .data$assembly, y = .data$n, fill = .data$category)) +
    ggplot2::geom_col(position = position) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = title,
      x = "Assembly",
      y = y_label,
      fill = "Category"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
