#' Singleton Utility Functions
#'
#' Functions for working with singleton proteins (proteins not assigned to any
#' orthogroup during clustering).
#'
#' @name singleton-functions
#' @keywords internal
NULL

#' Get singletons from clustering result
#'
#' Extracts the singletons tibble from an orthogroup_result object.
#'
#' @param orthogroups An `orthogroup_result` object.
#'
#' @return A tibble with columns `protein_id` and `assembly`.
#' @export
#'
#' @examples
#' \dontrun{
#' singletons <- get_singletons(clusters)
#' }
get_singletons <- function(orthogroups) {
  # Validate input
  if (!inherits(orthogroups, "orthogroup_result")) {
    cli::cli_abort("{.arg orthogroups} must be an {.cls orthogroup_result} object")
  }

  # Return singletons tibble, or empty tibble if none

  if (is.null(orthogroups$singletons)) {
    return(tibble::tibble(
      protein_id = character(),
      assembly = character()
    ))
  }

  orthogroups$singletons
}

#' Count singletons
#'
#' Returns the number of singleton proteins in a clustering result.
#'
#' @param orthogroups An `orthogroup_result` object.
#'
#' @return Integer count of singletons.
#' @export
#'
#' @examples
#' \dontrun{
#' n <- n_singletons(clusters)
#' }
n_singletons <- function(orthogroups) {
  # Validate input
  if (!inherits(orthogroups, "orthogroup_result")) {
    cli::cli_abort("{.arg orthogroups} must be an {.cls orthogroup_result} object")
  }

  if (is.null(orthogroups$singletons)) {
    return(0L)
  }

  nrow(orthogroups$singletons)
}

#' Summarize singletons by assembly
#'
#' Returns a summary of singleton counts per assembly, optionally including
#' the percentage of total proteins that are singletons.
#'
#' @param orthogroups An `orthogroup_result` object.
#' @param proteins Optional `protein_collection` to calculate percentages.
#'
#' @return A tibble with columns:
#'   - `assembly`: Assembly name
#'   - `n_singletons`: Number of singletons in assembly
#'   - `n_total`: Total proteins in assembly (if proteins provided)
#'   - `pct_singleton`: Percentage of proteins that are singletons (if proteins provided)
#' @export
#'
#' @examples
#' \dontrun{
#' summary <- singletons_by_assembly(clusters)
#' summary_with_pct <- singletons_by_assembly(clusters, proteins = proteins)
#' }
singletons_by_assembly <- function(orthogroups, proteins = NULL) {
  # Validate input
  if (!inherits(orthogroups, "orthogroup_result")) {
    cli::cli_abort("{.arg orthogroups} must be an {.cls orthogroup_result} object")
  }

  singletons <- get_singletons(orthogroups)

  # Handle empty singletons
  if (nrow(singletons) == 0) {
    if (is.null(proteins)) {
      return(tibble::tibble(
        assembly = character(),
        n_singletons = integer()
      ))
    } else {
      return(tibble::tibble(
        assembly = character(),
        n_singletons = integer(),
        n_total = integer(),
        pct_singleton = numeric()
      ))
    }
  }

  # Count singletons per assembly
  result <- singletons |>
    dplyr::count(.data$assembly, name = "n_singletons")

  # If proteins provided, add totals and percentages
 if (!is.null(proteins)) {
    if (!inherits(proteins, "protein_collection")) {
      cli::cli_abort("{.arg proteins} must be a {.cls protein_collection} object")
    }

    # Get total protein counts per assembly
    totals <- proteins$summary |>
      dplyr::select("assembly_name", "n_proteins") |>
      dplyr::rename(assembly = "assembly_name", n_total = "n_proteins")

    # Join and calculate percentage
    result <- result |>
      dplyr::left_join(totals, by = "assembly") |>
      dplyr::mutate(pct_singleton = .data$n_singletons / .data$n_total * 100)
  }

  result
}
