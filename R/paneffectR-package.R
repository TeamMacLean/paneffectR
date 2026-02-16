#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data .env
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  # Check for BiocManager (needed for Bioconductor deps)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    packageStartupMessage(
      "Note: BiocManager not installed. ",
      "Some features require Bioconductor packages.\n",
      "Install with: install.packages('BiocManager')"
    )
  }

  # Check for Bioconductor dependencies
  bioc_pkgs <- c("ComplexHeatmap", "Biostrings")
  missing <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    packageStartupMessage(
      "Note: Missing Bioconductor packages: ", paste(missing, collapse = ", "), "\n",
      "Install with: BiocManager::install(c('", paste(missing, collapse = "', '"), "'))"
    )
  }
}
