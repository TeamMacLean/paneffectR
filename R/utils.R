#' Utility Functions
#'
#' @name utils
#' @keywords internal
NULL

#' Find an external tool
#'
#' Locates an external bioinformatics tool using a hybrid resolution strategy:
#' 1. Explicit path if provided
#' 2. Conda environment if specified
#' 3. System PATH as fallback
#'
#' @param tool Character. Name of the tool (e.g., "diamond", "mmseqs").
#' @param tool_path Character. Optional explicit path to the tool.
#' @param conda_env Character. Optional conda environment name.
#' @param error Logical. Whether to error if tool not found (default: TRUE).
#'
#' @return Character path to the tool, or NULL if not found and error = FALSE.
#' @export
#'
#' @examples
#' \dontrun
#' find_tool("diamond")
#' find_tool("diamond", conda_env = "bioinf")
#' find_tool("diamond", tool_path = "/opt/diamond/diamond")
#' }
find_tool <- function(tool,
                      tool_path = NULL,
                      conda_env = NULL,
                      error = TRUE) {
  # 1. Explicit path provided
  if (!is.null(tool_path)) {
    if (file.exists(tool_path)) {
      return(tool_path)
    }
    if (error) {
      cli::cli_abort("Tool path not found: {.path {tool_path}}")
    }
    return(NULL)
  }

  # 2. Conda environment specified
  if (!is.null(conda_env)) {
    conda_base <- find_conda_base()
    if (!is.null(conda_base)) {
      tool_in_env <- file.path(conda_base, "envs", conda_env, "bin", tool)
      if (file.exists(tool_in_env)) {
        return(tool_in_env)
      }
    }
    if (error) {
      cli::cli_abort("Tool {.val {tool}} not found in conda env: {.val {conda_env}}")
    }
    return(NULL)
  }

  # 3. Fall back to PATH
  path <- Sys.which(tool)
  if (nzchar(path)) {
    return(unname(path))
  }

  if (error) {
    cli::cli_abort(c(
      "Tool {.val {tool}} not found",
      "i" = "Install with: {.code mamba install -c bioconda {tool}}",
      "i" = "Or specify path: {.arg tool_path} or {.arg conda_env}"
    ))
  }
  NULL
}

#' Check if an external tool is installed
#'
#' Convenience wrapper around find_tool() that returns TRUE/FALSE.
#'
#' @inheritParams find_tool
#'
#' @return Logical. TRUE if tool is found, FALSE otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' if (check_tool_installed("diamond")) {
#'   # use diamond
#' }
#' }
check_tool_installed <- function(tool,
                                 tool_path = NULL,
                                 conda_env = NULL) {
  !is.null(find_tool(tool, tool_path, conda_env, error = FALSE))
}

#' Find conda base directory
#'
#' Attempts to locate the conda/mamba base directory.
#'
#' @return Character path to conda base, or NULL if not found.
#' @keywords internal
find_conda_base <- function() {
  # Try CONDA_PREFIX first (set when conda is active)
  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = NA)
  if (!is.na(conda_prefix)) {
    # CONDA_PREFIX points to active env; base is typically parent of envs
    # If we're in base, CONDA_PREFIX is the base
    if (basename(dirname(conda_prefix)) == "envs") {
      return(dirname(dirname(conda_prefix)))
    }
    # Check if this is base itself
    if (dir.exists(file.path(conda_prefix, "envs"))) {
      return(conda_prefix)
    }
  }

  # Try MAMBA_ROOT_PREFIX (micromamba)
  mamba_root <- Sys.getenv("MAMBA_ROOT_PREFIX", unset = NA)
  if (!is.na(mamba_root) && dir.exists(mamba_root)) {
    return(mamba_root)
  }

  # Try CONDA_EXE and work backwards
  conda_exe <- Sys.getenv("CONDA_EXE", unset = NA)
  if (!is.na(conda_exe) && file.exists(conda_exe)) {
    # conda_exe is typically <base>/bin/conda or <base>/condabin/conda
    base <- dirname(dirname(conda_exe))
    if (dir.exists(file.path(base, "envs"))) {
      return(base)
    }
  }

  # Common locations
  home <- Sys.getenv("HOME")
  candidates <- c(
    file.path(home, "miniforge3"),
    file.path(home, "mambaforge"),
    file.path(home, "miniconda3"),
    file.path(home, "anaconda3"),
    "/opt/conda"
  )

  for (cand in candidates) {
    if (dir.exists(cand) && dir.exists(file.path(cand, "envs"))) {
      return(cand)
    }
  }

  NULL
}

#' Get number of available CPU threads
#'
#' Returns the number of CPU threads available for parallel processing.
#'
#' @param default Integer. Default value if detection fails (default: 1).
#'
#' @return Integer number of threads.
#' @keywords internal
get_threads <- function(default = 1L) {
  # Try parallel package
  if (requireNamespace("parallel", quietly = TRUE)) {
    n <- parallel::detectCores(logical = TRUE)
    if (!is.na(n) && n > 0) {
      return(as.integer(n))
    }
  }
  as.integer(default)
}
