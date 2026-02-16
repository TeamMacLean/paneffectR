# Check if an external tool is installed

Convenience wrapper around find_tool() that returns TRUE/FALSE.

## Usage

``` r
check_tool_installed(
  tool,
  tool_path = NULL,
  conda_prefix = NULL,
  conda_env = NULL
)
```

## Arguments

- tool:

  Character. Name of the tool (e.g., "diamond", "mmseqs").

- tool_path:

  Character. Optional explicit path to the tool binary.

- conda_prefix:

  Character. Optional direct path to a conda/mamba environment (e.g.,
  "./this_project_env"). The tool is expected at `<prefix>/bin/<tool>`.

- conda_env:

  Character. Optional conda environment name (looked up in standard
  conda envs directory).

## Value

Logical. TRUE if tool is found, FALSE otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
if (check_tool_installed("diamond")) {
  # use diamond
}
if (check_tool_installed("diamond", conda_prefix = "./this_project_env")) {
  # use diamond from prefix env
}
} # }
```
