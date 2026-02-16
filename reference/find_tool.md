# Find an external tool

Locates an external bioinformatics tool using a hybrid resolution
strategy:

1.  Explicit path if provided

2.  Conda prefix (direct path to environment) if specified

3.  Conda environment name if specified

4.  System PATH as fallback

## Usage

``` r
find_tool(
  tool,
  tool_path = NULL,
  conda_prefix = NULL,
  conda_env = NULL,
  error = TRUE
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

- error:

  Logical. Whether to error if tool not found (default: TRUE).

## Value

Character path to the tool, or NULL if not found and error = FALSE.

## Examples

``` r
if (FALSE) { # \dontrun{
find_tool("diamond")
find_tool("diamond", conda_prefix = "./this_project_env")
find_tool("diamond", conda_env = "bioinf")
find_tool("diamond", tool_path = "/opt/diamond/diamond")
} # }
```
