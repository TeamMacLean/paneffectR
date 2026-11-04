# Run DIAMOND reciprocal best hits clustering

Performs all-vs-all DIAMOND blastp search and identifies reciprocal best
hits to build orthogroups.

## Usage

``` r
run_diamond_rbh(
  proteins,
  min_identity = 30,
  min_coverage = 50,
  evalue = 1e-05,
  threads = NULL,
  mode = "fast",
  tool_path = NULL,
  conda_prefix = NULL,
  conda_env = NULL,
  keep_temp = FALSE,
  timeout = 3600,
  quiet = FALSE
)
```

## Arguments

- proteins:

  A protein_collection object

- min_identity:

  Minimum percent identity (default 30)

- min_coverage:

  Minimum query coverage (default 50)

- evalue:

  E-value threshold (default 1e-5)

- threads:

  Number of threads (default: auto-detect)

- mode:

  "fast" or "thorough" (default "fast")

- tool_path:

  Optional explicit path to diamond binary

- conda_prefix:

  Optional path to a conda/mamba environment

- conda_env:

  Optional name of a conda/mamba environment

- keep_temp:

  Keep temporary files for debugging (default FALSE)

- timeout:

  Seconds to allow each DIAMOND step before aborting

- quiet:

  Suppress progress messages (default FALSE)

## Value

An orthogroup_result object
