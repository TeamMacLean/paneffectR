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
  keep_temp = FALSE
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

- keep_temp:

  Keep temporary files for debugging (default FALSE)

## Value

An orthogroup_result object
