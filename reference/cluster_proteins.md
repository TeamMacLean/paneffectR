# Cluster proteins across assemblies

Groups proteins from a protein_collection into orthogroups using
sequence similarity. Multiple clustering methods are available.

## Usage

``` r
cluster_proteins(
  proteins,
  method = "diamond_rbh",
  mode = "fast",
  min_identity = 30,
  min_coverage = 50,
  evalue = 1e-05,
  threads = NULL,
  tool_path = NULL,
  conda_prefix = NULL,
  keep_temp = FALSE
)
```

## Arguments

- proteins:

  A protein_collection object

- method:

  Clustering method: "diamond_rbh" (default), "orthofinder", or
  "mmseqs2"

- mode:

  Speed/sensitivity trade-off: "fast" (default) or "thorough"

- min_identity:

  Minimum percent identity threshold (default 30)

- min_coverage:

  Minimum query coverage threshold (default 50)

- evalue:

  E-value threshold (default 1e-5)

- threads:

  Number of CPU threads (default: auto-detect)

- tool_path:

  Optional explicit path to the clustering tool binary

- conda_prefix:

  Optional direct path to a conda/mamba environment containing the tool
  (e.g., "./this_project_env"). Tool is expected at
  `<prefix>/bin/<tool>`.

- keep_temp:

  Keep temporary files for debugging (default FALSE)

## Value

An orthogroup_result object containing:

- orthogroups: tibble with orthogroup_id, assembly, protein_id

- method: the clustering method used

- parameters: list of parameters used

- singletons: proteins not assigned to any orthogroup

## Examples

``` r
if (FALSE) { # \dontrun{
# Load proteins
proteins <- load_proteins(fasta_dir = "assemblies/")

# Cluster with default settings (DIAMOND RBH)
result <- cluster_proteins(proteins)

# Use thorough mode with stricter thresholds
result <- cluster_proteins(
  proteins,
  mode = "thorough",
  min_identity = 70,
  min_coverage = 80
)

# Use a project-local conda environment
result <- cluster_proteins(
  proteins,
  conda_prefix = "./this_project_env"
)
} # }
```
