# Create an orthogroup_result object

Result of clustering proteins across assemblies.

## Usage

``` r
new_orthogroup_result(
  orthogroups,
  method,
  parameters = list(),
  singletons = NULL,
  stats = NULL
)
```

## Arguments

- orthogroups:

  A tibble with columns: orthogroup_id, assembly, protein_id.

- method:

  Character. Clustering method used.

- parameters:

  List. Clustering parameters used.

- singletons:

  A tibble. Proteins not assigned to any group.

- stats:

  A tibble. Orthogroup sizes and assembly coverage.

## Value

An `orthogroup_result` S3 object.
