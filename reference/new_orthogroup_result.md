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

## Examples

``` r
og <- tibble::tibble(
  orthogroup_id = c("OG0001", "OG0001", "OG0002", "OG0002"),
  assembly = c("isolate_A", "isolate_B", "isolate_A", "isolate_B"),
  protein_id = c("A_001", "B_001", "A_002", "B_002")
)
singles <- tibble::tibble(assembly = "isolate_A", protein_id = "A_003")
result <- new_orthogroup_result(og, method = "manual", singletons = singles)
result
#> -- orthogroup_result (manual) --
#> 2 orthogroups
#> 1 singleton
```
