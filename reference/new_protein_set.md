# Create a protein_set object

Represents proteins from a single assembly.

## Usage

``` r
new_protein_set(
  assembly_name,
  proteins,
  metadata = list(),
  source_files = list()
)
```

## Arguments

- assembly_name:

  Character. Unique identifier for the assembly.

- proteins:

  A tibble with columns: protein_id, sequence, and optionally score,
  rank, and other metadata.

- metadata:

  List. Optional assembly-level metadata.

- source_files:

  List. Paths to original .faa and .csv files.

## Value

A `protein_set` S3 object.

## Examples

``` r
ps <- new_protein_set(
  assembly_name = "isolate_A",
  proteins = tibble::tibble(
    protein_id = c("A_001", "A_002"),
    sequence = c("MKVLAAGIV", "MSTNPKPQR")
  )
)
ps
#> -- protein_set: isolate_A --
#> 2 proteins
```
