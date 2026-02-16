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
