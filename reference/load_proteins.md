# Load proteins from multiple assemblies

Main entry point for loading protein data. Discovers and loads FASTA
files, optionally with associated score and name mapping files.

## Usage

``` r
load_proteins(
  fasta_dir,
  score_dir = NULL,
  pattern = "*.faa",
  name_mapping = TRUE
)
```

## Arguments

- fasta_dir:

  Character. Directory containing FASTA files.

- score_dir:

  Character. Optional directory containing score CSV files.

- pattern:

  Character. Glob pattern for FASTA files (default: "\*.faa").

- name_mapping:

  Logical. Whether to load name mapping files.

## Value

A `protein_collection` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load from omnieff output
proteins <- load_proteins(
  fasta_dir = "reformatted/",
  score_dir = "scored/"
)

# Load just FASTAs
proteins <- load_proteins(fasta_dir = "my_fastas/")
} # }
```
