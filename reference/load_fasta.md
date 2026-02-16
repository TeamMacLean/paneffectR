# Load a single FASTA file

Parses a FASTA file into a tibble of protein IDs and sequences.

## Usage

``` r
load_fasta(path)
```

## Arguments

- path:

  Character. Path to FASTA file.

## Value

A tibble with columns: protein_id, sequence.
