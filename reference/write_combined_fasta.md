# Write protein_collection to single FASTA file

Combines all proteins from all assemblies into a single FASTA file
suitable for DIAMOND database creation and searching.

## Usage

``` r
write_combined_fasta(proteins, path)
```

## Arguments

- proteins:

  A protein_collection object

- path:

  Path to write the FASTA file

## Value

Invisibly returns the path
