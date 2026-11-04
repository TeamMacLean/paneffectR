# Derive an assembly name from a FASTA file path

Strips the directory, the extension, and any `.gz`/`.bz2`/`.xz` suffix.

## Usage

``` r
fasta_assembly_name(path)
```

## Arguments

- path:

  Character path to a FASTA file.

## Value

Character assembly name.
