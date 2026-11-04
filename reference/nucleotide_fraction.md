# Fraction of residues that are nucleotide letters

Used to detect nucleotide FASTA files mislabelled as protein. Looks at
up to `max_seqs` sequences.

## Usage

``` r
nucleotide_fraction(sequences, max_seqs = 1000L)
```

## Arguments

- sequences:

  Character vector of sequences.

- max_seqs:

  Maximum number of sequences to inspect.

## Value

Numeric fraction in `[0, 1]`, or `NA` if there are no residues.
