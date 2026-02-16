# Compute reciprocal best hits from BLAST results

Given a table of BLAST hits with bitscores, identifies pairs where each
protein is the other's best hit.

## Usage

``` r
compute_rbh(hits)
```

## Arguments

- hits:

  A tibble with columns: qseqid, sseqid, bitscore

## Value

A tibble with columns: protein_a, protein_b (RBH pairs)
