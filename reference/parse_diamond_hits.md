# Parse DIAMOND tabular output

Reads DIAMOND blastp output in tabular format and filters by identity,
coverage, and e-value thresholds.

## Usage

``` r
parse_diamond_hits(path, min_identity, min_coverage, evalue)
```

## Arguments

- path:

  Path to DIAMOND output file

- min_identity:

  Minimum percent identity threshold

- min_coverage:

  Minimum query coverage threshold

- evalue:

  Maximum e-value threshold

## Value

A tibble with columns: qseqid, sseqid, pident, length, qlen, evalue,
bitscore, qcov
