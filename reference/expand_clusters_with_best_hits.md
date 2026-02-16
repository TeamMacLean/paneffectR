# Expand orthogroups by adding singletons whose best hit is in a cluster

Iteratively adds singleton proteins to existing orthogroups when their
best BLAST hit is already a member of that orthogroup. This allows
orthogroups to grow beyond the strict RBH pairs.

## Usage

``` r
expand_clusters_with_best_hits(orthogroups, singletons, hits)
```

## Arguments

- orthogroups:

  A tibble with columns: orthogroup_id, protein_id

- singletons:

  Character vector of singleton protein IDs

- hits:

  A tibble with columns: qseqid, sseqid, bitscore

## Value

A list with:

- \$orthogroups: expanded tibble with orthogroup_id, protein_id

- \$singletons: character vector of remaining singletons
