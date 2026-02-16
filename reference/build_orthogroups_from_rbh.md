# Build orthogroups from RBH pairs using connected components

Groups proteins into orthogroups based on RBH relationships. Proteins
connected through any chain of RBH pairs form one orthogroup.

## Usage

``` r
build_orthogroups_from_rbh(rbh_pairs, all_protein_ids)
```

## Arguments

- rbh_pairs:

  A tibble with columns: protein_a, protein_b

- all_protein_ids:

  Character vector of all protein IDs (for identifying singletons)

## Value

A tibble with columns: orthogroup_id, protein_id Only proteins in
orthogroups are returned (singletons excluded)
