# Validate a protein_collection before clustering

Checks that there is something to cluster, that protein IDs are safe for
external tools, and that sequences look like protein.

## Usage

``` r
validate_collection_for_clustering(proteins)
```

## Arguments

- proteins:

  A protein_collection object.

## Value

Invisibly `TRUE`.
