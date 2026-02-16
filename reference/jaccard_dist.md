# Calculate Jaccard distance

Jaccard dissimilarity: 1 - \|A ∩ B\| / \|A ∪ B\| For binary data: 1 -
(number of shared 1s) / (number where either is 1)

## Usage

``` r
jaccard_dist(mat)
```

## Arguments

- mat:

  Numeric matrix (rows are observations)

## Value

A dist object
