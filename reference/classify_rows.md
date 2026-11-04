# Classify orthogroups by how many assemblies they occur in

Core: present in every column; Unique: present in exactly one;
Accessory: anything in between. Rows with no presence at all (possible
after subsetting) are labelled Absent. Unused levels are dropped so that
the heatmap split only shows classes that occur.

## Usage

``` r
classify_rows(mat)
```

## Arguments

- mat:

  Numeric matrix without NAs.

## Value

A factor with levels in the order Core, Accessory, Unique, Absent.
