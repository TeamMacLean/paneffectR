# Validate and process annotation data frame

Validate and process annotation data frame

## Usage

``` r
validate_annotation(annotation, expected_names, type)
```

## Arguments

- annotation:

  Data frame with row names

- expected_names:

  Character vector of expected row names

- type:

  "row" or "column" for error messages

## Value

Processed annotation data frame (reordered to match expected_names, with
NA for missing)
