# Load name mapping from CSV

Parses a name mapping CSV that maps original protein IDs to standardized
names.

## Usage

``` r
load_name_mapping(path)
```

## Arguments

- path:

  Character. Path to name mapping CSV file.

## Value

A tibble with columns: assembly, old_name, new_name, protein_number.
