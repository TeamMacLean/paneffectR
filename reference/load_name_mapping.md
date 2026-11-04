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

## Examples

``` r
csv <- system.file("testdata", "assembly1_name_mapping.csv", package = "paneffectR")
mapping <- load_name_mapping(csv)
head(mapping)
#> # A tibble: 6 × 4
#>   assembly  old_name                                     new_name protein_number
#>   <chr>     <chr>                                        <chr>             <int>
#> 1 assembly1 _NODE_23_length_1941_cov_4.403433_000001.1   assembl…             39
#> 2 assembly1 _NODE_24_length_1937_cov_6.784946_000001.1   assembl…             40
#> 3 assembly1 _NODE_77_length_1510_cov_5.166085_000001.1   assembl…             94
#> 4 assembly1 _NODE_105_length_1378_cov_343.911606_000001… assembl…            122
#> 5 assembly1 _NODE_152_length_1247_cov_3.528205_000001.1  assembl…            164
#> 6 assembly1 _NODE_183_length_1197_cov_3.154464_000001.1  assembl…            194
```
