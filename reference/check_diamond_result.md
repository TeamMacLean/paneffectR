# Check the result of a DIAMOND run

Aborts with an informative message when a DIAMOND step timed out or
exited with a non-zero status.

## Usage

``` r
check_diamond_result(result, step, timeout)
```

## Arguments

- result:

  Result of
  [`processx::run()`](http://processx.r-lib.org/reference/run.md).

- step:

  Character name of the DIAMOND subcommand (for messages).

- timeout:

  Timeout that was applied, in seconds.

## Value

Invisibly `TRUE` if the run succeeded.
