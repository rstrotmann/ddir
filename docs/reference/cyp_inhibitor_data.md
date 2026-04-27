# Make CYP inhibition table

Make CYP inhibition table

## Usage

``` r
cyp_inhibitor_data(target, ki, source = character(0))
```

## Arguments

- target:

  Name of the target.

- ki:

  Named list.

- source:

  Named list.

## Value

A data table.

## Examples

``` r
cyp_inhibitor_data(
"examplinib",
c(CYP3A4 = 1, CYP1A2 = 2, CYP2D6 = 3),
source = c(CYP3A4 = "source 1", CYP1A2 = "source 2"))
#> Error in enframe(source, name = "cyp", value = "source"): could not find function "enframe"
```
