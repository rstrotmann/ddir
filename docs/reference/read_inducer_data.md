# Read csv-formatted CYP inducer data

This function loads CYP inducer data from a csv file. The following
fields are expected:

- 'name' The name of the perpetrator compound as character.

- 'cyp' The CYP enzyme as (upper case) character.

- 'emax' The \\E\_{max}\\, i.e., the maximum induction effect determined
  in vitro as numeric.

- 'ec50' The \\EC\_{50}\\ in µM as numeric.

- 'maxc' The maximal concentration in µM tested in the in vitro assay as
  numeric.

- 'source' Optional source information as character.

## Usage

``` r
read_inducer_data(source)
```

## Arguments

- source:

  The connection to read from.

## Value

The data as data frame.

## Details

Comment lines must start with '#'.

A valid source is, e.g.,

    examplinib, CYP1A2,  1,    NA,   5,  study 007
    examplinib, CYP2B6,  1,    NA,   5,  study 007
    examplinib, CYP2C8,  NA,   NA,   NA,
    examplinib, CYP2C9,  NA,   NA,   NA,
    examplinib, CYP2C19, NA,   NA,   NA,
    examplinib, CYP2D6,  NA,   NA,   NA,
    examplinib, CYP3A4,   7.35, 1.64, 3,  study 007


    M1, CYP1A2,  1,    NA,   5,  study 007
    M1, CYP2B6,  6.98, 1.86, 5,  study 007
    M1, CYP2C8,  NA,   NA,   NA,
    M1, CYP2C9,  NA,   NA,   NA,
    M1, CYP2C19, NA,   NA,   NA,
    M1, CYP2D6,  NA,   NA,   NA,
    M1, CYP3A4,  22.7, 1.1,  5,  study 007

## Examples

``` r
read_inducer_data(textConnection(examplinib_cyp_induction_string))
#>          name     cyp  emax ec50 maxc    source
#> 1  examplinib  CYP1A2  1.00   NA    5 study 007
#> 2  examplinib  CYP2B6  1.00   NA    5 study 007
#> 3  examplinib  CYP2C8    NA   NA   NA          
#> 4  examplinib  CYP2C9    NA   NA   NA          
#> 5  examplinib CYP2C19    NA   NA   NA          
#> 6  examplinib  CYP2D6    NA   NA   NA          
#> 7  examplinib  CYP3A4  7.35 1.64    3 study 007
#> 8          M1  CYP1A2  1.00   NA    5 study 007
#> 9          M1  CYP2B6  6.98 1.86    5 study 007
#> 10         M1  CYP2C8    NA   NA   NA          
#> 11         M1  CYP2C9    NA   NA   NA          
#> 12         M1 CYP2C19    NA   NA   NA          
#> 13         M1  CYP2D6    NA   NA   NA          
#> 14         M1  CYP3A4 22.70 1.10    5 study 007
```
