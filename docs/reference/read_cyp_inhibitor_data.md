# Read CYP inhibition data

Read CYP inhibition data from a file or text connection.

## Usage

``` r
read_cyp_inhibitor_data(source)
```

## Arguments

- source:

  The file or text connection to read from.

## Value

A data frame.

## Details

The following, comma-separated fields are expected in the input (in this
order):

- 'name' The perpetrator compound name as character.

- 'cyp' The CYP enzyme as (upper case) character.

- 'ki' The \\k_i\\ in µM as numeric.

- 'source' Optional source information as character.

Lines starting with '#' are considered comments and are not evaluated.

The following is an example of a valid input:

        # PARENT
        # compound, CYP, ki, source
        examplinib, CYP1A2,  NA
        examplinib, CYP2B6,  NA,
        examplinib, CYP2C8,  11,   study 001
        examplinib, CYP2C9,  0.6, study 001
        examplinib, CYP2C19, 0.25,   study 001
        examplinib, CYP2D6,  NA,
        examplinib, CYP3A4,  12.5, study 001

        # METABOLITE
        M1,         CYP2C9,  4.4,  study 002

## Examples

``` r
read_cyp_inhibitor_data(textConnection(examplinib_cyp_inhibition_string))
#>         name     cyp   ki    source
#> 1 examplinib  CYP1A2   NA          
#> 2 examplinib  CYP2B6   NA          
#> 3 examplinib  CYP2C8   11 study 001
#> 4 examplinib  CYP2C9  0.6 study 001
#> 5 examplinib CYP2C19 0.25 study 001
#> 6 examplinib  CYP2D6   NA          
#> 7 examplinib  CYP3A4 12.5 study 001
#> 8         M1  CYP2C9  4.4 study 002
```
