# Read UGT inhibition data

Read UGT inhibition data from a file or text connection.

## Usage

``` r
read_ugt_inhibitor_data(source)
```

## Arguments

- source:

  The file or text connection to read from.

## Value

A data frame.

## Details

The following, comma-separated fields are expected in the input:

- 'name' The perpetrator compound name as character.

- 'ugt' The UGT enzyme as (upper case) character.

- 'ic50' The \\IC\_{50}\\ in µM as numeric.

- 'source' Optional source information as character.

Lines starting with '#' are considered comments and are not evaluated.

The following is an example of a valid input:

    # enzyme, IC50, source
    # parent compound

    examplinib, UGT1A1, 15, study 009
    examplinib, UGT1A3, 15, study 009
    examplinib, UGT1A4, 15, study 009
    examplinib, UGT1A6, 15, study 009
    examplinib, UGT1A9, 3.8, study 009
    examplinib, UGT2B7, 15, study 009
    examplinib, UGT2B15, 15, study 009
    examplinib, UGT2B17, 6.1, study 009

    # metabolite M1

    M1, UGT1A1, 1.1, study 009
    M1, UGT1A3, 5.8, study 009
    M1, UGT1A4, 6.2, study 009
    M1, UGT1A6, 15, study 009
    M1, UGT1A9, 3.6, study 009
    M1, UGT2B7, 15, study 009
    M1, UGT2B15, 9.6, study 009
    M1, UGT2B17, 2.2, study 009

## Examples

``` r
read_ugt_inhibitor_data(textConnection(examplinib_ugt_inhibition_string))
#>          name     ugt ic50    source
#> 1  examplinib  UGT1A1   15 study 009
#> 2  examplinib  UGT1A3   15 study 009
#> 3  examplinib  UGT1A4   15 study 009
#> 4  examplinib  UGT1A6   15 study 009
#> 5  examplinib  UGT1A9  3.8 study 009
#> 6  examplinib  UGT2B7   15 study 009
#> 7  examplinib UGT2B15   15 study 009
#> 8  examplinib UGT2B17  6.1 study 009
#> 9          M1  UGT1A1  1.1 study 009
#> 10         M1  UGT1A3  5.8 study 009
#> 11         M1  UGT1A4  6.2 study 009
#> 12         M1  UGT1A6   15 study 009
#> 13         M1  UGT1A9  3.6 study 009
#> 14         M1  UGT2B7   15 study 009
#> 15         M1 UGT2B15  9.6 study 009
#> 16         M1 UGT2B17  2.2 study 009
```
