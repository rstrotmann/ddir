# Read transporter inhibition data

Read transporter inhibition data from a file or text connection.

## Usage

``` r
read_transporter_inhibitor_data(source)
```

## Arguments

- source:

  The file or text connection to read from.

## Value

A data frame.

## Details

The following, comma-separated fields are expected (in this order):

- 'name' The perpetrator compound name as character.

- 'cyp' The UGT enzyme as (upper case) character.

- 'ic50' The \\IC\_{50}\\ of the inhibition in μM.

- 'source' Optional source information as character.

Lines starting with '#' are considered comments and are not evaluated.

The following is an example of a valid input:

    # name,     parameter, value, source
    examplinib, Pgp,       0.41,  study 005
    examplinib, BCRP,      1.9,  study 005
    examplinib, OCT1,      2.3,    study 006
    examplinib, OATP1B1,   177,   study 006
    examplinib, OATP1B3,   35,   study 006
    examplinib, OAT1,      271,
    examplinib, OAT3,      300,
    examplinib, BSEP,      12.8,
    examplinib, OCT2,      67,    study 006
    examplinib, MATE1,     3.6,    study 006
    examplinib, MATE2k,    1.1,    study 006

## Examples

``` r
read_transporter_inhibitor_data(textConnection(examplinib_transporter_inhibition_string))
#>          name transporter ic50    source
#> 1  examplinib         Pgp 0.41 study 005
#> 2  examplinib        BCRP  1.9 study 005
#> 3  examplinib        OCT1  2.3 study 006
#> 4  examplinib     OATP1B1  177 study 006
#> 5  examplinib     OATP1B3   35 study 006
#> 6  examplinib        OAT1  271          
#> 7  examplinib        OAT3  300          
#> 8  examplinib        BSEP 12.8          
#> 9  examplinib        OCT2   67 study 006
#> 10 examplinib       MATE1  3.6 study 006
#> 11 examplinib      MATE2k  1.1 study 006
```
