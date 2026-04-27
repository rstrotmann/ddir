# Read csv-formatted CYP TDI data

This function reads comma-separated CYP TDI data from a file or text
connection.

## Usage

``` r
read_tdi_data(source)
```

## Arguments

- source:

  The connection to read from.

## Value

A data frame.

## Details

The following fields are expected (in this order):

- 'name' The perpetrator compound name as character.

- 'cyp' The CYP enzyme as (upper case) character.

- 'ki' The \\K_I\\ in µM as numeric.

- 'kinact' The \\k\_{inact}\\ in \\1/h\\ as numeric.

- 'source' Optional source information as character. Lines starting with
  '#' are interpreted as comments and are not evaluated.

A valid data set is, e.g.,

          # compound, CYP, Ki, Kinact, source
          examplinib, CYP3A4, 0.17,   0.04, study 001
