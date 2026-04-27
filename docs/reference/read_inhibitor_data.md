# Read csv-formatted CYP inhibition data

**\[deprecated\]**

This function is deprecated in favor of
[`read_cyp_inhibitor_data()`](read_cyp_inhibitor_data.md),
[`read_ugt_inhibitor_data()`](read_ugt_inhibitor_data.md) or
[`read_transporter_inhibitor_data()`](read_transporter_inhibitor_data.md).

## Usage

``` r
read_inhibitor_data(source)
```

## Arguments

- source:

  The connection to read from.

## Value

A data frame.

## Details

This function loads CYP inhibition data from a csv file. The expected
fields are (in this order) the compound name, the CYP enzyme, the Ki and
the source information for the respective data. The latter field may
remain empty.

Comment lines must start with '#'.

A valid source is, e.g.,

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
