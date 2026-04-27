# Examplinib CYP inhibition data

This data frame is a typical input to the following functions:
[`basic_cyp_inhibition_risk()`](basic_cyp_inhibition_risk.md)
basic_cyp_inhibition_risk_table
[`mech_stat_cyp_risk()`](mech_stat_cyp_risk.md)
`mech_stat_cyp_risk_table()`

## Usage

``` r
examplinib_cyp_inhibition_data
```

## Format

### `examplinib_cyp_inhibition_data`

A data frame with the columns `name`, `param`, `value` and `source`,
where:

- name is the name of the compound for which the data is recorded

- param contains the respective CYP enzyme names.

- value contains the Ki values for the respective CYP enzyme

- source provides information for the source of the respective value,
  often the name of the DMPK study. This entry is optional.

## Source

Fictional data, made up for demo purposes.

## Details

CYP inhibition data can contain ki data for multiple compounds.

          name   param      value    source
            M1    name         M1
            M1  CYP2C9        4.4 study 002
    examplinib    name examplinib
    examplinib  CYP1A2         NA
    examplinib  CYP2B6         NA
    examplinib  CYP2C8         11 study 001
    examplinib  CYP2C9       13.5 study 001
    examplinib CYP2C19         15 study 001
    examplinib  CYP2D6         NA
    examplinib  CYP3A4       12.5 study 001
