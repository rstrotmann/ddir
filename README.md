# ddir

<!-- badges: start -->
<!-- badges: end -->

This package provides functions for the analysis of the clinical drug-drug
interaction (DDI) risk as per the relevant FDA and EMA guidelines.

## Installation

You can install the development version of `ddir` using:

``` r
devtools::install_github("rstrotmann/ddir")
```

## Example

The package provides sample data for a fictional drug, `examplinib`. The
following code calculates the perpetrator DDI risk for direct CYP inhibition
(basic method) and provides a markdown-formatted output table: 

``` r
library(ddir)

p <- examplinib_compounds[[1]]
property_table(p)
basic_cyp_inhibition_risk_table(p, examplinib_cyp_inhibition_data)
```

See also:
`kinetic_cyp_induction_risk_table()`,
`mech_stat_cyp_risk_table()`,
`basic_ugt_inhibition_risk_table()`,
`transporter_inhibition_risk_table()`

Full documentation can be found together with the source code on
[github](https://github.com/rstrotmann/ddir).
