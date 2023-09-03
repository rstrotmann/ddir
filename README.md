
# ddir

<!-- badges: start -->
<!-- badges: end -->

The goal of `ddir` is to facilitate the analysis of the drug-drug interaction risk
for drug compounds as per the relevant FDA and EMA guidelines.

## Installation

You can install the development version of ddir like so:

``` r
library(devtools)

devtools::install_github("rstrotmann/ddir")
```

## Example

The package provides sample data for a fictional drug, `examplinib`. The
following code calculates the DDI risk for direct CYP inhibition (basic method)
and provides a markdown-formatted output table: 

``` r
library(ddir)

p <- examplinib_compounds[[1]]
property_table(p)
basic_cyp_inhibition_risk_table(p, examplinib_cyp_inhibition_data)
```

See also:
* `kinetic_cyp_induction_risk_table()`
* `mech_stat_cyp_risk_table()`
* `basic_ugt_inhibition_risk_table()`
* `transporter_inhibition_risk_table()`
