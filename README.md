
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ddir

<!-- badges: start -->
<!-- badges: end -->

This package provides functions for the calculation and documentation of
the potential of drugs to act as perpetrators of clinically relevant
drug-drug interaction (DDI).

## Installation

You can install the development version of ddir like so:

``` r
devtools::install_github("rstrotmann/ddir")
```

## Example

The package provides sample data for a fictional drug, `examplinib`. The
following code shows the relevant drug properties for examplinib, and
calculates the perpetrator DDI risk for direct CYP inhibition (basic
method). Both are provided as markdown-formatted output tables:

``` r
library(ddir)

p <- examplinib_parent
property_table(p)
```

| parameter            | value   | source        |
|:---------------------|:--------|:--------------|
| oral                 | TRUE    |               |
| $MW$ (g/mol)         | 492.6   |               |
| $dose$ (mg)          | 450     | clinical dose |
| $C_{max,ss}$ (ng/ml) | 3530    | study 001     |
| $f_u$                | 0.023   | study 002     |
| $f_{u,mic}$          | 1       | default       |
| $R_B$                | 1       | study 003     |
| $F_a$                | 0.81    | study 003     |
| $F_g$                | 1       | default       |
| $k_a$ (1/min)        | 0.00267 | unknown       |
| $solubility$ (mg/l)  | Inf     | default       |

Compound parameters for examplinib

``` r
basic_cyp_inhibition_risk_table(p, examplinib_cyp_inhibition_data)
```

| CYP     | $K_{i}$ (µM) | $K_{i,u}$ (µM) | $R_1$ | risk (hepatic) | $R_{1,gut}$ | risk (intestinal) |
|:--------|-------------:|---------------:|:------|:---------------|:------------|:------------------|
| CYP1A2  |           NA |             NA | NA    | NA             | NA          | NA                |
| CYP2B6  |           NA |             NA | NA    | NA             | NA          | NA                |
| CYP2C8  |         11.0 |           11.0 | 1.015 | FALSE          | NA          | NA                |
| CYP2C9  |         13.5 |           13.5 | 1.012 | FALSE          | NA          | NA                |
| CYP2C19 |         15.0 |           15.0 | 1.011 | FALSE          | NA          | NA                |
| CYP2D6  |           NA |             NA | NA    | NA             | NA          | NA                |
| CYP3A4  |         12.5 |           12.5 | 1.013 | FALSE          | 293.3       | TRUE              |

Risk for direct CYP inhibition by examplinib (basic model)

See also: `kinetic_cyp_induction_risk_table()`,
`mech_stat_cyp_risk_table()`, `basic_ugt_inhibition_risk_table()`,
`transporter_inhibition_risk_table()`

## Full DDI perpetrator report

As an easy starting point for your own full DDI perpetrator report, load
the sample DDI report from the package source or
[github](https://github.com/rstrotmann/ddir/blob/main/R/DDI-report_examplinib.Rmd),
replace the compound-specific data with the data for your drug, and
`knit()` to a Word or pdf document.

Full documentation for `ddir` can be found together with the source code
on [github](https://github.com/rstrotmann/ddir).
