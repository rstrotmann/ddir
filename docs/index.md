# ddir

Functions for the calculation and documentation of the potential of
drugs to act as precipitants of clinically relevant drug-drug
interaction (DDI).

## Installation

You can install the development version of ddir like so:

``` r

devtools::install_github("rstrotmann/ddir")
```

## Example

The package provides sample data for a fictional drug, `examplinib`. The
following code shows the relevant drug properties for examplinib, and
calculates the precipitant DDI risk for direct CYP inhibition (basic
method). All are provided as markdown-formatted output tables:

``` r

library(ddir)
library(knitr)

p <- examplinib
print(p)
#> 
#> 
#> Table: Precipitant compound parameters for examplinib
#> 
#> |parameter            |value   |source |
#> |:--------------------|:-------|:------|
#> |oral                 |1       |       |
#> |$MW$ (g/mol)         |492.6   |       |
#> |$dose$ (mg)          |450     |       |
#> |$solubility$ (mg/l)  |Inf     |       |
#> |$C_{max,ss}$ (ng/ml) |3530    |       |
#> |$f_u$                |0.023   |       |
#> |$f_{u,mic}$          |1       |       |
#> |$R_B$                |1       |       |
#> |$F_a$                |0.81    |       |
#> |$F_g$                |1       |       |
#> |$k_a$ (1/min)        |0.00267 |       |
print(examplinib_cyp_inhibition)
```

| object  |    ki | source    |
|:--------|------:|:----------|
| CYP1A2  |    NA |           |
| CYP2B6  |    NA |           |
| CYP2C8  | 11.00 | study 001 |
| CYP2C9  |  0.60 | study 001 |
| CYP2C19 |  0.25 | study 001 |
| CYP2D6  |    NA |           |
| CYP3A4  | 12.50 | study 001 |

In vitro inhibition data for precipitant examplinib

``` r

print(basic_cyp_inhibition_risk(p, examplinib_cyp_inhibition))
```

| object  |    ki |   kiu | source    |    r | risk_hep |    r_gut | risk_intest |
|:--------|------:|------:|:----------|-----:|:---------|---------:|:------------|
| CYP1A2  |    NA |    NA | NA        |   NA | NA       |       NA | NA          |
| CYP2B6  |    NA |    NA | NA        |   NA | NA       |       NA | NA          |
| CYP2C8  | 11.00 | 11.00 | study 001 | 0.01 | FALSE    |       NA | NA          |
| CYP2C9  |  0.60 |  0.60 | study 001 | 0.27 | TRUE     |       NA | NA          |
| CYP2C19 |  0.25 |  0.25 | study 001 | 0.66 | TRUE     |       NA | NA          |
| CYP2D6  |    NA |    NA | NA        |   NA | NA       |       NA | NA          |
| CYP3A4  | 12.50 | 12.50 | study 001 | 0.01 | FALSE    | 292.3264 | TRUE        |

Direct CYP inhibition risk for examplinib

See also:
[`kinetic_cyp_induction_risk()`](reference/kinetic_cyp_induction_risk.md),
[`mech_stat_cyp_risk()`](reference/mech_stat_cyp_risk.md),
[`basic_ugt_inhibition_risk()`](reference/basic_ugt_inhibition_risk.md),
[`transporter_inhibition_risk()`](reference/transporter_inhibition_risk.md)

Full documentation for `ddir` can be found together with the source code
on [github](https://github.com/rstrotmann/ddir).
