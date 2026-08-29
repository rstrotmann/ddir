# ddir

Functions for the calculation and documentation of the potential of
drugs to act as perpetrators of clinically relevant drug-drug
interaction (DDI).

## Installation

You can install the development version of ddir like so:

``` r

devtools::install_github("rstrotmann/ddir")
```

## Example

The package provides sample data for a fictional drug, `examplinib`. The
following code shows the relevant drug properties for examplinib, and
calculates the perpetrator DDI risk for direct CYP inhibition (basic
method). All are provided as markdown-formatted output tables:

``` r

library(ddir)
library(knitr)

p <- examplinib_parent
print(p)
```

| parameter               | value   | source        |
|:------------------------|:--------|:--------------|
| oral                    | 1       |               |
| \\MW\\ (g/mol)          | 492.6   |               |
| \\dose\\ (mg)           | 450     | clinical dose |
| \\solubility\\ (mg/l)   | Inf     |               |
| \\C\_{max,ss}\\ (ng/ml) | 3530    | study 001     |
| \\f_u\\                 | 0.023   | study 002     |
| \\f\_{u,mic}\\          | 1       | default       |
| \\R_B\\                 | 1       | study 003     |
| \\F_a\\                 | 0.81    | study 003     |
| \\F_g\\                 | 1       | default       |
| \\k_a\\ (1/min)         | 0.00267 | unknown       |

Perpetrator compound parameters for examplinib

``` r

print(examplinib_cyp_inh_parent)
```

| target  | \\K\_{i}\\ (\\\mu M\\) | source      |
|:--------|-----------------------:|:------------|
| CYP2C8  |                  11.00 | (study 001) |
| CYP2C9  |                   0.60 | (study 001) |
| CYP2C19 |                   0.25 | (study 001) |
| CYP3A4  |                  12.50 | (study 001) |

In vitro inhibition data for examplinib

``` r

print(basic_cyp_inhibition_risk(p, examplinib_cyp_inh_parent))
```

| target | \\K\_{i}\\ (\\\mu M\\) | \\K\_{i,u}\\ (\\\mu M\\) | source | \\R\\ | risk (hepatic) | \\R\_{gut}\\ | risk (intestinal) |
|:---|---:|---:|:---|---:|:---|---:|:---|
| CYP2C8 | 11.00 | 11.00 | study 001 | 0.01 | No | NA |  |
| CYP2C9 | 0.60 | 0.60 | study 001 | 0.27 | Yes | NA |  |
| CYP2C19 | 0.25 | 0.25 | study 001 | 0.66 | Yes | NA |  |
| CYP3A4 | 12.50 | 12.50 | study 001 | 0.01 | No | 292.33 | Yes |

Direct CYP inhibition risk for examplinib

See also:
[`kinetic_cyp_induction_risk()`](reference/kinetic_cyp_induction_risk.md),
[`mech_stat_cyp_risk()`](reference/mech_stat_cyp_risk.md),
[`basic_ugt_inhibition_risk()`](reference/basic_ugt_inhibition_risk.md),
[`transporter_inhibition_risk()`](reference/transporter_inhibition_risk.md)

Full documentation for `ddir` can be found together with the source code
on [github](https://github.com/rstrotmann/ddir).
