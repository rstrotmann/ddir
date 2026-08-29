# DDI assessment

This vignette shows how to assess the potential of a drug to act as a
perpetrator of drug-drug interactions with CYP and UGT enzymes and drug
transporters. The live code contained in this document relies on the
following R packages:

``` r
library(ddir)
library(tidyverse)
```

## Regulatory background

In May-2024, ICH published a new harmonized guidance on the assessment
of enzyme- or transporter-mediated drug interactions that will be
adopted by the ICH-abiding regulatory agencies, including FDA, EMA and
PMDA.

Currently (Nov-2024), the following guidance documents are provided by
the individual regulatory authorities:

- [M12 Drug Interaction Studies (FDA,
  Aug-2024)](https://www.fda.gov/media/161199/download)
- [M12 Questions and Answers document
  (Aug-2024)](https://www.fda.gov/media/180488/download?attachment)
- [ICH M12 Drug Interaction Studies (EMA,
  May-2024)](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf)

## Drug properties

As a general concept, the risk for enzyme or transporter interactions is
evaluated for a given drug exposure level that corresponds to the
maximal unbound clinical exposure in the relevant pharmacokinetic
compartment. For interactions with hepatic enzymes, the unbound maximal
plasma concentration at steady state (\\I\_{max,ss,u}\\) is considered,
and for intestinal enzymes as victim of orally administered drugs, the
maximal clinical dose, dissolved in a volume of 250 mg (\\I\_{gut}\\).

For interactions with basolateral hepatic transporters (OATP1B1,
OATP1B3), the unbound hepatic inlet concentration (\\I\_{max,inlet,u}\\)
is considered that is composed of the \\I\_{max,ss,u}\\ and a portal
venous term reflecting the intestinally absorbed drug escaping gut
metabolism.

These concentrations can be automatically derived from the following set
of drug-specific parameters:

| Parameter                                      | Parameter name | Default | Unit  |
|:-----------------------------------------------|:--------------:|:-------:|:-----:|
| Molar weight                                   |       mw       |   NA    | g/mol |
| Clinical dose                                  |      dose      |   NA    |  mg   |
| Maximal total plasma concentration (\\C_max\\) |     imaxss     |   NA    | ng/ml |
| Fraction unbound (\\f_u\\)                     |       fu       |    1    |       |
| Microsomal unbound fraction (\\f\_{u,mic}\\)   |     fumic      |    1    |       |
| Blood-to-plasma concentration ratio (\\R_b\\)  |       rb       |    1    |       |
| Fraction absorbed (\\f_a\\)                    |       fa       |    1    |       |
| Fraction escapting gut metabolism (\\f_g\\)    |       fg       |    1    |       |
| Absorption rate constant (\\k_a\\)             |       ka       |   0.1   | 1/min |
| Solubility                                     |   solubility   |   Inf   | mg/l  |

The specification of the key compound parameters is therefore the first
step in the formal analysis. In the scope of this package, these key
parameters are aggregated into a ‘perpetrator’ object. The package
contains a sample perpetrator object for the ficitional drug
‘examplinib’:

``` r
examplinib_parent
#> ----- DDI perpetrator object -----
#> name         examplinib                       
#> oral         TRUE                             
#> mw           492.6 g/mol                      
#> dose         450 mg         (clinical dose)   
#> imaxss       3530 ng/ml     (study 001)       
#> fu           0.023          (study 002)       
#> fumic        1              (default)         
#> rb           1              (study 003)       
#> fa           0.81           (study 003)       
#> fg           1              (default)         
#> ka           0.00267 /min   (unknown)         
#> solubility   Inf mg/ml
```

The `ddir` package provides multiple options to create perpetrator
objects. In the simplest case, the function `make_perpetrator()` can be
used:

``` r
perp <- perpetrator(
  name = "test",
  dose = 100,
  imaxss = 1000,
  mw = 500,
  oral = TRUE,
  fu = 1,
  fumic = 1,
  rb = 1,
  fa = 1,
  fg = 1,
  ka = 0.1,
  solubility = Inf
)
```

Other options are presented in a dedicated vignette (see
`vignette("load-data")`).

To facilitate the generation of R notebook-based reports,
`property_table()` can be used to print a markdown-formatted table of
the key compound parameters:

``` r
print(examplinib_parent)
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

## Perpetrator concentrations

The relevant perpetrator concentrations for a perpetrator compound can
be determined using [`key_conc()`](../reference/key_conc.md) - please
see the documentation to this function for details about the
calculations.

``` r
key_conc(examplinib_parent)
```

| parameter               | value (\\ng/ml\\) | value (\\\mu M\\) |
|:------------------------|------------------:|------------------:|
| \\I\_{gut}\\            |        1800000.00 |           3654.08 |
| \\I\_{max,ss,u}\\       |             81.19 |              0.16 |
| \\I\_{max,inlet,u}\\    |             95.04 |              0.19 |
| \\I\_{max,intestinal}\\ |           3244.05 |              6.59 |

Key perpetrator concentrations for examplinib

## Direct enzyme inhibition

### CYP enzymes

``` r
print(basic_cyp_inhibition_risk(examplinib_parent, examplinib_cyp_inh_parent))
```

| target | \\K\_{i}\\ (\\\mu M\\) | \\K\_{i,u}\\ (\\\mu M\\) | source | \\R\\ | risk (hepatic) | \\R\_{gut}\\ | risk (intestinal) |
|:---|---:|---:|:---|---:|:---|---:|:---|
| CYP2C8 | 11.00 | 11.00 | study 001 | 0.01 | No | NA |  |
| CYP2C9 | 0.60 | 0.60 | study 001 | 0.27 | Yes | NA |  |
| CYP2C19 | 0.25 | 0.25 | study 001 | 0.66 | Yes | NA |  |
| CYP3A4 | 12.50 | 12.50 | study 001 | 0.01 | No | 292.33 | Yes |

Direct CYP inhibition risk for examplinib

### UGT enzymes

``` r
print(basic_ugt_inhibition_risk(examplinib_parent, examplinib_ugt_inh_parent))
```

| target  | \\K\_{i}\\ (\\\mu M\\) | \\K\_{i,u}\\ (\\\mu M\\) | source    | \\R\\ | risk |
|:--------|-----------------------:|-------------------------:|:----------|------:|:-----|
| UGT1A1  |                   7.50 |                     7.50 | study 009 |  0.02 | Yes  |
| UGT1A3  |                   7.50 |                     7.50 | study 009 |  0.02 | Yes  |
| UGT1A4  |                   7.50 |                     7.50 | study 009 |  0.02 | Yes  |
| UGT1A6  |                   7.50 |                     7.50 | study 009 |  0.02 | Yes  |
| UGT1A9  |                   1.90 |                     1.90 | study 009 |  0.09 | Yes  |
| UGT2B7  |                   7.50 |                     7.50 | study 009 |  0.02 | Yes  |
| UGT2B15 |                   7.50 |                     7.50 | study 009 |  0.02 | Yes  |
| UGT2B17 |                   3.05 |                     3.05 | study 009 |  0.05 | Yes  |

UGT inhibition risk for examplinib

## Time-dependent enzyme inhibition

``` r
print(basic_cyp_tdi_risk(examplinib_parent, examplinib_cyp_tdi_parent))
```

| target | \\K\_{i}\\ (\\\mu M\\) | \\f_u\\ | \\k\_{inact}\\ (1/h) | \\k\_{deg}\\ (1/h) | source | \\R\\ | risk |
|:---|---:|---:|---:|---:|:---|---:|:---|
| CYP3A4 | 0.17 | 0.023 | 0.04 | 0.0193 | study 001 | 2.72 | Yes |

Time-dependent CYP inhibition risk for examplinib

## CYP induction

### Fold-change method

``` r
print(static_cyp_induction_risk(examplinib_parent, examplinib_cyp_ind_parent))
```

| target | \\E\_{max}\\ | maxc | source | \\max c\\ (\\\mu M\\) | risk | note |
|:---|---:|---:|:---|---:|:---|:---|
| CYP1A2 | 1.00 | 5 | study 007 | 5 | No | Not tested up to 50-fold Cmax,u |
| CYP2B6 | 1.00 | 5 | study 007 | 5 | No | Not tested up to 50-fold Cmax,u |
| CYP3A4 | 7.35 | 3 | study 007 | 3 | Yes | Not tested up to 50-fold Cmax,u |

Static CYP induction risk for examplinib

### Basic kinetic model

``` r
print(kinetic_cyp_induction_risk(examplinib_parent, examplinib_cyp_ind_parent))
```

| target | \\E\_{max}\\ | \\EC\_{50}\\ | maxc | source | \\max c\\ (\\\mu M\\) | \\R\\ | risk |
|:---|---:|---:|---:|:---|---:|---:|:---|
| CYP1A2 | 1.00 | NA | 5 | study 007 | 5 | NA |  |
| CYP2B6 | 1.00 | NA | 5 | study 007 | 5 | NA |  |
| CYP3A4 | 7.35 | 1.64 | 3 | study 007 | 3 | 0.21 | Yes |

Kinetic CYP induction risk for examplinib

## Mechanistic-static assessment of CYP-related DDI

``` r
print(mech_stat_cyp_risk(examplinib_parent, examplinib_cyp_inh_parent,
                   examplinib_cyp_ind_parent))
```

| target | substrate | \\K\_{i,u}\\ (\\\mu M\\) | \\f\_{gut}\\ | fm | \\fm\_{CYP}\\ | Ag | Ah | Bg | Bh | Cg | Ch | \\R\_{AUC}\\ | risk |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| CYP2C8 | repaglinide | 11.00 | 1.00 | 1.00 | 0.61 | 0.63 | 0.98 | 1 | 1 | 1.00 | 1.00 | 1.01 | No |
| CYP2C9 | S-warfarin | 0.60 | 1.00 | 1.00 | 0.91 | 0.08 | 0.76 | 1 | 1 | 1.00 | 1.00 | 1.28 | Yes |
| CYP2C19 | omeprazole | 0.25 | 1.00 | 1.00 | 0.87 | 0.04 | 0.56 | 1 | 1 | 1.00 | 1.00 | 1.61 | Yes |
| CYP3A4 | midazolam | 12.50 | 0.57 | 0.96 | 1.00 | 0.65 | 0.98 | 1 | 1 | 6.88 | 1.77 | 0.23 | Yes |

Mechanistic-static CYP inhibition risk for examplinib

## Transporter inhibition

``` r
print(transporter_inhibition_risk(examplinib_parent, examplinib_transporter_inh_parent))
```

| target   | \\IC\_{50}\\ | source    | \\i\\      |   \\R\\ | threshold | risk |
|:---------|-------------:|:----------|:-----------|--------:|----------:|:-----|
| Pgp_int  |         0.41 | study 005 | igut       | 8912.39 |     10.00 | Yes  |
| Pgp_sys  |         0.41 | study 005 | imaxssu    |    0.40 |      0.02 | Yes  |
| BCRP_int |         1.90 | study 005 | igut       | 1923.20 |     10.00 | Yes  |
| BCRP_sys |         1.90 | study 005 | imaxssu    |    0.09 |      0.02 | Yes  |
| OATP1B1  |       177.00 | study 006 | imaxinletu |    0.00 |      0.10 | No   |
| OATP1B3  |        35.00 | study 006 | imaxinletu |    0.01 |      0.10 | No   |
| OAT1     |       271.00 |           | imaxssu    |    0.00 |      0.10 | No   |
| OAT3     |       300.00 |           | imaxssu    |    0.00 |      0.10 | No   |
| BSEP     |        12.80 |           | imaxssu    |    0.01 |      0.10 | No   |
| OCT1     |         2.30 | study 006 | NA         |      NA |        NA |      |
| OCT2     |        67.00 | study 006 | imaxssu    |    0.00 |      0.10 | No   |
| MATE1    |         3.60 | study 006 | imaxssu    |    0.05 |      0.02 | Yes  |
| MATE2k   |         1.10 | study 006 | imaxssu    |    0.15 |      0.02 | Yes  |

Transporter inhibition risk for examplinib
