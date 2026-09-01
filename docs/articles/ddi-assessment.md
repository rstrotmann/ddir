# DDI assessment

This vignette shows how to assess the potential of a drug to act as a
precipitant of drug-drug interactions with CYP and UGT enzymes and drug
transporters.

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
maximal clinical dose, dissolved in a volume of 250 ml (\\I\_{gut}\\).

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
parameters are aggregated into a ‘precipitant’ object. The package
contains a sample precipitant object for the ficitional drug
‘examplinib’:

``` r

examplinib
```

| parameter               | value   | source |
|:------------------------|:--------|:-------|
| oral                    | 1       |        |
| \\MW\\ (g/mol)          | 492.6   |        |
| \\dose\\ (mg)           | 450     |        |
| \\solubility\\ (mg/l)   | Inf     |        |
| \\C\_{max,ss}\\ (ng/ml) | 3530    |        |
| \\f_u\\                 | 0.023   |        |
| \\f\_{u,mic}\\          | 1       |        |
| \\R_B\\                 | 1       |        |
| \\F_a\\                 | 0.81    |        |
| \\F_g\\                 | 1       |        |
| \\k_a\\ (1/min)         | 0.00267 |        |

Precipitant compound parameters for examplinib {.table}

The function [`precipitant()`](../reference/precipitant.md) can be used
to create custom precipitant objects:

``` r

perp <- precipitant(
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

## Precipitant concentrations

The relevant precipitant concentrations for a precipitant compound can
be determined using
[`key_conc_table()`](../reference/key_conc_table.md) - please see the
documentation to this function for details about the calculations.

``` r

key_conc_table(examplinib)
```

| parameter               | value (\\ng/ml\\) | value (\\\mu M\\) |
|:------------------------|------------------:|------------------:|
| \\I\_{gut}\\            |        1800000.00 |           3654.08 |
| \\I\_{max,ss,u}\\       |             81.19 |              0.16 |
| \\I\_{max,inlet,u}\\    |             95.04 |              0.19 |
| \\I\_{max,intestinal}\\ |           3244.05 |              6.59 |

Key precipitant concentrations for examplinib {.table}

## Direct enzyme inhibition

### CYP enzymes

In vitro CYP inhibition data is expected as `inhibition_data` object
that can be constructed like so:

``` r

cyp_inh <- inhibition_data(
  tibble::tribble(
    ~object,  ~ki,     ~source,
   "CYP1A2",   NA,          NA,
   "CYP2B6",   NA,          NA,
   "CYP2C8",   11, "study 001",
   "CYP2C9",  0.6, "study 001",
  "CYP2C19", 0.25, "study 001",
   "CYP2D6",   NA,          NA,
   "CYP3A4", 12.5, "study 001"
  ),
  precipitant = "examplinib"
)
```

Printing this object yields a convenient table view:

``` r

print(cyp_inh)
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

In vitro inhibition data for precipitant examplinib {.table}

The function to assess the clinical risk for direct CYP inhibition,
[`ddir::basic_cyp_inhibition_risk()`](../reference/basic_cyp_inhibition_risk.md),
takes a precipitant object and a CYP inhibition data object:

``` r

print(basic_cyp_inhibition_risk(examplinib, cyp_inh))
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

Direct CYP inhibition risk for examplinib {.table}

### UGT enzymes

``` r

print(basic_ugt_inhibition_risk(examplinib, examplinib_ugt_inhibition))
```

| object  |   ki |  kiu | source    |    r | risk  |
|:--------|-----:|-----:|:----------|-----:|:------|
| UGT1A1  | 15.0 | 15.0 | study 009 | 0.01 | FALSE |
| UGT1A3  | 15.0 | 15.0 | study 009 | 0.01 | FALSE |
| UGT1A4  | 15.0 | 15.0 | study 009 | 0.01 | FALSE |
| UGT1A6  | 15.0 | 15.0 | study 009 | 0.01 | FALSE |
| UGT1A9  |  3.8 |  3.8 | study 009 | 0.04 | TRUE  |
| UGT2B7  | 15.0 | 15.0 | study 009 | 0.01 | FALSE |
| UGT2B15 | 15.0 | 15.0 | study 009 | 0.01 | FALSE |
| UGT2B17 |  6.1 |  6.1 | study 009 | 0.03 | TRUE  |

UGT inhibition risk for examplinib {.table}

## Time-dependent enzyme inhibition

``` r

print(basic_cyp_tdi_risk(examplinib, examplinib_cyp_tdi))
```

| object |   ki |    fu | kinact |   kdeg | source    |    r | risk |
|:-------|-----:|------:|-------:|-------:|:----------|-----:|:-----|
| CYP3A4 | 0.17 | 0.023 |   0.04 | 0.0193 | study 001 | 2.72 | TRUE |

Time-dependent CYP inhibition risk for examplinib {.table}

## CYP induction

### Fold-change method

``` r

print(static_cyp_induction_risk(examplinib, examplinib_cyp_induction))
```

| object  | emax | max_c | source    | risk  | note                            |
|:--------|-----:|------:|:----------|:------|:--------------------------------|
| CYP1A2  | 1.00 |     5 | study 007 | FALSE | Not tested up to 50-fold Cmax,u |
| CYP2B6  | 1.00 |     5 | study 007 | FALSE | Not tested up to 50-fold Cmax,u |
| CYP2C8  |   NA |    NA | NA        | NA    |                                 |
| CYP2C9  |   NA |    NA | NA        | NA    |                                 |
| CYP2C19 |   NA |    NA | NA        | NA    |                                 |
| CYP2D6  |   NA |    NA | NA        | NA    |                                 |
| CYP3A4  | 7.35 |     3 | study 007 | TRUE  | Not tested up to 50-fold Cmax,u |

Static CYP induction risk for examplinib {.table}

### Basic kinetic model

``` r

print(kinetic_cyp_induction_risk(examplinib, examplinib_cyp_induction))
```

| object  | emax | ec50 | max_c | source    |    r | risk |
|:--------|-----:|-----:|------:|:----------|-----:|:-----|
| CYP1A2  | 1.00 |   NA |     5 | study 007 |   NA | NA   |
| CYP2B6  | 1.00 |   NA |     5 | study 007 |   NA | NA   |
| CYP2C8  |   NA |   NA |    NA | NA        |   NA | NA   |
| CYP2C9  |   NA |   NA |    NA | NA        |   NA | NA   |
| CYP2C19 |   NA |   NA |    NA | NA        |   NA | NA   |
| CYP2D6  |   NA |   NA |    NA | NA        |   NA | NA   |
| CYP3A4  | 7.35 | 1.64 |     3 | study 007 | 0.21 | TRUE |

Kinetic CYP induction risk for examplinib {.table}

## Mechanistic-static assessment of CYP-related DDI

``` r

print(mech_stat_cyp_risk(examplinib, examplinib_cyp_inhibition, examplinib_cyp_induction, examplinib_cyp_tdi))
```

| object  | substrate   |   kiu | fgut |   fm | fmcyp |   Ag |   Ah |   Bg |   Bh |   Cg |   Ch | aucr | risk  |
|:--------|:------------|------:|-----:|-----:|------:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|:------|
| CYP1A2  | tizanidine  |    NA | 1.00 | 0.95 |  0.98 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | FALSE |
| CYP2B6  | NA          |    NA |   NA |   NA |    NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |   NA | NA    |
| CYP2C8  | repaglinide | 11.00 | 1.00 | 1.00 |  0.61 | 0.63 | 0.98 | 1.00 | 1.00 | 1.00 | 1.00 | 1.01 | FALSE |
| CYP2C9  | S-warfarin  |  0.60 | 1.00 | 1.00 |  0.91 | 0.08 | 0.76 | 1.00 | 1.00 | 1.00 | 1.00 | 1.28 | TRUE  |
| CYP2C19 | omeprazole  |  0.25 | 1.00 | 1.00 |  0.87 | 0.04 | 0.56 | 1.00 | 1.00 | 1.00 | 1.00 | 1.61 | TRUE  |
| CYP2D6  | desipramine |    NA | 1.00 | 1.00 |  0.85 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | FALSE |
| CYP3A4  | midazolam   | 12.50 | 0.57 | 0.96 |  1.00 | 0.65 | 0.98 | 0.43 | 0.48 | 6.88 | 1.77 | 0.84 | FALSE |

Mechanistic-static risk assessment for examplinib {.table
style="width:100%;"}

## Transporter inhibition

``` r

print(transporter_inhibition_risk(examplinib, examplinib_transporter_inhibition))
```

| object   |   ic50 | source    | i          |       r | threshold | risk  |
|:---------|-------:|:----------|:-----------|--------:|----------:|:------|
| Pgp_int  |   0.41 | study 005 | igut       | 8912.39 |     10.00 | TRUE  |
| Pgp_sys  |   0.41 | study 005 | imaxssu    |    0.40 |      0.02 | TRUE  |
| BCRP_int |   1.90 | study 005 | igut       | 1923.20 |     10.00 | TRUE  |
| BCRP_sys |   1.90 | study 005 | imaxssu    |    0.09 |      0.02 | TRUE  |
| OATP1B1  | 177.00 | study 006 | imaxinletu |    0.00 |      0.10 | FALSE |
| OATP1B3  |  35.00 | study 006 | imaxinletu |    0.01 |      0.10 | FALSE |
| OAT1     | 271.00 | NA        | imaxssu    |    0.00 |      0.10 | FALSE |
| OAT3     | 300.00 | NA        | imaxssu    |    0.00 |      0.10 | FALSE |
| BSEP     |  12.80 | NA        | imaxssu    |    0.01 |      0.10 | FALSE |
| OCT1     |   2.30 | study 006 | NA         |      NA |        NA | NA    |
| OCT2     |  67.00 | study 006 | imaxssu    |    0.00 |      0.10 | FALSE |
| MATE1    |   3.60 | study 006 | imaxssu    |    0.05 |      0.02 | TRUE  |
| MATE2k   |   1.10 | study 006 | imaxssu    |    0.15 |      0.02 | TRUE  |

Transporter inhibition risk for examplinib {.table}
