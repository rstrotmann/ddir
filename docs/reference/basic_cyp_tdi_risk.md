# CYP time-dependent inhibition risk

CYP time-dependent inhibition risk

## Usage

``` r
basic_cyp_tdi_risk(perp, cyp_tdi, cyp_kdeg = cyp_turnover)
```

## Arguments

- perp:

  Perpetrator object.

- cyp_tdi:

  inhibition_data object.

- cyp_kdeg:

  CYP turnover rates as data frame, defaults to
  [cyp_turnover](cyp_turnover.md).

## Value

DDI risk object.

## Examples

``` r
basic_cyp_tdi_risk(examplinib, examplinib_cyp_tdi)
#> ──────── Clinical DDI risk assessment ──────── 
#> Time-dependent CYP inhibition risk for examplinib 
#> 
#> object   ki     fu      kinact   kdeg     source      r      risk   
#> CYP3A4   0.17   0.023   0.04     0.0193   study 001   2.72   TRUE
```
