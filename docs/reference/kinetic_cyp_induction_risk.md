# Kinetic assessment of CYP induction risk

Kinetic assessment of CYP induction risk

## Usage

``` r
kinetic_cyp_induction_risk(perp, cyp_ind, d = 1)
```

## Arguments

- perp:

  Precipitant object

- cyp_ind:

  induction_data object.

- d:

  Scaling factor, defaults to 1.

## Value

DDI risk object

## Examples

``` r
kinetic_cyp_induction_risk(examplinib, examplinib_cyp_induction)
#> ──────── Clinical DDI risk assessment ──────── 
#> Kinetic CYP induction risk for examplinib 
#> 
#> object    emax   ec50   max_c   source      r       risk   
#> CYP1A2    1      NA     5       study 007   NA      NA     
#> CYP2B6    1      NA     5       study 007   NA      NA     
#> CYP2C8    NA     NA     NA      NA          NA      NA     
#> CYP2C9    NA     NA     NA      NA          NA      NA     
#> CYP2C19   NA     NA     NA      NA          NA      NA     
#> CYP2D6    NA     NA     NA      NA          NA      NA     
#> CYP3A4    7.35   1.64   3       study 007   0.213   TRUE
```
