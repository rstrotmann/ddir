# Static risk assessment for CYP induction

Static risk assessment for CYP induction

## Usage

``` r
static_cyp_induction_risk(perp, cyp_ind)
```

## Arguments

- perp:

  Precipitant object.

- cyp_ind:

  Induction object.

## Value

DDI risk object

## Examples

``` r
static_cyp_induction_risk(examplinib, examplinib_cyp_induction)
#> ──────── Clinical DDI risk assessment ──────── 
#> Static CYP induction risk for examplinib 
#> 
#> object    emax   max_c   source      risk    note                              
#> CYP1A2    1      5       study 007   FALSE   Not tested up to 50-fold Cmax,u   
#> CYP2B6    1      5       study 007   FALSE   Not tested up to 50-fold Cmax,u   
#> CYP2C8    NA     NA      NA          NA                                        
#> CYP2C9    NA     NA      NA          NA                                        
#> CYP2C19   NA     NA      NA          NA                                        
#> CYP2D6    NA     NA      NA          NA                                        
#> CYP3A4    7.35   3       study 007   TRUE    Not tested up to 50-fold Cmax,u
```
