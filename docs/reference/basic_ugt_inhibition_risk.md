# UGT inhibition risk

Evaluate the clinical risk for reversible inhibition of UGT enzymes
according to the [ICH M12
guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).

## Usage

``` r
basic_ugt_inhibition_risk(perp, ugt_inh)
```

## Arguments

- perp:

  The precipitant object.

- ugt_inh:

  UGT inhibition data object.

## Value

DDI risk object.

## Details

This function assumes that the UGT inhibition data is provided as
\\IC\_{50}\\. According to [Cheng, Prusoff
1973](https://doi.org/10.1016/0006-2952(73)90196-2)), \\K_i\\ can be
assumed to be \\IC\_{50}/2\\ at the experimental conditions commonly
used in the in vitro inhibition studies where substrate concentrations
are close to \\K_M\\.

The relevant metric for basic modeling of the UGT inhibition risk is
\\R=C\_{max,ss,u}/K\_{i,u}\\

Refer to Section 2.1.2.1 of the [ICH M12 guidance
document](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf))
for details.

\\R\>0.02\\ are considered to indicate a potential UGT inhibition risk.

## Examples

``` r
basic_ugt_inhibition_risk(examplinib, examplinib_ugt_inhibition)
#> ──────── Clinical DDI risk assessment ──────── 
#> UGT inhibition risk for examplinib 
#> 
#> object    ki    kiu   source      r      risk    
#> UGT1A1    15    15    study 009   0.01   FALSE   
#> UGT1A3    15    15    study 009   0.01   FALSE   
#> UGT1A4    15    15    study 009   0.01   FALSE   
#> UGT1A6    15    15    study 009   0.01   FALSE   
#> UGT1A9    3.8   3.8   study 009   0.04   TRUE    
#> UGT2B7    15    15    study 009   0.01   FALSE   
#> UGT2B15   15    15    study 009   0.01   FALSE   
#> UGT2B17   6.1   6.1   study 009   0.03   TRUE
```
