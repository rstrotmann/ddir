# Basic CYP inhibition risk

Evaluate the clinical risk for direct (reversible) CYP inhibition
according to the basic model defined in the [ICH M12
guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).

## Usage

``` r
basic_cyp_inhibition_risk(perp, cyp_inh)
```

## Arguments

- perp:

  The precipitant object.

- cyp_inh:

  CYP inhibition data object.

## Value

DDI risk object.

## Details

For the basic modeling of direct (reversible) CYP enzyme inhibition, the
ratio of the relevant inhibitor concentration to the \\K_i\\ of the
respective CYP enzyme is considered, i.e., \\R\\ for hepatic enzymes and
\\R\_{gut}\\ for intestinal enzymes (refer to Section 2.1.2.1 of the
[ICH M12 guidance
document](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf)).

### Liver

\$\$R=\frac{C\_{max,ss,u}}{K\_{i,u}}\$\$

\\R\\ values \> 0.02, i.e., maximal unbound precipitant concentrations
50-fold over \\K_i\\ are considered to indicate a potential clinical CYP
inhibition risk using this method.

### Intestine

\$\$R\_{gut}=\frac{I\_{gut}}{K\_{i,u}}\$\$

where

\$\$I\_{gut}=\frac{Dose}{250\\ ml}\$\$

\\R\\ values \> 10 are considered to indicate a clinical risk for
intestinal CYP3A inhibition.

In the output, the columns `risk_hep` and `risk_intest` indicate whether
the regulatory threshold is reached.

## Examples

``` r
basic_cyp_inhibition_risk(examplinib, examplinib_cyp_inhibition)
#> ──────── Clinical DDI risk assessment ──────── 
#> Direct CYP inhibition risk for examplinib 
#> 
#> object    ki     kiu    source      r      risk_hep   r_gut      risk_intest   
#> CYP1A2    NA     NA     NA          NA     NA         NA         NA            
#> CYP2B6    NA     NA     NA          NA     NA         NA         NA            
#> CYP2C8    11     11     study 001   0.01   FALSE      NA         NA            
#> CYP2C9    0.6    0.6    study 001   0.27   TRUE       NA         NA            
#> CYP2C19   0.25   0.25   study 001   0.66   TRUE       NA         NA            
#> CYP2D6    NA     NA     NA          NA     NA         NA         NA            
#> CYP3A4    12.5   12.5   study 001   0.01   FALSE      292.3264   TRUE
```
