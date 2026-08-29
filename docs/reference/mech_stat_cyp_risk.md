# Mechanistic-static risk assessment

Mechanistic-static risk assessment

## Usage

``` r
mech_stat_cyp_risk(
  perp,
  cyp_inh,
  cyp_ind,
  cyp_tdi = NULL,
  d = 1,
  include_induction = TRUE,
  substr = cyp_reference_substrates,
  cyp_kdeg = cyp_turnover,
  qh = 1.616,
  qent = 18/60
)
```

## Arguments

- perp:

  Precipitant object.

- cyp_inh:

  Inhibitor object.

- cyp_ind:

  Inducer object.

- cyp_tdi:

  Inhibitor object.

- d:

  Hepatocyte batch scaling factor, defaults to 1.

- include_induction:

  Logical.

- substr:

  The CYP probe substrates to be used as data frame, defaults to
  [`ddir::cyp_reference_substrates`](cyp_reference_substrates.md).

- cyp_kdeg:

  The CYP turnover data as data frame. Defaults to the built-in
  reference data [`ddir::cyp_turnover`](cyp_turnover.md).

- qh:

  Hepatic blood flow in l/min, defaults to 1.616 l/min.

- qent:

  Enteric blood flow in l/min, defaults to 0.3 l/min (18 l/h).

## Value

DDI risk object

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

\$\$I\_{gut}=\frac{Dose}{250\\ mg}\$\$

\\R\\ values \> 10 are considered to indicate a clinical risk for
intestinal CYP3A inhibition.

In the output, the columns `risk_hep` and `risk_intest` indicate whether
the regulatory threshold is reached for the respective enzyme.

## Examples

``` r
mech_stat_cyp_risk(
  examplinib,
  examplinib_cyp_inhibition,
  examplinib_cyp_induction,
  examplinib_cyp_tdi
 )
#> ──────── Clinical DDI risk assessment ──────── 
#> Mechanistic-static risk assessment for examplinib 
#> 
#> object    substrate     kiu    fgut   fm     fmcyp   Ag     Ah     Bg     Bh     Cg     Ch     aucr   risk    
#> CYP1A2    tizanidine    NA     1      0.95   0.98    1      1      1      1      1      1      1      FALSE   
#> CYP2B6    NA            NA     NA     NA     NA      1      1      1      1      1      1      NA     NA      
#> CYP2C8    repaglinide   11     1      1      0.61    0.63   0.98   1      1      1      1      1.01   FALSE   
#> CYP2C9    S-warfarin    0.6    1      1      0.91    0.08   0.76   1      1      1      1      1.28   TRUE    
#> CYP2C19   omeprazole    0.25   1      1      0.87    0.04   0.56   1      1      1      1      1.61   TRUE    
#> CYP2D6    desipramine   NA     1      1      0.85    1      1      1      1      1      1      1      FALSE   
#> CYP3A4    midazolam     12.5   0.57   0.96   1       0.65   0.98   0.43   0.48   6.88   1.77   0.84   FALSE
```
