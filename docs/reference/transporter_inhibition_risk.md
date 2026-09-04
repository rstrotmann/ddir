# Drug transporter inhibition risk

Drug transporter inhibition risk

## Usage

``` r
transporter_inhibition_risk(
  perp,
  transporter_inh,
  transporter_ref = transporter_reference_data,
  qh = 1.616
)
```

## Arguments

- perp:

  Precipitant object.

- transporter_inh:

  Inhibitor object.

- transporter_ref:

  Data frame.

- qh:

  Hepatic blood flow in l/min, defaults to 1.616 l/min.

## Value

DDI risk object.

## Details

The metric for the assessment of transporter interactions is
\\R=\[I\]/IC\_{50}\\.

The relevant precipitant concentrations \\\[I\]\\ and regulatory
thresholds of concern are:

|  |  |  |
|----|----|----|
| I | transporter | threshold |
| \\I\_{gut}\\ | P-gp and BRCR when drugs are orally administered | 10 |
| \\C\_{max,ss,u}\\ | P-gp and BRCR when drugs are administered parenterally or for drug metabolites | 0.02 |
| \\I\_{max,inlet,u}\\ | hepatic basolateral transporters OCT1, OATP1B1 and OATP1B3 | 0.1 |
| \\C\_{max,ss,u}\\ | renal basolateral transporters OAT1, OAT3 and OCT2 | 0.1 |
| \\C\_{max,ss,u}\\ | apical transporters MATE1 and MATE2-K | 0.02 |

## Examples

``` r
transporter_inhibition_risk(examplinib, examplinib_transporter_inhibition)
#> ──────── Clinical DDI risk assessment ──────── 
#> Transporter inhibition risk for examplinib 
#> 
#> object     ic50   source      i            r          threshold   risk    
#> Pgp_int    0.41   study 005   igut         8910       10          TRUE    
#> Pgp_sys    0.41   study 005   imaxssu      0.402      0.02        TRUE    
#> BCRP_int   1.9    study 005   igut         1920       10          TRUE    
#> BCRP_sys   1.9    study 005   imaxssu      0.0867     0.02        TRUE    
#> OATP1B1    177    study 006   imaxinletu   0.00109    0.1         FALSE   
#> OATP1B3    35     study 006   imaxinletu   0.00551    0.1         FALSE   
#> OAT1       271    NA          imaxssu      0.000608   0.1         FALSE   
#> OAT3       300    NA          imaxssu      0.000549   0.1         FALSE   
#> BSEP       12.8   NA          imaxssu      0.0129     0.1         FALSE   
#> OCT1       2.3    study 006   imaxssu      0.0717     0.1         FALSE   
#> OCT2       67     study 006   imaxssu      0.00246    0.1         FALSE   
#> MATE1      3.6    study 006   imaxssu      0.0458     0.02        TRUE    
#> MATE2k     1.1    study 006   imaxssu      0.15       0.02        TRUE
```
