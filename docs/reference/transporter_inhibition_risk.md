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

  Perpetrator object.

- transporter_inh:

  Inhibitor object.

- transporter_ref:

  Data frame.

## Value

DDI risk object.

## Details

The metric for the assessment of transporter interactions is
\\R=\[I\]/IC\_{50}\\.

The relevant perpetrator concentrations \\\[I\]\\ and regulatory
thresholds of concern are:

|  |  |  |
|----|----|----|
| I | transporter | threshold |
| \\I\_{gut}\\ | P-gp and BRCR when drugs are orally administered | 10 |
| \\C\_{max,ss,u}\\ | P-gp and BRCR when drugs are administered parenterally or for drug metabolites | 0.02 |
| \\I\_{max,inlet,u}\\ | hepatic basolateral transporters OCT1, OATP1B1 and OATP1B3 | 0.1 |
| \\C\_{max,ss,u}\\ | renal basolateral transporters OAT1, OAT3 and OCT2 | 0.1 |
| \\C\_{max,ss,u}\\ | apical transporters MATE1 and MATE2-K | 0.02 |
