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

  The perpetrator object.

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

\\R\\ values \> 0.02, i.e., maximal unbound perpetrator concentrations
50-fold over \\K_i\\ are considered to indicate a potential clinical CYP
inhibition risk using this method.

### Intestine

\$\$R\_{gut}=\frac{I\_{gut}}{K\_{i,u}}\$\$

where

\$\$I\_{gut}=\frac{Dose}{250\\ mg}\$\$

\\R\\ values \> 10 are considered to indicate a clinical risk for
intestinal CYP3A inhibition.

In the output, the columns `risk_hep` and `risk_intest` indicate whether
the regulatory threshold is reached.
