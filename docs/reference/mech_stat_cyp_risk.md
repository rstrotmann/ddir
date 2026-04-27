# Title

Title

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

  Perpetrator object

- cyp_inh:

  Inhibitor object

- cyp_ind:

  Inducer object

- cyp_tdi:

  Inhibitor object

- d:

  Numeric

- include_induction:

  Logical

- substr:

  Data frame

- cyp_kdeg:

  Data frame

## Value

DDI risk object
