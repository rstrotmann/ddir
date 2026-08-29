# Plot in vitro CYP induction data

Plot in vitro CYP induction data

## Usage

``` r
induction_plot(x, type = "REL", log = TRUE)
```

## Arguments

- x:

  The in vitro data set as data frame. The following fields are
  expected:

  - DONOR The hepatocyte donor batch as character.

  - PRECIPITANT The DDI precipitant as character.

  - CONC The precipitant concentration.

  - OBJECT The DDI target object as character.

  - FOLD The mRNA or enzyme activity change at the respective
    precipitant concentration.

  - REL The mRNA or enzyme activity induction, as percent of positive
    control at the respective precipitant concentration.

  FOLD and REL are optional but a column that corresponds to the 'type'
  argument must be present.

- type:

  The DDI metric, e.g. FOLD or REL, as character. A different type can
  be given but a column of that name must be present in the data set.

- log:

  Plot logarithmic x axis, as logical.

## Value

A ggplot object.

## Examples

``` r
induction_plot(examplinib_in_vitro_ind)

induction_plot(examplinib_in_vitro_ind, type = "FOLD")
```
