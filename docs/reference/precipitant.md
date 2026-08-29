# Precipitant class definition

Precipitant class definition

## Usage

``` r
precipitant(
  name = "",
  oral = TRUE,
  mw = 0,
  dose = 0,
  imaxss = 0,
  fu = 1,
  fumic = 1,
  rb = 1,
  fa = 1,
  fg = 1,
  ka = 0.1,
  solubility = Inf,
  source = character(0)
)
```

## Arguments

- name:

  character.

- oral:

  logical.

- mw:

  numeric.

- dose:

  numeric.

- imaxss:

  numeric.

- fu:

  numeric.

- fumic:

  numeric.

- rb:

  numeric.

- fa:

  numeric.

- fg:

  numeric.

- ka:

  numeric.

- solubility:

  numeric.

- source:

  character

## Value

A precipitant class object.

## Examples

``` r
precipitant(
  name = "examplinib",
  oral = TRUE,
  mw = 492.6,
  dose = 450,
  imaxss = 3530,
  fu = 0.023,
  fumic = 1,
  rb = 1,
  fa = 0.81,
  fg = 1,
  ka = 0.00267,
  solubility = Inf,
  source = c(dose = "clinical dose", imaxss = "study 001", fu = "study 002")
)
#> ──────── DDI precipitant examplinib ──────── 
#> param        value     source          
#> oral         1                         
#> mw           492.6                     
#> dose         450       clinical dose   
#> solubility   Inf                       
#> imaxss       3530      study 001       
#> fu           0.023     study 002       
#> fumic        1                         
#> rb           1                         
#> fa           0.81                      
#> fg           1                         
#> ka           0.00267
```
