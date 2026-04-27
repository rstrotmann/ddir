# perpetrator class constructor

perpetrator class constructor

## Usage

``` r
perpetrator(
  name,
  oral,
  mw,
  dose,
  imaxss,
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

  Character.

- oral:

  Oral administration, as logical.

- mw:

  Molecular weight in g/mol.

- dose:

  Administered dose in mg.

- imaxss:

  Cmax at steady state in ng/ml.

- fu:

  Fraction unbound.

- fumic:

  Fraction unbound in microsomes.

- rb:

  blood cell distribution.

- fa:

  Fraction absorbed.

- fg:

  Fraction escaping gut metabolism.

- ka:

  Absorption constant in 1/min.

- solubility:

  Aqueous solubility in mg/l.

- source:

  Source information for parameters as named character vector

## Value

perpetrator object

## Examples

``` r
perpetrator(
"examplinib", TRUE, 492.6, 450, 3530, fu = 0.023, fa = 0.81, ka = .00267,
source = c(dose = "clinical dose", imaxss = "study 001", fu = "study 002")
)
#> ----- DDI perpetrator object -----
#> name         examplinib                       
#> oral         TRUE                             
#> mw           492.6 g/mol                      
#> dose         450 mg         (clinical dose)   
#> imaxss       3530 ng/ml     (study 001)       
#> fu           0.023          (study 002)       
#> fumic        1                                
#> rb           1                                
#> fa           0.81                             
#> fg           1                                
#> ka           0.00267 /min                     
#> solubility   Inf mg/ml
```
