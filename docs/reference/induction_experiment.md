# induction experiment data set constructor

induction experiment data set constructor

## Usage

``` r
induction_experiment(data = NULL, precipitant = "")
```

## Arguments

- data:

  A data frame.

- precipitant:

  The precipitant.

## Value

induction_experiment object.

## Examples

``` r
induction_experiment(examplinib_in_vitro_ind, "examplinib")
#> ──────── CYP induction experiment ──────── 
#> Experiental data for in vitro CYP induction by examplinib 
#> 
#> DONOR   SAMPLE             CONC   OBJECT   FOLD     REL     SOURCE                 
#> A       test               0.01   CYP3A4   0.904    -0.16   Study examplinib_ind   
#> A       test               0.03   CYP3A4   0.94     -0.1    Study examplinib_ind   
#> A       test               0.1    CYP3A4   1.066    0.11    Study examplinib_ind   
#> A       test               0.3    CYP3A4   1.182    0.3     Study examplinib_ind   
#> A       test               0.5    CYP3A4   1.451    0.74    Study examplinib_ind   
#> A       test               1      CYP3A4   1.943    1.55    Study examplinib_ind   
#> A       test               3      CYP3A4   6.331    8.78    Study examplinib_ind   
#> A       test               5      CYP3A4   NA       NA      Study examplinib_ind   
#> A       positive_control   NA     CYP3A4   61.699   100     Study examplinib_ind   
#> B       test               0.01   CYP3A4   0.728    -2.58   Study examplinib_ind   
#> B       test               0.03   CYP3A4   0.881    -1.13   Study examplinib_ind   
#> B       test               0.1    CYP3A4   0.796    -1.93   Study examplinib_ind   
#> B       test               0.3    CYP3A4   0.655    -3.27   Study examplinib_ind   
#> B       test               0.5    CYP3A4   0.751    -2.36   Study examplinib_ind   
#> B       test               1      CYP3A4   0.844    -1.48   Study examplinib_ind   
#> B       test               3      CYP3A4   1.325    3.08    Study examplinib_ind   
#> B       test               5      CYP3A4   2.31     12.41   Study examplinib_ind   
#> B       positive_control   NA     CYP3A4   11.552   100     Study examplinib_ind   
#> C       test               0.01   CYP3A4   0.892    -1.31   Study examplinib_ind   
#> C       test               0.03   CYP3A4   0.796    -2.48   Study examplinib_ind   
#> C       test               0.1    CYP3A4   0.754    -2.99   Study examplinib_ind   
#> C       test               0.3    CYP3A4   0.711    -3.51   Study examplinib_ind   
#> C       test               0.5    CYP3A4   0.648    -4.28   Study examplinib_ind   
#> C       test               1      CYP3A4   0.728    -3.31   Study examplinib_ind   
#> C       test               3      CYP3A4   1.2      2.43    Study examplinib_ind   
#> C       test               5      CYP3A4   1.585    7.11    Study examplinib_ind   
#> C       positive_control   NA     CYP3A4   9.223    100     Study examplinib_ind   
#> A       test               0.3    CYP1A2   1.38     5.8     Study examplinib_ind   
#> A       test               1      CYP1A2   2.07     8.7     Study examplinib_ind   
#> A       test               3      CYP1A2   3.08     12.94   Study examplinib_ind   
#> A       test               10     CYP1A2   7.37     30.97   Study examplinib_ind   
#> A       positive_control   NA     CYP1A2   23.8     100     Study examplinib_ind   
#> B       test               0.3    CYP1A2   1.24     4.28    Study examplinib_ind   
#> B       test               1      CYP1A2   1.94     6.69    Study examplinib_ind   
#> B       test               3      CYP1A2   2.86     9.86    Study examplinib_ind   
#> B       test               10     CYP1A2   5.09     17.55   Study examplinib_ind   
#> B       positive_control   NA     CYP1A2   29       100     Study examplinib_ind   
#> C       test               0.3    CYP1A2   1.56     4.47    Study examplinib_ind   
#> C       test               1      CYP1A2   2.25     6.45    Study examplinib_ind   
#> C       test               3      CYP1A2   3.21     9.2     Study examplinib_ind   
#> C       test               10     CYP1A2   4.15     11.89   Study examplinib_ind   
#> C       positive_control   NA     CYP1A2   34.9     100     Study examplinib_ind
```
