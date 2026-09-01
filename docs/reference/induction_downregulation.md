# Flag concentration-dependent CYP mRNA down-regulation

**\[experimental\]**

For each CYP isoform (and, by default, each hepatocyte donor), decide
whether vehicle-normalized mRNA fold-change falls as precipitant
concentration rises, rather than being induced.

## Usage

``` r
induction_downregulation(x, fold_threshold = 0.5, by_donor = TRUE)
```

## Arguments

- x:

  An `induction_experiment` object. A `FOLD` column is required.

- fold_threshold:

  Maximum fold-change at the highest concentration that still counts as
  suppression, as numeric. Defaults to 0.5.

- by_donor:

  If `TRUE` (default), one row per isoform and donor. If `FALSE`, donors
  are pooled and there is one row per isoform.

## Value

A data frame with one row per isoform, or per isoform and donor, and the
columns:

- `OBJECT` CYP isoform.

- `DONOR` Hepatocyte donor (omitted if `by_donor = FALSE`).

- `n` Number of concentrations used.

- `min_fold`, `max_fold` Range of `FOLD`.

- `fold_maxc` `FOLD` at the highest `CONC`.

- `rho` Spearman correlation of `CONC` and `FOLD`.

- `downregulation` `TRUE` if the curve meets both criteria above.

## Details

Only `SAMPLE == "test"` rows with `CONC > 0` and non-missing `FOLD` are
used. A curve is flagged when both of the following hold:

- fold-change at the highest tested concentration is below
  `fold_threshold` (default 0.5, i.e. ≤ 50% of vehicle);

- Spearman's rank correlation of `CONC` with `FOLD` is negative, so
  fold-change tends to decrease as concentration increases.

`FOLD` is fold vs vehicle (1 = no change, \< 1 = suppression). Percent
of positive control (`REL`) is not used. Isolated `FOLD < 1` points at
low concentration, which occur in ordinary induction curves, do not
trigger the flag.

The rank correlation uses order only, not a linear fit on the µM scale.
No p-value is required; typical induction experiments have too few
concentrations for a Spearman test to be informative.

## Examples

``` r
induction_downregulation(examplinib_in_vitro_ind)
#> # A tibble: 6 × 8
#>   OBJECT DONOR     n min_fold max_fold fold_maxc   rho downregulation
#>   <chr>  <chr> <int>    <dbl>    <dbl>     <dbl> <dbl> <lgl>         
#> 1 CYP1A2 A         4    1.38      7.37      7.37 1     FALSE         
#> 2 CYP1A2 B         4    1.24      5.09      5.09 1     FALSE         
#> 3 CYP1A2 C         4    1.56      4.15      4.15 1     FALSE         
#> 4 CYP3A4 A         7    0.904     6.33      6.33 1     FALSE         
#> 5 CYP3A4 B         8    0.655     2.31      2.31 0.619 FALSE         
#> 6 CYP3A4 C         8    0.648     1.58      1.58 0.238 FALSE         
induction_downregulation(examplinib_in_vitro_ind1)
#> # A tibble: 9 × 8
#>   OBJECT DONOR     n min_fold max_fold fold_maxc   rho downregulation
#>   <chr>  <chr> <int>    <dbl>    <dbl>     <dbl> <dbl> <lgl>         
#> 1 CYP1A2 A         4    0.087    0.951     0.087  -1   TRUE          
#> 2 CYP1A2 B         4    0.188    1.54      0.188  -0.8 TRUE          
#> 3 CYP1A2 C         4    0.072    1.09      0.072  -1   TRUE          
#> 4 CYP2B6 A         4    0.188    1.07      0.188  -0.8 TRUE          
#> 5 CYP2B6 B         4    0.341    1.36      0.341  -0.8 TRUE          
#> 6 CYP2B6 C         4    0.107    0.978     0.107  -1   TRUE          
#> 7 CYP3A4 A         4    0.285    0.973     0.285  -0.2 TRUE          
#> 8 CYP3A4 B         4    0.293    1.79      0.293  -0.8 TRUE          
#> 9 CYP3A4 C         4    0.286    1.06      0.286  -0.8 TRUE          
```
