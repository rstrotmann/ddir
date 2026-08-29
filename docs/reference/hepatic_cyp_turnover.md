# Hepatic CYP turnover data based on various publications

**\[deprecated\]**

Please use [`ddir::cyp_turnover`](cyp_turnover.md)instead.

## Usage

``` r
hepatic_cyp_turnover
```

## Format

A data frame with 6 columns:

- cyp: The CYP enzymne

- method: The experimental method used (below)

- mean_hl: The mean CYP degradation half-life measured,

- in_vivo: Study was conducted in vivo

- reference: The source publication (PMID or DOI)

- kdeg: The CYP degradation constant, i.e., \\log(2)/mean half-life\\.

## Source

This data set is taken from:

Yang J, Liao M, Shou M, Jamei M, Yeo KR, Tucker GT, Rostami-Hodjegan A.
Cytochrome p450 turnover: regulation of synthesis and degradation,
methods for determining rates, and implications for the prediction of
drug interactions. Curr Drug Metab. 2008 Jun;9(5):384-94. doi:
10.2174/138920008784746382. PMID: 18537575.

## Details

The following experimental methods were used in the original
publication:

- in vitro method 1: Radio-labeling of enzyme (‘pulse-chase’ method)

- in vitro method 2: Degradation of enzyme in cultured hepatocytes or
  liver slices

- in vivo method 1: Recovery of enzyme activity after enzyme induction

- in vivo method 2: Recovery of enzyme activity after mechanism-based
  inhibition (MBI)

- in vivo method 3: Pharmacokinetic modeling of auto-induction

These are the first few lines of the data frame:


       cyp            method mean_hl in_vivo                       reference   kdeg
    CYP1A2 In vitro Method 1      51   FALSE                   PMID: 2136526 0.0136
    CYP1A2 In vitro Method 2      43   FALSE   DOI: 10.1007/3-540-29804-5_25 0.0161
    CYP1A2 In vitro Method 2      36   FALSE                  PMID: 10997941 0.0193
    CYP1A2  In vivo Method 1      39    TRUE DOI: 10.1016/j.clpt.2004.04.003 0.0178
    CYP1A2  In vivo Method 3     105    TRUE    DOI: 10.1038/sj.clpt.6100431 0.0066
    CYP2A6 In vitro Method 2     226   FALSE                  PMID: 10997941 0.0031
