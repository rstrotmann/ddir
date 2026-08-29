# CYP reference substrate data

CYP reference substrates commonly used in the mechanistic static
assessment of the CYP DDI potential of drugs.

## Usage

``` r
cyp_reference_substrates
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 6
rows and 5 columns.

## Source

FDA and EMA guidelines.

## Details

The CYP reference substrates currently implemented include:


    cyp     , substrate   , fgut , fm   , fmcyp
    CYP1A2  , tizanidine  , 1    , 0.95 , 0.98
    CYP2C8  , repaglinide , 1    , 1    , 0.61
    CYP2C9  , S-warfarin  , 1    , 1    , 0.91
    CYP2C19 , omeprazole  , 1    , 1    , 0.87
    CYP3A4  , midazolam   , 0.57 , 0.96 , 1

## See also

[`mech_stat_cyp_risk()`](mech_stat_cyp_risk.md)
