# Drug transporter reference data

Drug transporter reference data

## Usage

``` r
transporter_reference_data
```

## Format

- 'transporter' The name of the drug transporter protein.

- 'threshold' The regulatory threshold for clinically relevant
  interactions.

- 'i' The precipitant concentration metric applicable for the
  interaction.

## Source

Section 2 of the [ICH M12
guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).

## Details


    transporter , threshold , i
    Pgp_int     , 10        , igut
    Pgp_sys     , 0.02      , imaxssu
    BCRP_int    , 10        , igut
    BCRP_sys    , 0.02      , imaxssu
    OATP1B1     , 0.1       , imaxinletu
    OATP1B3     , 0.1       , imaxinletu
    OAT1        , 0.1       , imaxssu
    OAT3        , 0.1       , imaxssu
    BSEP        , 0.1       , imaxssu
    OCT2        , 0.1       , imaxssu
    MATE1       , 0.02      , imaxssu
    MATE2k      , 0.02      , imaxssu

## See also

[`key_conc_table()`](key_conc_table.md)

[`transporter_inhibition_risk()`](transporter_inhibition_risk.md)
