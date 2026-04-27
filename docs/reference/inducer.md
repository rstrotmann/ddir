# DDI CYP inducer object constructor

This class constructor expects as the as the data argument a data frome
with the following columns:

## Usage

``` r
inducer(data, object = "")
```

## Arguments

- data:

  Data frame with the columns 'target', 'emax', 'ec50', 'max_c' and
  'source'.

- object:

  The name of the DDI perpetrator as character.

## Value

Inducer object

## Details

- target: The CYP induction target. Must be one of "CYP1A1", "CYP1A2",
  "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C18", "CYP2C19",
  "CYP2D6", "CYP2E1", "CYP2J2", "CYP3A4", "CYP3A5", or "CYP3A7".

- emax: The maximal induction effect.

- ec50: The EC50, i.e., the concentration causing the half-maximal
  induction effect.

- max_c: The maximal concentration tested in the assay in uM.

- source: Source information for the parameter.
