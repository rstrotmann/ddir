# DDI inhibition data object constructor.

DDI inhibition data object constructor.

## Usage

``` r
inhibitor(data, object = "")
```

## Arguments

- data:

  Data frame with the following fields:

  ### Direct CYP or UGT inhibition

  - target: The CYP or UGT enzyme target as character.

  - ki: The ki in uM.

  - source: Source information (e.g., study name) as character.

  ### Time-dependent CYP inhibition

  - target: The CYP or UGT enzyme target as character.

  - ki: The ki in uM.

  - kinact: The kinact in 1/h.

  - source: Source information (e.g., study name) as character.

  ### Transporter inhibition

  - target: The transporter target as character.

  - ic50: The IC50 in uM.

  - source: Source information (e.g., study name) as character.

- object:

  Character

## Value

Inhibitor object
