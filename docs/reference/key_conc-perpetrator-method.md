# Key perpetrator concentrations

This function prints a markdown-formatted table of the relevant
perpetrator concentrations in \\\mu M\\ and ng/ml for a DDI perpetrator
compound.

## Usage

``` r
# S4 method for class 'perpetrator'
key_conc(x, round = 2, qh = 1.616, qent = 18/60)
```

## Arguments

- x:

  DDI perpetrator object.

- round:

  Number of decimal places.

- qh:

  Hepatic blood flow in L/min.

- qent:

  Enteric blood flow in L/min.

## Value

Markdown-formatted text.

## Details

### Gut concentration

\$\$I\_{gut} = \frac{D} {250}\$\$

If the above exceeds the aqueous solubility of the drug, \\I\_{gut}\\ is
set to its solubility.

### Unbound systemic concentration

\$\$I\_{max,ss,u}=I\_{max,ss} \* f_u\$\$

### Unbound hepatic inlet concentration

For orally administered (parent) compounds, the hepatic inlet
concentration is the systemic concentration plus a portal term:

\$\$portal\\ term = D\*\frac{F_a\*F_g\*k_a}{Q_h\*R_B}\*1000\\ ng/ml\$\$

where \\D\\ is the administered dose in mg, \\F_a\\ the fraction
absorbed after oral administration, \\F_g\\ the fraction available after
gut metabolism, \\k_a\\ the absorption rate, \\Q_h\\ the hepatic blood
flow and \\R_B\\ the blood-to-plasma ratio.

The relevant hepatic inlet concentration (\\I\_{max,inlet,u}\\, also
called \\I_h\\ in the mechanistic static modeling equations)
concentration is the sum of the maximal systemic plasma concentration
and the portal contribution:

\$\$I\_{max,inlet,u}=(I\_{max,ss} + portal\\ term) \* f_u\$\$

### Enteric (villous) concentration

For orally administered (parent) compounds, the villous concentration in
the gut (\\I\_{enteric}\\, also called \\I_g\\ in the mechanistic static
modeling equations) is calculated as:

\$\$I\_{enteric,u} = D \* \frac{F_a\*k_a}{Q\_{ent}} \*1000\\ ng/ml\$\$

where \\F_a\\ is the fraction absorbed after oral administration,
\\k_a\\ the absorption rate, \\Q\_{ent}\\ the enteric villous blood flow
and \\R_B\\ the blood-to-plasma distribution ratio of the compound.

Note that as per the FDA guideline (refer to FDA, 2020, Fig. 7, and
Rostami-Hodjegan and Tucker, 2004) the blood-to-plasma ratio and the
plasma binding of the drug are not applicable.
