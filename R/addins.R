#' Insert a precipitant object constructor
#'
#' @export
precipitant_addin <- function() {
  rstudioapi::insertText(
    r"(precipitant(
  name = "name",
  dose = 100,
  imaxss = 1000,
  mw = 500,
  oral = TRUE,
  fu = 1,
  fumic = 1,
  rb = 1,
  fa = 1,
  fg = 1,
  ka = 0.1,
  solubility = Inf
))"
  )
}


#' Insert a CYP inhibitor object constructor
#'
#' @export
cyp_inh_addin <- function() {
  rstudioapi::insertText(
    r"(inhibition_data(
  data = tibble::tribble(
    ~object, ~ki, ~source,
   "CYP1A2",  NA,      "",
   "CYP2B6",  NA,      "",
   "CYP2C8",  NA,      "",
   "CYP2C9",  NA,      "",
  "CYP2C19",  NA,      "",
   "CYP2D6",  NA,      "",
   "CYP3A4",  NA,      ""
  ),
  precipitant = "name"
))"
  )
}


#' Insert a CYP inhibitor object constructor
#'
#' @export
ugt_inh_addin <- function() {
  rstudioapi::insertText(
    r"(inhibition_data(
  data = tibble::tribble(
     ~object, ~ki, ~source,
    "UGT1A1",  NA,      "",
    "UGT1A3",  NA,      "",
    "UGT1A4",  NA,      "",
    "UGT1A6",  NA,      "",
    "UGT1A9",  NA,      "",
    "UGT2B7",  NA,      "",
   "UGT2B15",  NA,      "",
   "UGT2B17",  NA,      ""
   ),
  precipitant = "name"
))"
  )
}


#' Insert a transporter inhibitor object constructor
#'
#' @export
transp_inh_addin <- function() {
  rstudioapi::insertText(
    r"(inhibition_data(
    data = tibble::tribble(
    ~object, ~ki, ~source,
      "Pgp",  NA,      "",
     "BCRP",  NA,      "",
  "OATP1B1",  NA,      "",
  "OATP1B3",  NA,      "",
     "OAT1",  NA,      "",
     "OAT3",  NA,      "",
     "BSEP",  NA,      "",
     "OCT1",  NA,      "",
     "OCT2",  NA,      "",
    "MATE1",  NA,      "",
   "MATE2k",  NA,      ""
  ),
  precipitant = "name"
))"
  )
}


#' Insert a CYP inducer object constructor
#'
#' @export
cyp_ind_addin <- function() {
  rstudioapi::insertText(
    r"(induction_data(
  data = tibble::tribble(
     ~object, ~emax, ~ec50, ~max_c, ~source,
    "CYP3A4",    NA,    NA,     NA,      "",
    "CYP2B6",    NA,    NA,     NA,      "",
    "CYP1A2",    NA,    NA,     NA,      "",
    "CYP2C8",    NA,    NA,     NA,      "",
    "CYP2C9",    NA,    NA,     NA,      "",
   "CYP2C19",    NA,    NA,     NA,      ""
  ),
  precipitant = "name"
))"
  )
}


#' Insert the standard DDI assessment report appendix as mardown
#'
#' @export
appendix_addin <- function() {
  rstudioapi::insertText(
    r"(\newpage
# Appendix 1: Calculations and formulas

## Relevant precipitant drug concentrations

### Gut concentration {-}

The maximal gut concentration ($I_{gut}$) for the
orally administered compounds is the administered dose dissolved in 250 ml.

$$I_{gut} = \frac{D} {250}$$

### Systemic concentration {-}

The unbound systemic ($C_{max,ss,u}$) concentration is considered the relevant
precipitant concentration for hepatic enzyme inhibition and induction:

$$C_{max,ss,u}=I_{max,ss} * f_u$$

### Hepatic inlet concentration {-}

The hepatic inlet concentration is considered the relevant precipitant
concentration for inhibition of the hepatic uptake transporters OATP1B1 and
OATP1B3, and for the hepatic terms in the mechanistic static modeling equation
(refer to Section '[Mechanistic static modeling of CYP inhibition/induction]').

The hepatic inlet concentration is composed of the systemic concentration and
the portal contribution. For orally administered drugs, the portal term is
calculated as:

$$portal\ term = D * \frac{F_a * F_g * k_a}{Q_h * R_B} * 1000\ ng/ml$$

with

* $D$ the administered dose in mg
* $F_a$ the fraction absorbed after oral administration
* $F_g$ the fraction available after gut metabolism
* $k_a$ the absorption rate
* $Q_h$ the hepatic blood flow
* $R_B$ the blood-to-plasma ratio.

The standard hepatic blood flow is assumed as 97 l/h/70 kg or 1.61 l/min/70 kg.

The relevant hepatic inlet ($I_{max,inlet,u}$, also called $I_h$ in the
mechanistic static modeling equations) concentration is the sum of the maximal
systemic plasma concentration and the portal contribution:

$$I_{max,inlet,u}=(C_{max,ss} + portal\ term) * f_u$$

### Enteric concentration {-}

For the parent compound, the villous concentration in the gut ($I_{enteric}$,
also called $I_g$ in the mechanistic static modeling equations) is calculated
as:

$$I_{enteric,u} = D * \frac{F_a*k_a}{Q_{ent}} *1000\ ng/ml$$

with

* $F_a$ the fraction absorbed after oral administration
* $k_a$ the absorption rate constant
* $Q_{ent}$ the enteric villous blood flow

Note that as per the ICH M12 guideline and [Rostami-Hodjegan and Tucker,
2004](https://doi.org/10.1016/j.ddtec.2004.10.002) the blood-to-plasma
distribution ratio and the plasma binding of the drug are not applicable for the
calculation of the villous concentration.

The standard villous blood flow is assumed as 18 l/h/70 kg or 0.3 l/min/70 kg.

## Basic modeling of enzyme inhibition

### Reversible inhibition

For the basic modeling of direct (reversible) enzyme inhibition, the ratios of
the relevant inhibitor concentration to the $K_{i,u}$ are considered (refer to
Section 2.1.2.1 of the [ICH M12 guidance document](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf)).

For in vitro studies conducted using human liver microsomes, the microsomal
unbound fraction, $f_{u,mic}$ is used to calculate $K_{i,u}$. If unknown, a
default of 1 is assumed.

$R$ values larger than 0.02 (liver) or 10 (gut), are considered to indicate a
potential clinical enzyme inhibition risk using this method.

#### Liver {-}

$$R=\frac{C_{max,ss,u}}{K_{i,u}}$$

#### Gut wall {-}

$$R_{gut}=\frac{I_{gut}}{K_{i,u}}$$

### Time-dependent CYP inhibition

For the basic modeling of the potential for time-dependent CYP inhibition (TDI),
the following metric is considered:

$$R=\frac {k_{obs} + k_{deg}}{k_{deg}}$$

with

$$k_{obs}=\frac {5 * k_{inact}*C_{max,u}}{K_{I,u} + 5 * C_{max,u}}$$

The CYP degradation constant, $k_{deg}$ is a physiological constant that should
be derived from the scientific literature. In this DDI assessment report,
standard values are used unless otherwise indicated.

Values of $R \ge 1.25$ is considered to indicate a clinically relevant TDI
potential and suggest the need for further investigation.

## Basic mRNA fold-change method to assess CYP induction

This basic risk assessment evaluates the mRNA induction for a set of hepatocyte
batches from different donors. Increases of CYP enzyme mRNA $\ge$ 2-fold at
at concentrations up to 50-fold above $C_{max,u}$ is considered to indicate a
clinical risk for CYP induction.

In the context of this assessment only the worst-case donor data is considered.

## Basic kinetic modeling of CYP induction

For the basic kinetic modeling of the CYP induction potential, the following
metric is considered (refer to Section 2.1.4.3 of the ICH M12 guideline):

$$R = \frac {1}{1+d* \frac {E_{max}*10*C_{max,ss,u}}{EC_{50,u} + 10*C_{max,ss,u}}}$$

with $d$ a scaling factor that has a standard value of 1. A different value can
be used if warranted by prior experience with the experimental conditions.

$R \le 0.8$ suggest a relevant in vivo CYP induction potential.


## Mechanistic static modeling of CYP modulation

In this approach, AUC ratios for specific DDI object substrates are projected
based on their known intestinal and hepatic metabolism. Both direct
(competitive) and time-dependent inhibition, as well as enzyme induction are
considered. AUC ratios are calculated according to the below formula (refer to
Section 7.5.1.2 of the ICH M12 guideline):

$$AUCR = \frac{1}{A_g*B_g*C_g*\left(1-F_g\right)+F_g} * \frac{1}{A_h*B_h*C_h*f_m+\left(1-f_m\right)}$$

This calculation is applied for typical probe substrates for which $F_g$, i.e.,
the fraction escaping gut metabolism and $f_m$, i.e., the fraction metabolized
are known.

Note that the $f_m$ is composed of the overall fraction metabolized for the
respective probe substrate, and the fraction metabolized by the CYP enzyme in
questions:

$$f_m=f_{m,overall} * f_{m,CYP}$$

The individual terms in the AUC calculation are:

### Reversible inhibition {-}

$$A_g = \frac{1}{1+\frac{I_g}{K_i}}$$

$$A_h = \frac{1}{1+\frac{I_h}{K_i}}$$

### Time-dependent inhibition {-}

$$B_g = \frac{k_{deg,g}}{k_{deg,g} + \frac{I_g*k_{inact}}{I_g+K_I}}$$

$$B_h = \frac{k_{deg,h}}{k_{deg,h} + \frac{I_h*k_{inact}}{I_h+K_I}}$$

### Induction {-}

$$C_g = 1 + \frac{d*E_{max}*I_g}{I_g+EC_{50}}$$

$$C_h = 1 + \frac{d*E_{max}*I_h}{I_h+EC_{50}}$$

with the hepatic inlet concentration $I_h=I_{max,inlet,u}$ and the intestinal
concentration $I_g=I_{enteric}$ (see above). $d$ is an induction scaling factor
(assumed to be 1 but can be adjusted based on the experimental conditions).

If the predicted AUC ratio is outside of the 0.8 to 1.25 interval, further
evaluation is required.

## Inhibition of drug transporters

As per the M12 guideline, the metric for the assessment of the drug transporter
inhibition risk is:

$$R=[I]/IC_{50,u}$$

In the respective in vitro studies, the substrate concentration is usually very
low, so that $K_i \approx IC_{50}$ can be assumed. Under common assay conditions,
no protein is added to the medium so that the fraction unbound can be assumed 1,
i.e. $IC_{50} = IC_{50,u}$.

The following relevant precipitant concentrations $[I]$ and regulatory
thresholds of concern apply for the transporters:

| $I$                     | transporter                                                                    | threshold |
| ---                   | ---                                                                            | ---       |
| $I_{gut}$         | P-gp and BRCR when drugs are orally administered                               | 10        |
| $C_{max,ss,u}$    | P-gp and BRCR when drugs are administered parenterally or for drug metabolites | 0.02      |
| $I_{max,inlet,u}$ | hepatic basolateral transporters OCT1, OATP1B1 and OATP1B3                     | 0.1       |
| $C_{max,ss,u}$    | renal basolateral transporters OAT1, OAT3 and OCT2                             | 0.1       |
| $C_{max,ss,u}$    | apical transporters MATE1 and MATE2-K                                          | 0.02      |


Refer to Section '[Relevant precipitant drug concentrations]' for the calculation
of the relevant precipitant concentrations.
    )"
  )
}





