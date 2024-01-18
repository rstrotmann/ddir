## code to prepare `examplinib` dataset goes here
library(tidyverse)
library(usethis)


## EXAMPLINIB COMPOUNDS

examplinib_compounds_string <- "
# PARENT\n\n# compound,  param,    value,     source
examplinib,  oral,     TRUE,
examplinib,  mw,       492.6,
examplinib,  dose,     450,       clinical dose
examplinib,  imaxss,   3530,      study 001
examplinib,  fu,       0.023,     study 002
examplinib,  fumic,    1,         default
examplinib,  rb,       1,         study 003
examplinib,  fa,       0.81,      study 003
examplinib,  fg,       1,         default
examplinib,  ka,       0.00267,   unknown

# METABOLITE
# compound,  param,    value,     source
M1,          oral,     FALSE,
M1,          mw,       506.56,
M1,          dose,     NA,
M1,          imaxss,   1038,      study 001
M1,          fu,       0.012,     study 002
M1,          fumic,    1,         default
M1,          rb,       1,         study 002
M1,          fa,       NA,
M1,          fg,       NA,
M1,          ka,
"

examplinib_compounds <- read_perpetrators(textConnection(
  examplinib_compounds_string))

examplinib_parent <- examplinib_compounds[[1]]
examplinib_metabolite <- examplinib_compounds[[2]]

usethis::use_data(examplinib_compounds_string, overwrite = TRUE)
usethis::use_data(examplinib_compounds, overwrite = TRUE)
usethis::use_data(examplinib_parent, overwrite = TRUE)
usethis::use_data(examplinib_metabolite, overwrite = TRUE)


## REVERSIBLE CYP INHIBITION

examplinib_cyp_inhibition_string <- "
# PARENT\n\n
# compound, CYP, ki, source\n\n
examplinib, CYP1A2,  NA,\n\n
examplinib, CYP2B6,  NA,\n\n
examplinib, CYP2C8,  11,   study 001\n\n
examplinib, CYP2C9,  13.5, study 001\n\n
examplinib, CYP2C19, 15,   study 001\n\n
examplinib, CYP2D6,  NA,\n\n
examplinib, CYP3A4,  12.5, study 001\n\n
\n\n
# METABOLITE\n\n
M1,         CYP2C9,  4.4,  study 002\n\n"

# examplinib_cyp_inhibition_data <- read_inhibitor_data(textConnection(
#   examplinib_cyp_inhibition_string))
examplinib_cyp_inhibition_data <- read_cyp_inhibitor_data(textConnection(
  examplinib_cyp_inhibition_string))

usethis::use_data(examplinib_cyp_inhibition_string, overwrite = TRUE)
usethis::use_data(examplinib_cyp_inhibition_data, overwrite = TRUE)


## TIME-DEPENDENT CYP INHIBITION

examplinib_cyp_tdi_string <- "
# compound, CYP,    ki,   kinact, source\n\n
examplinib, CYP3A4, 0.17, 0.04, study 001\n\n
"

examplinib_cyp_tdi_data <- read_tdi_data(textConnection(
  examplinib_cyp_tdi_string))

usethis::use_data(examplinib_cyp_tdi_string, overwrite = TRUE)
usethis::use_data(examplinib_cyp_tdi_data, overwrite = TRUE)


## CYP INDUCTION

examplinib_cyp_induction_string <- "
# PARENT\n\n
# compound, CYP, Emax, EC50, max c, source\n\n
examplinib, CYP1A2,  1,    NA,   5,  study 007\n\n
examplinib, CYP2B6,  1,    NA,   5,  study 007\n\n
examplinib, CYP2C8,  NA,   NA,   NA,\n\n
examplinib, CYP2C9,  NA,   NA,   NA,\n\n
examplinib, CYP2C19, NA,   NA,   NA,\n\n
examplinib, CYP2D6,  NA,   NA,   NA,\n\n
examplinib, CYP3A4,   7.35, 1.64, 3,  study 007\n\n
\n\n
# METABOLITE\n\n
# compound, CYP, ki, source\n\n
M1, CYP1A2,  1,    NA,   5,  study 007\n\n
M1, CYP2B6,  6.98, 1.86, 5,  study 007\n\n
M1, CYP2C8,  NA,   NA,   NA, \n\n
M1, CYP2C9,  NA,   NA,   NA, \n\n
M1, CYP2C19, NA,   NA,   NA, \n\n
M1, CYP2D6,  NA,   NA,   NA, \n\n
M1, CYP3A4,  22.7, 1.1,  5,  study 007\n\n"

examplinib_cyp_induction_data <- read_inducer_data(textConnection(
  examplinib_cyp_induction_string))

usethis::use_data(examplinib_cyp_induction_string, overwrite = TRUE)
usethis::use_data(examplinib_cyp_induction_data, overwrite = TRUE)

## UGT INHIBITION

examplinib_ugt_inhibition_string <- "
# PARENT\n\n
# compound, enzyme, IC50, source\n\n
examplinib, UGT1A1, 15, study 009\n\n
examplinib, UGT1A3, 15, study 009\n\n
examplinib, UGT1A4, 15, study 009\n\n
examplinib, UGT1A6, 15, study 009\n\n
examplinib, UGT1A9, 3.8, study 009\n\n
examplinib, UGT2B7, 15, study 009\n\n
examplinib, UGT2B15, 15, study 009\n\n
examplinib, UGT2B17, 6.1, study 009\n\n
\n\n
# METABOLITE\n\n
# compound, enzyme, IC50, source\n\n
M1, UGT1A1, 1.1, study 009\n\n
M1, UGT1A3, 5.8, study 009\n\n
M1, UGT1A4, 6.2, study 009\n\n
M1, UGT1A6, 15, study 009\n\n
M1, UGT1A9, 3.6, study 009\n\n
M1, UGT2B7, 15, study 009\n\n
M1, UGT2B15, 9.6, study 009\n\n"

# examplinib_ugt_inhibition_data <- read_inhibitor_data(textConnection(
#   examplinib_ugt_inhibition_string))

examplinib_ugt_inhibition_data <- read_ugt_inhibitor_data(textConnection(
  examplinib_ugt_inhibition_string))

usethis::use_data(examplinib_ugt_inhibition_string, overwrite = TRUE)
usethis::use_data(examplinib_ugt_inhibition_data, overwrite = TRUE)

## TRANSPORTER INHIBITION

examplinib_transporter_inhibition_string <- "
# PARENT\n\n
# name,     transporter, IC50, source\n\n
examplinib, Pgp,       0.41,  study 005\n\n
examplinib, BCRP,      1.9,  study 005\n\n
examplinib, OCT1,      2.3,    study 006\n\n
examplinib, OATP1B1,   177,   study 006\n\n
examplinib, OATP1B3,   35,   study 006\n\n
examplinib, OAT1,      271,\n\n
examplinib, OAT3,      300,  \n\n
examplinib, BSEP,      12.8,\n\n
examplinib, OCT2,      67,    study 006\n\n
examplinib, MATE1,     3.6,    study 006\n\n
examplinib, MATE2k,    1.1,    study 006\n\n"

examplinib_transporter_inhibition_data <- read_transporter_inhibitor_data(
  textConnection(examplinib_transporter_inhibition_string))

usethis::use_data(examplinib_transporter_inhibition_string, overwrite = TRUE)
usethis::use_data(examplinib_transporter_inhibition_data, overwrite = TRUE)









