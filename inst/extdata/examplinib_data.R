compounds <-read_perpetrators(textConnection("
  # PARENT
  # compound,  param,    value,     source
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
  M1,          ka,       NA,
"))


cyp_inhibition_data <- read_cyp_inhibitor_data(textConnection("
  # PARENT
  # compound, CYP, ki, source
  examplinib, CYP1A2,  NA
  examplinib, CYP2B6,  NA,
  examplinib, CYP2C8,  11,   study 001
  examplinib, CYP2C9,  13.5, study 001
  examplinib, CYP2C19, 15,   study 001
  examplinib, CYP2D6,  NA,
  examplinib, CYP3A4,  12.5, study 001

  # METABOLITE
  M1,         CYP2C9,  4.4,  study 002
"))


cyp_tdi_data <- read_tdi_data(textConnection("
  # compound, CYP, Ki, Kinact, source
  examplinib, CYP3A4, 0.17,   0.04, study 001
"))


cyp_induction_data <- read_inducer_data(textConnection("
  # PARENT
  # compound, CYP, Emax, EC50, max c, source
  examplinib, CYP1A2,  1,    NA,   5,  study 007
  examplinib, CYP2B6,  1,    NA,   5,  study 007
  examplinib, CYP2C8,  NA,   NA,   NA,
  examplinib, CYP2C9,  NA,   NA,   NA,
  examplinib, CYP2C19, NA,   NA,   NA,
  examplinib, CYP2D6,  NA,   NA,   NA,
  examplinib, CYP3A4,   7.35, 1.64, 3,  study 007

  # METABOLITE
  # compound, CYP, ki, source
  M1, CYP1A2,  1,    NA,   5,  study 007
  M1, CYP2B6,  6.98, 1.86, 5,  study 007
  M1, CYP2C8,  NA,   NA,   NA,
  M1, CYP2C9,  NA,   NA,   NA,
  M1, CYP2C19, NA,   NA,   NA,
  M1, CYP2D6,  NA,   NA,   NA,
  M1, CYP3A4,  22.7, 1.1,  5,  study 007
"))


ugt_inhibition_data <- read_ugt_inhibitor_data(textConnection("
# PARENT
# compound, enzyme, IC50, source
examplinib, UGT1A1, 15, study 009
examplinib, UGT1A3, 15, study 009
examplinib, UGT1A4, 15, study 009
examplinib, UGT1A6, 15, study 009
examplinib, UGT1A9, 3.8, study 009
examplinib, UGT2B7, 15, study 009
examplinib, UGT2B15, 15, study 009
examplinib, UGT2B17, 6.1, study 009

# METABOLITE
# compound, enzyme, IC50, source
M1, UGT1A1, 1.1, study 009
M1, UGT1A3, 5.8, study 009
M1, UGT1A4, 6.2, study 009
M1, UGT1A6, 15, study 009
M1, UGT1A9, 3.6, study 009
M1, UGT2B7, 15, study 009
M1, UGT2B15, 9.6, study 009
"))


transporter_inhibition_data <- read_transporter_inhibitor_data(textConnection("
# PARENT
# name,     transporter, IC50, source
examplinib, Pgp,       0.41,  study 005
examplinib, BCRP,      1.9,  study 005
examplinib, OCT1,      2.3,    study 006
examplinib, OATP1B1,   177,   study 006
examplinib, OATP1B3,   35,   study 006
examplinib, OAT1,      271,
examplinib, OAT3,      300,
examplinib, BSEP,      12.8,
examplinib, OCT2,      67,    study 006
examplinib, MATE1,     3.6,    study 006
examplinib, MATE2k,    1.1,    study 006
"))
