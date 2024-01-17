## code to prepare `transporter` dataset goes here
library(devtools)
library(usethis)

transporter_reference_data <- tribble(
  ~transporter, ~fda_thld, ~ema_thld, ~i,
  "Pgp_int",    10.0,      10.00,     "igut",
  "Pgp_sys",    0.1,       0.02,      "imaxssu",
  "BCRP_int",   10.0,      10.00,     "igut",
  "BCRP_sys",   0.1,       0.02,      "imaxssu",
  "OCT1",       NA,        0.04,      "imaxinletu",
  "OATP1B1",    0.1,       0.04,      "imaxinletu",
  "OATP1B3",    0.1,       0.04,      "imaxinletu",
  "OAT1",       0.1,       0.04,      "imaxssu",
  "OAT3",       0.1,       0.04,      "imaxssu",
  "BSEP",       0.1,       0.02,      "imaxssu",
  "OCT2",       0.1,       0.02,      "imaxssu",
  "MATE1",      0.1,       0.02,      "imaxssu",
  "MATE2k",     0.1,       0.02,      "imaxssu"
)

usethis::use_data(transporter_reference_data, overwrite = TRUE)


