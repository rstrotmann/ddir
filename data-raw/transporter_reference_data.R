## code to prepare `transporter` dataset goes here
library(devtools)
library(usethis)

transporter_reference_data <- tribble(
  ~transporter, ~rank, ~fda_thld, ~ema_thld, ~i,
  "Pgp_int",    1,     10.0,      10.00,     "igut",
  "Pgp_sys",    2,     0.1,       0.02,      "imaxssu",
  "BCRP_int",   3,     10.0,      10.00,     "igut",
  "BCRP_sys",   4,     0.1,       0.02,      "imaxssu",
  "OCT1",       5,     NA,        0.04,      "imaxinletu",
  "OATP1B1",    6,     0.1,       0.04,      "imaxinletu",
  "OATP1B3",    7,     0.1,       0.04,      "imaxinletu",
  "OAT1",       8,     0.1,       0.04,      "imaxssu",
  "OAT3",       9,     0.1,       0.04,      "imaxssu",
  "BSEP",       10,    0.1,       0.02,      "imaxssu",
  "OCT2",       11,    0.1,       0.02,      "imaxssu",
  "MATE1",      12,    0.1,       0.02,      "imaxssu",
  "MATE2k",     13,    0.1,       0.02,      "imaxssu"
)

usethis::use_data(transporter_reference_data, overwrite = TRUE)


