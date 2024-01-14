## code to prepare `cyp_reference_data` dataset goes here
library(tidyverse)
library(usethis)

cyp_reference_substrates <- tribble(
  ~cyp,      ~substrate,    ~fgut, ~fm,  ~fmcyp,
  "CYP1A2",  "tizanidine",  1.00,  0.95, 0.98,
  "CYP2C8",  "repaglinide", 1.00,  1.00, 0.61,
  "CYP2C9",  "S-warfarin", 1.00,  1.00, 0.91,
  "CYP2C19", "omeprazole",  1.00,  1.00, 0.87,
  "CYP3A4",  "midazolam",   0.57,  0.96, 1.00
)

usethis::use_data(cyp_reference_substrates, overwrite = TRUE)
