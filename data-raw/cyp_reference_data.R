## code to prepare `cyp_reference_data` dataset goes here
library(tidyverse)
library(usethis)

cyp_reference_substrates <- tibble::tribble(
  ~cyp,    ~substrate, ~fgut,  ~fm, ~fmcyp,
  "CYP1A2",  "tizanidine",     1, 0.95,   0.98,
  "CYP2C8", "repaglinide",     1,    1,   0.61,
  "CYP2C9",  "S-warfarin",     1,    1,   0.91,
  "CYP2C19",  "omeprazole",    1,    1,   0.87,
  "CYP3A4",   "midazolam",  0.57, 0.96,      1,
  "CYP2D6", "desipramine",     1,    1,   0.85
)

usethis::use_data(cyp_reference_substrates, overwrite = TRUE)


cyp_turnover <- tribble(
  ~cyp,      ~kdeg_hepatic, ~kdeg_intestinal,
  "CYP1A1",  0.0183, NA,
  "CYP1A2",  0.0183, NA,
  "CYP2A6",  0.0267, NA,
  "CYP2B6",  0.0217, NA,
  "CYP2C8",  0.0301, NA,
  "CYP2C9",  0.0067, 0.03,
  "CYP2C18", 0.0267, NA,
  "CYP2C19", 0.0267, 0.03,
  "CYP2D6",  0.0099, 0.03,
  "CYP2E1",  0.0176, NA,
  "CYP2J2",  0.0194, 0.03,
  "CYP3A4",  0.0193, 0.03,
  "CYP3A5",  0.0193, 0.03,
  "CYP3A7",  0.0193, NA,
)

usethis::use_data(cyp_turnover, overwrite = TRUE)
