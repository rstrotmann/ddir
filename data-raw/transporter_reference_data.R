## code to prepare `transporter` dataset goes here
library(devtools)
library(usethis)

transporter_reference_data <- tibble::tribble(
  ~transporter, ~threshold,           ~i,
     "Pgp_int",         10,       "igut",
     "Pgp_sys",       0.02,    "imaxssu",
    "BCRP_int",         10,       "igut",
    "BCRP_sys",       0.02,    "imaxssu",
     "OATP1B1",        0.1, "imaxinletu",
     "OATP1B3",        0.1, "imaxinletu",
        "OAT1",        0.1,    "imaxssu",
        "OAT3",        0.1,    "imaxssu",
        "BSEP",        0.1,    "imaxssu",
        "OCT1",         NA,           NA,
        "OCT2",        0.1,    "imaxssu",
       "MATE1",       0.02,    "imaxssu",
      "MATE2k",       0.02,    "imaxssu"
  )


usethis::use_data(transporter_reference_data, overwrite = TRUE)


