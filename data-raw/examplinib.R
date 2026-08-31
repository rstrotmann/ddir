examplinib_raw <- tibble::tribble(
    ~param,  ~value,         ~source,
    "oral",       1,              NA,
      "mw",   492.6,              NA,
    "dose",     450, "clinical dose",
  "imaxss",    3530,     "study 001",
      "fu",   0.023,     "study 002",
   "fumic",       1,       "default",
      "rb",       1,     "study 003",
      "fa",    0.81,     "study 003",
      "fg",       1,              NA,
      "ka", 0.00267,       "unknown"
  )

examplinib <- do.call(
  precipitant,
  append(
    setNames(as.list(examplinib_raw$value), examplinib_raw$param),
    list(name = "examplinib")
  )
)
usethis::use_data(examplinib, overwrite = TRUE)


examplinib_cyp_inhibition <- tibble::tribble(
    ~object,  ~ki,     ~source,
   "CYP1A2",   NA,          NA,
   "CYP2B6",   NA,          NA,
   "CYP2C8",   11, "study 001",
   "CYP2C9",  0.6, "study 001",
  "CYP2C19", 0.25, "study 001",
   "CYP2D6",   NA,          NA,
   "CYP3A4", 12.5, "study 001"
  ) |>
  inhibition_data(precipitant = "examplinib")
usethis::use_data(examplinib_cyp_inhibition, overwrite = TRUE)


examplinib_cyp_tdi <- tibble::tribble(
   ~object,  ~ki, ~kinact,     ~source,
  "CYP3A4", 0.17,    0.04, "study 001"
  )|>
  inhibition_data(precipitant = "examplinib")
usethis::use_data(examplinib_cyp_tdi, overwrite = TRUE)


examplinib_cyp_induction <- tibble::tribble(
    ~object, ~emax, ~ec50, ~max_c,     ~source,
   "CYP1A2",     1,    NA,      5, "study 007",
   "CYP2B6",     1,    NA,      5, "study 007",
   "CYP2C8",    NA,    NA,     NA,          NA,
   "CYP2C9",    NA,    NA,     NA,          NA,
  "CYP2C19",    NA,    NA,     NA,          NA,
   "CYP2D6",    NA,    NA,     NA,          NA,
   "CYP3A4",  7.35,  1.64,      3, "study 007"
  ) |>
  induction_data(precipitant = "examplinib")
usethis::use_data(examplinib_cyp_induction, overwrite = TRUE)

examplinib_ugt_inhibition <-tibble::tribble(
    ~object, ~ki,     ~source,
   "UGT1A1",  15, "study 009",
   "UGT1A3",  15, "study 009",
   "UGT1A4",  15, "study 009",
   "UGT1A6",  15, "study 009",
   "UGT1A9", 3.8, "study 009",
   "UGT2B7",  15, "study 009",
  "UGT2B15",  15, "study 009",
  "UGT2B17", 6.1, "study 009"
  ) |>
  inhibition_data(precipitant = "examplinib")
usethis::use_data(examplinib_ugt_inhibition, overwrite = TRUE)

examplinib_transporter_inhibition <-tibble::tribble(
    ~object,  ~ic50,     ~source,
      "Pgp", 0.41, "study 005",
     "BCRP",  1.9, "study 005",
     "OCT1",  2.3, "study 006",
  "OATP1B1",  177, "study 006",
  "OATP1B3",   35, "study 006",
     "OAT1",  271,          NA,
     "OAT3",  300,          NA,
     "BSEP", 12.8,          NA,
     "OCT2",   67, "study 006",
    "MATE1",  3.6, "study 006",
   "MATE2k",  1.1, "study 006"
  ) |>
  inhibition_data(precipitant = "examplinib")
usethis::use_data(examplinib_transporter_inhibition, overwrite = TRUE)

