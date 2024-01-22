pull_data_string <- function(filename) {
  paste(readLines(
    fs::path_package(paste0("inst/extdata/", filename),
                     package = "ddir")),
    collapse = "\n")
}


## EXAMPLINIB COMPOUNDS

examplinib_compounds_string <- pull_data_string("examplinib_compounds.csv")
examplinib_compounds_single_string <- pull_data_string("examplinib_compounds_single.csv")

examplinib_compounds <- read_perpetrators(textConnection(
  examplinib_compounds_string))

examplinib_parent <- examplinib_compounds[[1]]
examplinib_metabolite <- examplinib_compounds[[2]]

usethis::use_data(examplinib_compounds_string, overwrite = TRUE)
usethis::use_data(examplinib_compounds_single_string, overwrite = TRUE)
usethis::use_data(examplinib_compounds, overwrite = TRUE)
usethis::use_data(examplinib_parent, overwrite = TRUE)
usethis::use_data(examplinib_metabolite, overwrite = TRUE)


## REVERSIBLE CYP INHIBITION

examplinib_cyp_inhibition_string <- pull_data_string("examplinib_cyp_inhibition.csv")

examplinib_cyp_inhibition_data <-
  read_cyp_inhibitor_data(textConnection(examplinib_cyp_inhibition_string))

usethis::use_data(examplinib_cyp_inhibition_string, overwrite = TRUE)
usethis::use_data(examplinib_cyp_inhibition_data, overwrite = TRUE)


## TIME-DEPENDENT CYP INHIBITION

examplinib_cyp_tdi_string <- pull_data_string("examplinib_cyp_tdi.csv")

examplinib_cyp_tdi_data <- read_tdi_data(textConnection(
  examplinib_cyp_tdi_string))

usethis::use_data(examplinib_cyp_tdi_string, overwrite = TRUE)
usethis::use_data(examplinib_cyp_tdi_data, overwrite = TRUE)


## CYP INDUCTION

examplinib_cyp_induction_string <- pull_data_string("examplinib_cyp_induction.csv")

examplinib_cyp_induction_data <- read_inducer_data(textConnection(
  examplinib_cyp_induction_string))

usethis::use_data(examplinib_cyp_induction_string, overwrite = TRUE)
usethis::use_data(examplinib_cyp_induction_data, overwrite = TRUE)


## UGT INHIBITION

examplinib_ugt_inhibition_string <- pull_data_string("examplinib_ugt_inhibition.csv")

examplinib_ugt_inhibition_data <-
  read_ugt_inhibitor_data(textConnection(examplinib_ugt_inhibition_string))

usethis::use_data(examplinib_ugt_inhibition_string, overwrite = TRUE)
usethis::use_data(examplinib_ugt_inhibition_data, overwrite = TRUE)


## TRANSPORTER INHIBITION

examplinib_transporter_inhibition_string <-
  pull_data_string("examplinib_transporter_inhibition.csv")

examplinib_transporter_inhibition_data <-
  read_transporter_inhibitor_data(textConnection(
    examplinib_transporter_inhibition_string))

usethis::use_data(examplinib_transporter_inhibition_string, overwrite = TRUE)
usethis::use_data(examplinib_transporter_inhibition_data, overwrite = TRUE)
