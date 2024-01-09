#' Drug transporter reference data
#'
#' @source FDA and EMA guidelines.
"transporter_reference_data"



#' CYP reference substrate data
#'
#' @source FDA and EMA guidelines.
"cyp_reference_substrates"


#' Fictional substrate data
#'
#' A list of two: `perpetrator` objects:
#' * examplinib, a fictional drug.
#' * M1, a fictional metabolite of examplinib.
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_compounds"


#' Fictional examplnib parent compound data
#'
#' A `perpetrator` object for 'examplinib'.
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_parent"


#' Fictional examplinib metabolite compound data
#'
#' A `perpetrator` object for 'examplinib-M1'.
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_metabolite"


#' Fictional substrate data as string
#'
#' A text string containing data for two: `perpetrator` objects:
#' * examplinib, a fictional drug.
#' * M1, a fictional metabolite of examplinib.
#'
#' @source Fictional data, made up for demo purposes.
#' @seealso [read_perpetrators()]
"examplinib_compounds_string"


#' Fictional CYP inhibition data
#'
#' This data frame is a typical input to the following functions:
#'
#' [basic_cyp_inhibition_risk()]
#' [basic_cyp_inhibition_risk_table]
#' [mech_stat_cyp_risk()]
#' [mech_stat_cyp_risk_table()]
#'
#' @details
#' CYP inhibition data can contain ki data for multiple compounds.
#'
#' @format ## `examplinib_cyp_inhibition_data`
#' A data frame with the columns `name`, `param`, `value` and `source`, where:
#' * name is the name of the compound for which the data is recorded
#' * param contains the respective CYP enzyme names.
#' * value contains the Ki values for the respective CYP enzyme
#' * source provides information for the source of the respective value, often
#'     the name of the DMPK study. This entry is optional.
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_inhibition_data"


#' Fictional CYP induction data
#'
"examplinib_cyp_induction_data"


#' Fictional transporter inhibition data
#'
"examplinib_transporter_inhibition_data"

#' Fictional UGT inhibition data
#'
"examplinib_ugt_inhibition_data"

