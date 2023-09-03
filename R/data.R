#' Drug transporter reference data
#'
#' @format ## `transporter_reference_data`
#' @source FDA and EMA guidelines.
"transporter_reference_data"



#' CYP reference substrate data
#'
#' @format ## `cyp_reference_substrates`
#' @source FDA and EMA guidelines.
"cyp_reference_substrates"


#' Fictional substrate data
#'
#' A list of two: `perpetrator`objects:
#' * examplinib, a fictional drug.
#' * M1, a fictional metabolite of examplinib.
#'
#' @source Fictional data, made up for demo purposes.
#' @seealso [perpetrator()]
#' @seealso [load_perpetrators()]
"examplinib_compounds"


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


