#' Drug transporter reference data
#'
#' @source FDA and EMA guidelines.
"transporter_reference_data"



#' CYP reference substrate data
#'
#' CYP reference substrates commonly used in the mechanistic static assessment
#' of the CYP DDI perpetrator potential of drugs.
#'
#' @details
#' The CYP reference substrates currently implemented include:
#'
#'\preformatted{
#'       cyp   substrate fgut   fm fmcyp
#' 1  CYP1A2  tizanidine 1.00 0.95  0.98
#' 2  CYP2C8 repaglinide 1.00 1.00  0.61
#' 3  CYP2C9  S-warfarin 1.00 1.00  0.91
#' 4 CYP2C19  omeprazole 1.00 1.00  0.87
#' 5  CYP3A4   midazolam 0.57 0.96  1.00
#' }
#'
#' @source FDA and EMA guidelines.
#' @seealso [mech_stat_cyp_risk()]
#' @seealso [mech_stat_cyp_risk_table]
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
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_induction_data"


#' Fictional transporter inhibition data
#'
#' @source Fictional data, made up for demo purposes.
#' "examplinib_transporter_inhibition_data"


#' Fictional UGT inhibition data
#'
#' @format ## `examplinib_ugt_inhibition_data`
#' A data frame with the columns 'name', 'param', 'value' and 'source'
#'
#' @source Fictional data, made up for demo purposes.
#' "examplinib_ugt_inhibition_data"

