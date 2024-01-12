#' Drug transporter reference data
#'
#' @details
#' \preformatted{
#' #'    param rank fda_thld ema_thld          i
#'  Pgp_int    1     10.0    10.00       igut
#'  Pgp_sys    2      0.1     0.02    imaxssu
#' BCRP_int    3     10.0    10.00       igut
#' BCRP_sys    4      0.1     0.02    imaxssu
#'     OCT1    5       NA     0.04 imaxinletu
#'  OATP1B1    6      0.1     0.04 imaxinletu
#'  OATP1B3    7      0.1     0.04 imaxinletu
#'     OAT1    8      0.1     0.04    imaxssu
#'     OAT3    9      0.1     0.04    imaxssu
#'     BSEP   10      0.1     0.02    imaxssu
#'     OCT2   11      0.1     0.02    imaxssu
#'    MATE1   12      0.1     0.02    imaxssu
#'   MATE2k   13      0.1     0.02    imaxssu
#' }
#' @source FDA and EMA guidelines.
"transporter_reference_data"



#' CYP reference substrate data
#'
#' CYP reference substrates commonly used in the mechanistic static assessment
#' of the CYP DDI perpetrator potential of drugs.
#'
#' @details
#' The CYP reference substrates currently implemented include:
#' \preformatted{
#'     cyp   substrate fgut   fm fmcyp
#'  CYP1A2  tizanidine 1.00 0.95  0.98
#'  CYP2C8 repaglinide 1.00 1.00  0.61
#'  CYP2C9  S-warfarin 1.00 1.00  0.91
#' CYP2C19  omeprazole 1.00 1.00  0.87
#'  CYP3A4   midazolam 0.57 0.96  1.00
#' }
#' @source FDA and EMA guidelines.
#' @seealso [mech_stat_cyp_risk()]
#' @seealso [mech_stat_cyp_risk_table]
"cyp_reference_substrates"


#' Perpetrator compound data for examplinib and its metabolite
#'
#' A list of two: `perpetrator` objects:
#' * examplinib, a fictional drug.
#' * M1, a fictional metabolite of examplinib.
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_compounds"


#' Perpetrator compound data for examplinib
#'
#' @details
#' \preformatted{
#'                 param      value        source
#' name             name examplinib
#' oral             oral       TRUE
#' mw                 mw      492.6
#' dose             dose        450 clinical dose
#' imaxss         imaxss       3530     study 001
#' fu                 fu      0.023     study 002
#' fumic           fumic          1       default
#' rb                 rb          1     study 003
#' fa                 fa       0.81     study 003
#' fg                 fg          1       default
#' ka                 ka    0.00267       unknown
#' solubility solubility        Inf       default
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_parent"


#' Perpetrator compound data for examplinib
#'
#' @details
#' \preformatted{
#'                 param  value    source
#' name             name     M1
#' oral             oral  FALSE
#' mw                 mw 506.56
#' dose             dose     NA
#' imaxss         imaxss   1038 study 001
#' fu                 fu  0.012 study 002
#' fumic           fumic      1   default
#' rb                 rb      1 study 002
#' fa                 fa     NA
#' fg                 fg     NA
#' ka                 ka     NA
#' solubility solubility    Inf   default
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_metabolite"


#' Examplinib compound data as string
#'
#' A character string containing data for two perpetrator objects:
#' * examplinib, a fictional drug.
#' * M1, a fictional metabolite of examplinib.
#'
#' @source Fictional data, made up for demo purposes.
#' @seealso [read_perpetrators()]
"examplinib_compounds_string"


#' Examplinib CYP inhibition data
#'
#' This data frame is a typical input to the following functions:
#' [basic_cyp_inhibition_risk()]
#' [basic_cyp_inhibition_risk_table]
#' [mech_stat_cyp_risk()]
#' [mech_stat_cyp_risk_table()]
#'
#' @details
#' CYP inhibition data can contain ki data for multiple compounds.
#' \preformatted{
#'       name   param      value    source
#'         M1    name         M1
#'         M1  CYP2C9        4.4 study 002
#' examplinib    name examplinib
#' examplinib  CYP1A2         NA
#' examplinib  CYP2B6         NA
#' examplinib  CYP2C8         11 study 001
#' examplinib  CYP2C9       13.5 study 001
#' examplinib CYP2C19         15 study 001
#' examplinib  CYP2D6         NA
#' examplinib  CYP3A4       12.5 study 001
#' }
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


#' Examplinib CYP induction data
#'
#' @details
#' \preformatted{
#'       name     cyp  emax ec50 maxc    source
#' examplinib  CYP1A2  1.00   NA    5 study 007
#' examplinib  CYP2B6  1.00   NA    5 study 007
#' examplinib  CYP2C8    NA   NA   NA
#' examplinib  CYP2C9    NA   NA   NA
#' examplinib CYP2C19    NA   NA   NA
#' examplinib  CYP2D6    NA   NA   NA
#' examplinib  CYP3A4  7.35 1.64    3 study 007
#'         M1  CYP1A2  1.00   NA    5 study 007
#'         M1  CYP2B6  6.98 1.86    5 study 007
#'         M1  CYP2C8    NA   NA   NA
#'         M1  CYP2C9    NA   NA   NA
#'         M1 CYP2C19    NA   NA   NA
#'         M1  CYP2D6    NA   NA   NA
#'         M1  CYP3A4 22.70 1.10    5 study 007
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_induction_data"


#' Examplinib transporter inhibition data
#'
#' @details
#' \preformatted{
#'       name   param      value    source
#' examplinib    name examplinib
#' examplinib     Pgp       0.41 study 005
#' examplinib    BCRP        1.9 study 005
#' examplinib    OCT1        2.3 study 006
#' examplinib OATP1B1        177 study 006
#' examplinib OATP1B3         35 study 006
#' examplinib    OAT1        271
#' examplinib    OAT3        300
#' examplinib    BSEP       12.8
#' examplinib    OCT2         67 study 006
#' examplinib   MATE1        3.6 study 006
#' examplinib  MATE2k        1.1 study 006
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_transporter_inhibition_data"


#' Examplinib UGT inhibition data
#'
#' @format
#' A data frame with the columns 'name', 'param', 'value' and 'source'
#' @details
#' \preformatted{
#'       name   param      value    source
#'         M1    name         M1
#'         M1  UGT1A1        1.1 study 009
#'         M1  UGT1A3        5.8 study 009
#'         M1  UGT1A4        6.2 study 009
#'         M1  UGT1A6         15 study 009
#'         M1  UGT1A9        3.6 study 009
#'         M1  UGT2B7         15 study 009
#'         M1 UGT2B15        9.6 study 009
#'         M1 UGT2B17        2.2 study 009
#' examplinib    name examplinib
#' examplinib  UGT1A1         15 study 009
#' examplinib  UGT1A3         15 study 009
#' examplinib  UGT1A4         15 study 009
#' examplinib  UGT1A6         15 study 009
#' examplinib  UGT1A9        3.8 study 009
#' examplinib  UGT2B7         15 study 009
#' examplinib UGT2B15         15 study 009
#' examplinib UGT2B17        6.1 study 009
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_ugt_inhibition_data"

