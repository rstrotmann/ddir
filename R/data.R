#' Drug transporter reference data
#'
#' @format
#' * 'transporter' The name of the drug transporter protein.
#' * 'fda_thld' The threshold for relevant interactions as defined by the
#' [FDA guidance](https://www.fda.gov/media/134582/download).
#' * 'ema_thld' the threshold for relevant interactions as defined in the
#' [EMA guidance](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-drug-interactions-revision-1_en.pdf).
#' * 'i' The perpetrator concentration applicable for the interaction.
#' @details
#' \preformatted{
#' transporter   fda_thld   ema_thld    i
#' "Pgp_int",    10.0,      10.00,     "igut",
#' "Pgp_sys",    0.1,       0.02,      "imaxssu",
#' "BCRP_int",   10.0,      10.00,     "igut",
#' "BCRP_sys",   0.1,       0.02,      "imaxssu",
#' "OCT1",       NA,        0.04,      "imaxinletu",
#' "OATP1B1",    0.1,       0.04,      "imaxinletu",
#' "OATP1B3",    0.1,       0.04,      "imaxinletu",
#' "OAT1",       0.1,       0.04,      "imaxssu",
#' "OAT3",       0.1,       0.04,      "imaxssu",
#' "BSEP",       0.1,       0.02,      "imaxssu",
#' "OCT2",       0.1,       0.02,      "imaxssu",
#' "MATE1",      0.1,       0.02,      "imaxssu",
#' "MATE2k",     0.1,       0.02,      "imaxssu"
#' }
#' @source [FDA guideline: In Vitro Drug Interaction Studies — Cytochrome P450
#' Enzyme- and Transporter-Mediated Drug Interactions, 2020](https://www.fda.gov/media/134582/download).
#' @source [EMA Guideline on the investigation of drug interactions, 2012](https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-drug-interactions-revision-1_en.pdf)
#' @seealso [key_concentrations()]
#' @seealso [transporter_inhibition_risk()]
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
#' @noRd
#' @keywords internal
"examplinib_compounds_string"


#' Examplinib CYP inhibitor data as string
#'
#' A character string containing CYP inhibition data for examplinib and M1.
#' @source Fictional data for demo purposes.
#' @seealso [read_inhibitor_data()]
#' @noRd
"examplinib_cyp_inhibition_string"


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


#' Examplinib CYP TDI data as string
#'
#' A character string containing CYP TDI data for examplinib.
#' @source Fictional data for demo purposes.
#' @seealso [read_tdi_data()]
#' @noRd
"examplinib_cyp_tdi_string"


#' Examplinib CYP TDI data
#'
#' @details
#' \preformatted{
#'       name    cyp   ki kinact    source
#' examplinib CYP3A4 30.7   0.04 study 001
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_tdi_data"


#' Examplinib CYP induction data as string
#'
#' A character string containing CYP induction data for examplinib and M1.
#' @source Fictional data for demo purposes.
#' @seealso [read_inducer_data()]
#' @noRd
"examplinib_cyp_induction_string"


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


#' Examplinib transporter inhibition data as string
#'
#' A character string containing transporter inhibition data for examplinib and M1.
#' @source Fictional data for demo purposes.
#' @seealso [read_transporter_inhibitor_data()]
#' @noRd
"examplinib_transporter_inhibition_string"


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


#' Examplinib UGT inhibition data as string
#'
#' A character string containing UGT inhibition data for examplinib and M1.
#' @source Fictional data for demo purposes.
#' @seealso [read_ugt_inhibitor_data()]
#' @noRd
"examplinib_ugt_inhibition_string"


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


#' Hepatic CYP turnover data based on various publications
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Please use [`cyp_turnover`] instead.
#'
#' @format
#' A data frame with 6 columns:
#' * cyp: The CYP enzymne
#' * method: The experimental method used (below)
#' * mean_hl: The mean CYP degradation half-life measured,
#' * in_vivo: Study was conducted in vivo
#' * reference: The source publication (PMID or DOI)
#' * kdeg: The CYP degradation constant, i.e., \eqn{log(2)/mean half-life}.
#'
#' @details
#' The following experimental methods were used in the original publication:
#' * in vitro method 1: Radio-labeling of enzyme (‘pulse-chase’ method)
#' * in vitro method 2: Degradation of enzyme in cultured hepatocytes or liver
#' slices
#' * in vivo method 1: Recovery of enzyme activity after enzyme induction
#' * in vivo method 2: Recovery of enzyme activity after mechanism-based
#' inhibition (MBI)
#' * in vivo method 3: Pharmacokinetic modeling of auto-induction
#' @details
#' These are the first few lines of the data frame:
#' \preformatted{
#'    cyp            method mean_hl in_vivo                       reference   kdeg
#' CYP1A2 In vitro Method 1      51   FALSE                   PMID: 2136526 0.0136
#' CYP1A2 In vitro Method 2      43   FALSE   DOI: 10.1007/3-540-29804-5_25 0.0161
#' CYP1A2 In vitro Method 2      36   FALSE                  PMID: 10997941 0.0193
#' CYP1A2  In vivo Method 1      39    TRUE DOI: 10.1016/j.clpt.2004.04.003 0.0178
#' CYP1A2  In vivo Method 3     105    TRUE    DOI: 10.1038/sj.clpt.6100431 0.0066
#' CYP2A6 In vitro Method 2     226   FALSE                  PMID: 10997941 0.0031
#' }
#' @source This data set is taken from:
#'
#' Yang J, Liao M, Shou M, Jamei M, Yeo KR, Tucker GT, Rostami-Hodjegan A.
#' Cytochrome p450 turnover: regulation of synthesis and degradation, methods
#' for determining rates, and implications for the prediction of drug
#' interactions.
#' Curr Drug Metab. 2008 Jun;9(5):384-94.
#' doi: 10.2174/138920008784746382. PMID: 18537575.
"hepatic_cyp_turnover"


#' CYP turnover rate constants
#'
#' @description
#' Mean degradation constants in 1/h for hepatic and intestinal CYP enzymes.
#' @seealso [basic_cyp_tdi_risk()]
#' @seealso [mech_stat_cyp_risk()]
#'
#' @format
#' A data frame with 3 columns
#'
#' @details
#' \preformatted{
#'     cyp kdeg_hepatic kdeg_intestinal
#'  CYP1A1       0.0183              NA
#'  CYP1A2       0.0183              NA
#'  CYP2A6       0.0267              NA
#'  CYP2B6       0.0217              NA
#'  CYP2C8       0.0301              NA
#'  CYP2C9       0.0067            0.03
#' CYP2C18       0.0267              NA
#' CYP2C19       0.0267            0.03
#'  CYP2D6       0.0099            0.03
#'  CYP2E1       0.0176              NA
#'  CYP2J2       0.0194            0.03
#'  CYP3A4       0.0193            0.03
#'  CYP3A5       0.0193            0.03
#'  CYP3A7       0.0193              NA
#' }
#' @source This data set is taken from various literature sources
"cyp_turnover"


