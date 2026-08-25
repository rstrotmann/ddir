#' Drug transporter reference data
#'
#' @format
#' * 'transporter' The name of the drug transporter protein.
#' * 'threshold' The regulatory threshold for clinically relevant interactions.
#' * 'i' The precipitant concentration metric applicable for the interaction.
#' @details
#' \preformatted{
#' transporter , threshold , i
#' Pgp_int     , 10        , igut
#' Pgp_sys     , 0.02      , imaxssu
#' BCRP_int    , 10        , igut
#' BCRP_sys    , 0.02      , imaxssu
#' OATP1B1     , 0.1       , imaxinletu
#' OATP1B3     , 0.1       , imaxinletu
#' OAT1        , 0.1       , imaxssu
#' OAT3        , 0.1       , imaxssu
#' BSEP        , 0.1       , imaxssu
#' OCT2        , 0.1       , imaxssu
#' MATE1       , 0.02      , imaxssu
#' MATE2k      , 0.02      , imaxssu
#' }
#' @source Section 2 of the [ICH M12 guideline](https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf).
#' @seealso [ddir::key_conc_table()]
#' @seealso [ddir::transporter_inh_risk()]
"transporter_reference_data"


#' CYP reference substrate data
#'
#' CYP reference substrates commonly used in the mechanistic static assessment
#' of the CYP DDI potential of drugs.
#'
#' @details
#' The CYP reference substrates currently implemented include:
#' \preformatted{
#' cyp     , substrate   , fgut , fm   , fmcyp
#' CYP1A2  , tizanidine  , 1    , 0.95 , 0.98
#' CYP2C8  , repaglinide , 1    , 1    , 0.61
#' CYP2C9  , S-warfarin  , 1    , 1    , 0.91
#' CYP2C19 , omeprazole  , 1    , 1    , 0.87
#' CYP3A4  , midazolam   , 0.57 , 0.96 , 1
#' }
#' @source FDA and EMA guidelines.
#' @seealso [ddir::mech_stat_cyp_risk()]
"cyp_reference_substrates"


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


#' Examplinib CYP inhibition data
#'
#' @details
#' CYP inhibition data for examplinib
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_inh_parent"


#' Examplinib metabolite CYP inhibition data
#'
#' @details
#' CYP inhibition data for examplinib metabolite
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_inh_metabolite"


#' Examplinib CYP TDI data
#'
#' @details
#' \preformatted{
#'       name    cyp   ki kinact    source
#' examplinib CYP3A4 30.7   0.04 study 001
#' }
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_tdi_parent"


#' Examplinib CYP induction data
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_ind_parent"


#' Examplinib metabolite CYP induction data
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_cyp_ind_metabolite"


#' Examplinib transporter inhibition data
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_transporter_inh_parent"


#' Examplinib UGT inhibition data
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_ugt_inh_parent"

#' Examplinib metabolite UGT inhibition data
#'
#' @source Fictional data, made up for demo purposes.
"examplinib_ugt_inh_metabolite"


#' Hepatic CYP turnover data based on various publications
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Please use [ddir::cyp_turnover()] instead.
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
#' @seealso [ddir::basic_cyp_tdi_risk()]
#' @seealso [ddir::mech_stat_cyp_risk()]
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


#' in vitro sample data for CYP3A4 TDI
#'
#' @description
#' CYP3A4 activity data at different pre-incubation timees and
#' different concentrations of a time-dependent inhibitor.
#'
#' @format
#' A data frame with 4 columns
#' * TIME Pre-incubation time in minutes.
#' * CONC Precipitant concentration in uM.
#' * ACT Enzyme activity in percent of initial.
#' * SOURCE Source study report
#'
#' @source Fully synthetic data set to demonstrate TDI-related calculations
"examplinib_in_vitro_tdi"


#' in vitro sample data for CYP3A4 mRNA induction
#'
#' @description
#' CYP3A4 mRNA fold-change data at different precipitant concentrations
#'
#' @format
#' A data frame with 7 columns
#' * DONOR Hepatocyte donor
#' * TYPE Test  or positive control
#' * CONC Precipitant concentration
#' * OBJECT DDI target
#' * FOLD mRNA fold change
#' * REL mRNA induction in percent relative to positive control
#' * SOURCE Source study report
#'
#' @source Fully synthetic data set to demonstrate CYP mRNA induction data
#' evaluation.
"examplinib_in_vitro_ind"


#' in vitro sample data for CYP3 mRNA induction
#'
#' @description
#' CYP3A4 mRNA fold-change data at different precipitant concentrations
#'
#' @format
#' A data frame with 7 columns
#' * DONOR Hepatocyte donor
#' * TYPE Test  or positive control
#' * CONC Precipitant concentration
#' * OBJECT DDI target
#' * FOLD mRNA fold change
#' * REL mRNA induction in percent relative to positive control
#' * SOURCE Source study report
#'
#' @source Fully synthetic data set to demonstrate CYP mRNA induction data
#' evaluation.
"examplinib_in_vitro_ind1"
