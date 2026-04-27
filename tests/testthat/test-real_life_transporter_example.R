# Transporter inhibition analogue script (3 drugs)
# ------------------------------------------------
# Case study: Baricitinib (OAT3 victim) with:
# - Probenecid 1000 mg BID
# - Ibuprofen 400 mg QD
# - Diclofenac 100 mg BID
#
# Literature basis:
# - Gomez-Mantilla et al., Clin Pharmacokinet 2023 (PMC10042977), Table 17
# - Posada et al., Clin Transl Sci 2017 (baricitinib transporter DDI context)
#
# IMPORTANT:
# This script uses a consistency-check workflow:
# 1) Take reported MSM AUCR (ft = 0.44) from the paper
# 2) Back-calculate Iu/IC50
# 3) Feed those into ddir::transporter_inh_risk() as perpetrator profiles
# So this verifies internal consistency with ddir's transporter ratio calculations.
suppressPackageStartupMessages({
  library(pkgload)
  library(dplyr)
  library(tibble)
})
pkgload::load_all(".")
# ---- helpers ----
infer_i_over_ic50 <- function(aucr, ft = 0.44) {
  # AUCR = 1 / ((1-ft) + ft/(1 + I/Ki))
  d <- 1 / aucr
  ft / (d - (1 - ft)) - 1
}
make_perp_from_Iu <- function(name, mw, Iu, fu = 1) {
  # Set imaxss so that imaxssu(perp, molar=TRUE) == Iu
  # imaxssu(molar=TRUE) = (imaxss * fu) / mw
  imaxss_required <- Iu * mw / fu
  perpetrator(
    name       = name,
    oral       = FALSE,       # focus on systemic OAT3
    mw         = mw,
    dose       = 1,
    imaxss     = imaxss_required,
    fu         = fu,
    fumic      = 1,
    rb         = 1,
    fa         = 1,
    fg         = 1,
    ka         = 0.1,
    solubility = Inf
  )
}
# ---- literature inputs ----
# AUCR predicted (MSM, ft=0.44) from Gomez-Mantilla 2023 Table 17
# IC50 values from transporter DDI literature context used in that section
drug_inputs <- tribble(
  ~drug,         ~mw,     ~ic50_oat3_uM, ~paper_aucr_msm_ft044,
  "probenecid",  285.36,  4.4,           1.67,
  "ibuprofen",   206.28,  4.4,           1.11,
  "diclofenac",  296.15,  3.8,           1.00
)
# ---- run consistency checks with ddir ----
results <- lapply(seq_len(nrow(drug_inputs)), function(i) {
  row <- drug_inputs[i, ]
  i_over_ic50 <- infer_i_over_ic50(row$paper_aucr_msm_ft044, ft = 0.44)
  Iu <- i_over_ic50 * row$ic50_oat3_uM
  perp <- make_perp_from_Iu(
    name = row$drug,
    mw   = row$mw,
    Iu   = Iu,
    fu   = 1
  )
  trans_inh <- inhibitor(
    data.frame(
      target = "OAT3",
      ic50   = row$ic50_oat3_uM,
      source = "literature"
    )
  )
  risk_tbl <- transporter_inh_risk(perp, trans_inh)@table |>
    filter(target == "OAT3")
  tibble(
    drug = row$drug,
    paper_aucr_msm_ft044 = row$paper_aucr_msm_ft044,
    inferred_Iu_uM = Iu,
    inferred_I_over_IC50 = i_over_ic50,
    ddir_r = risk_tbl$r[[1]],                # should match inferred_I_over_IC50
    ddir_threshold = risk_tbl$threshold[[1]],
    ddir_risk = risk_tbl$risk[[1]]
  )
})
summary_tbl <- bind_rows(results) |>
  mutate(
    across(where(is.numeric), ~ round(.x, 4)),
    delta_r = round(ddir_r - inferred_I_over_IC50, 6)
  ) |>
  select(
    drug,
    paper_aucr_msm_ft044,
    inferred_Iu_uM,
    inferred_I_over_IC50,
    ddir_r,
    delta_r,
    ddir_threshold,
    ddir_risk
  )
print(summary_tbl, n = Inf)
