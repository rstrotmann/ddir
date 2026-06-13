# Clean, reproducible script for 3 literature-based DDI examples using ddir
# + compact comparison table: paper vs ddir
# -------------------------------------------------------------------------
# Sources:
# 1) Gomez-Mantilla et al. Clin Pharmacokinet 2023 (PMC10042977)
# 2) Bergagnini-Kolev et al. AAPS J 2023 (DOI: 10.1208/s12248-023-00828-z)
suppressPackageStartupMessages({
  library(pkgload)
  library(dplyr)
  library(tibble)
})
pkgload::load_all(".")
# helper for safe relative delta
rel_delta <- function(ddir_value, paper_value) {
  if (is.na(paper_value) || paper_value == 0) return(NA_real_)
  (ddir_value - paper_value) / paper_value
}
results <- list()

# -------------------------------------------------------------------------
# Example 1: Ivosidenib (CYP3A induction; no-gut scenario)
# Paper (MSM): AUCR ~ 0.46
# -------------------------------------------------------------------------
perp_ivo <- perpetrator(
  name       = "ivosidenib",
  oral       = FALSE,        # no-gut scenario
  mw         = 583,
  dose       = 500,
  imaxss     = 9.54 * 583,   # from Cavg,ss 9.54 uM
  fu         = 0.057,
  fumic      = 1,
  rb         = 0.56,
  fa         = 0.69,
  fg         = 1,
  ka         = 1.38 / 60,    # h^-1 to min^-1
  solubility = Inf
)

res_ivo <- mech_stat_cyp_risk(
  perp = perp_ivo,
  cyp_inh = inhibitor(data.frame(object = "CYP3A4", ki = NA_real_, source = "none")),
  cyp_ind = inducer(data.frame(object = "CYP3A4", emax = 21.25, ec50 = 8.6, max_c = 24.59, source = "paper")),
  cyp_tdi = inhibitor(NULL),
  include_induction = TRUE
)@table

ivo_ddir <- res_ivo$aucr[[1]]
ivo_paper <- 0.46
results[[length(results) + 1]] <- tibble(
  example = "Ivosidenib -> Midazolam (AUCR, no-gut)",
  metric = "AUCR",
  paper_value = ivo_paper,
  ddir_value = ivo_ddir,
  abs_delta = ddir_value - paper_value,
  rel_delta = rel_delta(ddir_value, paper_value)
)

# -------------------------------------------------------------------------
# Example 2: Voxelotor (CYP3A inhibition, scenario 2; no-gut)
# Paper (MSM): AUCR ~ 2.27
# -------------------------------------------------------------------------
perp_vox <- perpetrator(
  name       = "voxelotor",
  oral       = FALSE,         # no-gut scenario
  mw         = 337,
  dose       = 1500,
  imaxss     = 44.25 * 337,   # from Cavg,ss 44.25 uM
  fu         = 0.002,
  fumic      = 1,
  rb         = 1,
  fa         = 1,
  fg         = 1,
  ka         = 0.1,
  solubility = Inf
)

res_vox <- mech_stat_cyp_risk(
  perp = perp_vox,
  cyp_inh = inhibitor(data.frame(object = "CYP3A4", ki = 0.06, source = "scenario2")),
  cyp_ind = inducer(data.frame(object = "CYP3A4", emax = NA_real_, ec50 = NA_real_, max_c = NA_real_, source = "none")),
  cyp_tdi = inhibitor(NULL),
  include_induction = FALSE
)@table

vox_ddir <- res_vox$aucr[[1]]
vox_paper <- 2.27
results[[length(results) + 1]] <- tibble(
  example = "Voxelotor -> Midazolam (AUCR, scenario2 no-gut)",
  metric = "AUCR",
  paper_value = vox_paper,
  ddir_value = vox_ddir,
  abs_delta = ddir_value - paper_value,
  rel_delta = rel_delta(ddir_value, paper_value)
)

# -------------------------------------------------------------------------
# Example 3: PUR1900 inhaled itraconazole basic static
# Paper: combined R1 ~ 1.35 (ITZ + OH-ITZ)
# -------------------------------------------------------------------------
perp_pur <- perpetrator(
  name       = "PUR1900_itraconazole",
  oral       = TRUE,
  mw         = 705.6,
  dose       = 35,
  imaxss     = 0.0215 * 705.6,  # from Cmax ITZ 0.0215 uM
  fu         = 0.016,
  fumic      = 1,
  rb         = 0.58,
  fa         = 0.159,
  fg         = 0.58,
  ka         = 0.1,
  solubility = Inf
)
res_pur <- basic_cyp_inhibition_risk(
  perp_pur,
  inhibitor(data.frame(object = "CYP3A4", ki = 0.0013, source = "paper"))
)@table

# Parent-only R1 from ddir
R1_parent <- 1 + res_pur$r[[1]]
# Add OH-ITZ term manually according to paper method:
# R1 = 1 + (Imax,u/Ki,u)_ITZ + (Imax,u/Ki,u)_OH-ITZ
oh_term <- (0.0120 * 0.016) / 0.0023
R1_combined_ddir_style <- 1 + res_pur$r[[1]] + oh_term
R1_combined_paper <- 1.35
results[[length(results) + 1]] <- tibble(
  example = "PUR1900/ITZ -> Midazolam (combined R1)",
  metric = "R1",
  paper_value = R1_combined_paper,
  ddir_value = R1_combined_ddir_style,
  abs_delta = ddir_value - paper_value,
  rel_delta = rel_delta(ddir_value, paper_value)
)
# optional: include parent-only informational row
results[[length(results) + 1]] <- tibble(
  example = "PUR1900/ITZ -> Midazolam (parent-only R1)",
  metric = "R1",
  paper_value = NA_real_,
  ddir_value = R1_parent,
  abs_delta = NA_real_,
  rel_delta = NA_real_
)

# -------------------------------------------------------------------------
# Final compact summary
# -------------------------------------------------------------------------
summary_tbl <- bind_rows(results) |>
  mutate(
    paper_value = round(paper_value, 4),
    ddir_value  = round(ddir_value, 4),
    abs_delta   = round(abs_delta, 4),
    rel_delta_pct = round(100 * rel_delta, 2)
  ) |>
  select(example, metric, paper_value, ddir_value, abs_delta, rel_delta_pct)
print(summary_tbl, n = Inf)
