#' CYP induction by calcitriol
#'
#' Source:
#' Chae YJ, Kim MS, Chung SJ, Lee MK, Lee KR, Maeng HJ. Pharmacokinetic
#' Estimation Models-based Approach to Predict Clinical Implications for CYP
#' Induction by Calcitriol in Human Cryopreserved Hepatocytes and HepaRG Cells.
#' Pharmaceutics. 2021 Jan 29;13(2):181. doi: 10.3390/pharmaceutics13020181.
#' PMID: 33572963; PMCID: PMC7911399.
#'

cyp_ind_calcitriol <- tibble::tribble(
  ~OBJECT,  ~DONOR,          ~CONC, ~FOLD,  ~SD,
  # (a) CYP2B6  — HepaRG from text; PHH digitized
  "CYP2B6", "Hepatocyte #1",    1,  1.24,  0.17,
  "CYP2B6", "Hepatocyte #1",    5,  0.97,  0.57,
  "CYP2B6", "Hepatocyte #1",   10,  1.24,  0.24,
  "CYP2B6", "Hepatocyte #1",  100,  1.64,  0.18,
  "CYP2B6", "Hepatocyte #2",    1,  0.98,    NA,
  "CYP2B6", "Hepatocyte #2",    5,  1.05,    NA,
  "CYP2B6", "Hepatocyte #2",   10,  1.49,  0.02,  # text
  "CYP2B6", "Hepatocyte #2",  100,  1.21,    NA,
  "CYP2B6", "Hepatocyte #3",    1,  1.13,    NA,
  "CYP2B6", "Hepatocyte #3",    5,  1.40,    NA,
  "CYP2B6", "Hepatocyte #3",   10,  1.52,    NA,
  "CYP2B6", "Hepatocyte #3",  100,  1.46,    NA,
  "CYP2B6", "HepaRG",           1,  1.70,  0.15,  # text
  "CYP2B6", "HepaRG",           5,  2.76,  0.56,
  "CYP2B6", "HepaRG",          10,  2.96,  0.39,
  "CYP2B6", "HepaRG",         100,  6.51,  1.14,
  # (b) CYP3A4  — digitized from Fig 3b
  "CYP3A4", "Hepatocyte #1",    1,  5.56,    NA,
  "CYP3A4", "Hepatocyte #1",    5,  6.33,    NA,
  "CYP3A4", "Hepatocyte #1",   10,  9.31,    NA,
  "CYP3A4", "Hepatocyte #1",  100, 33.9,   2.1,
  "CYP3A4", "Hepatocyte #2",    1,  1.8,     NA,
  "CYP3A4", "Hepatocyte #2",    5,  1.7,     NA,
  "CYP3A4", "Hepatocyte #2",   10,  1.9,     NA,
  "CYP3A4", "Hepatocyte #2",  100,  3.1,     NA,
  "CYP3A4", "Hepatocyte #3",    1,  1.8,     NA,
  "CYP3A4", "Hepatocyte #3",    5,  3.7,     NA,
  "CYP3A4", "Hepatocyte #3",   10,  4.0,     NA,
  "CYP3A4", "Hepatocyte #3",  100,  7.9,     NA,
  "CYP3A4", "HepaRG",           1,  6.7,     NA,
  "CYP3A4", "HepaRG",           5, 11.7,    4.8,
  "CYP3A4", "HepaRG",          10, 13.4,    2.9,
  "CYP3A4", "HepaRG",         100, 27.6,    1.5,
  # (c) CYP2C8  — digitized; 100 nM HepaRG = 4.31 in text
  "CYP2C8", "Hepatocyte #1",    1,  1.4,     NA,
  "CYP2C8", "Hepatocyte #1",    5,  2.31,   0.6,
  "CYP2C8", "Hepatocyte #1",   10,  2.3,     NA,
  "CYP2C8", "Hepatocyte #1",  100,  2.65,   0.2,
  "CYP2C8", "Hepatocyte #2",    1,  1.6,     NA,
  "CYP2C8", "Hepatocyte #2",    5,  1.80,    NA,
  "CYP2C8", "Hepatocyte #2",   10,  2.02,    NA,
  "CYP2C8", "Hepatocyte #2",  100,  2.10,    NA,
  "CYP2C8", "Hepatocyte #3",    1,  1.1,     NA,
  "CYP2C8", "Hepatocyte #3",    5,  1.46,    NA,
  "CYP2C8", "Hepatocyte #3",   10,  1.57,    NA,
  "CYP2C8", "Hepatocyte #3",  100,  1.35,    NA,
  "CYP2C8", "HepaRG",           1,  1.4,     NA,
  "CYP2C8", "HepaRG",           5,  1.8,     NA,
  "CYP2C8", "HepaRG",          10,  2.29,    NA,
  "CYP2C8", "HepaRG",         100,  4.31,   0.8,  # text 4.31
) |>
  mutate(SAMPLE = "test", SOURCE = "Chae et al., Pharmaceutics 2021") |>
  mutate(REL = NA)


cyp_ind_calcitriol |>
  ggplot(aes(x = CONC, y = FOLD, color = DONOR)) +
  geom_point() +
  geom_errorbar(aes(ymin = FOLD - SD, ymax = FOLD + SD), width = 0.1) +
  geom_line() +
  scale_x_log10(limits = c(0.1, 100)) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~OBJECT, scales = "free") +
  theme_bw()


x <- induction_experiment(cyp_ind_calcitriol, "calcitriol")

data <- mutate(x@data, ID = paste0(.data$OBJECT, "_", .data$DONOR))

sigm <- function(c, emax, ec50) {
  1 + (emax / (1 + exp(log(ec50)-log(c))))
}

out <- list()

temp <- data |>
  filter(.data$SAMPLE == "test") |>
  nest_by(DONOR, OBJECT, ID) #|>
  # mutate(emax_obs = max(data$FOLD, na.rm = TRUE))

out$data <- temp |>
  mutate(mod = list(
    nlsLM(
      FOLD ~ sigm(CONC, emax, ec50),
      data = data,
      start = list(ec50 = .1, emax = 1),
      lower = c(ec50 = 0, emax = 0),
      upper = c(ec50 = 100, emax = 100),
      control = nls.lm.control(maxiter = 1000)
    )
  )) |>
  mutate(modpar = list(broom::tidy(mod)))

pred <- data.frame(
  CONC = 10^seq(
    log10(min(data$CONC, na.rm = TRUE)),
    log10(max(data$CONC, na.rm = TRUE)),
    length.out = 100)
)

temp <- lapply(out$data$mod, function(x) predict(x, newdata = pred))
names(temp) <- out$data$ID
pred <- bind_cols(pred, temp) |>
  pivot_longer(cols = -1, names_to = "ID", values_to = "FOLD") |>
  separate(ID, c("OBJECT", "DONOR"), "_", remove = FALSE)

out$fold_plot <- ggplot(data = NULL, aes(x = CONC, y = FOLD, color = OBJECT)) +
  geom_line(data = pred) +
  geom_point(data = filter(data, SAMPLE == "test", !is.na(FOLD)), size = 2) +
  facet_wrap(~ID, scales = "free") +
  scale_x_log10() +
  expand_limits(y = 0) +
  labs(title = paste("In vitro CYP induction by", x@precipitant)) +
  theme_bw() +
  theme(legend.position = "none")

out$ind_param <- out$data |>
  unnest(modpar) |>
  ungroup() |>
  select(-c("ID", "data", "mod"))

# make inducer object from the donor that has the respective highers emax
max_c <- data |>
  reframe(
    max_c = max(CONC, na.rm = TRUE), .by = c("ID", "SOURCE"))

temp <- out$data |>
  unnest(modpar) |>
  select(-c("data", "mod")) |>
  group_by(OBJECT) |>
  # filter(emax_obs == max(emax_obs, na.rm = TRUE)) |>
  filter(term == "ec50") |>
  left_join(max_c, by = "ID") |>
  select(object = OBJECT, emax = emax, ec50 = estimate, max_c, source = SOURCE) |>
  mutate(ec50 = round(ec50, 2))

out$inducer <- inducer(temp, precipitant = x@precipitant)

out
