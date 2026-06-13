## in vitro CYP3A4 activity after varying pre-incubation times
## fully synthetic data

examplinib_in_vitro_tdi <- tibble::tribble(
  ~TIME, ~CONC_0.2, ~CONC_0.66, ~CONC_2, ~CONC_6.66, ~CONC_20, ~CONC_50, ~CONC_100,
      0,      89.5,       93.5,   112.7,       98.5,     87.4,     91.1,      68.2,
    2.5,      87.7,       97.6,   112.1,       98.8,     82.6,     85.1,      56.1,
      5,      87.3,       95.7,   114.7,      106.4,     87.7,     82.3,      58.3,
     10,      88.6,      101.4,   120.4,      106.4,     85.4,     70.4,      48.8,
     20,      97.4,      100.4,   113.6,       95.7,     70.8,     53.5,      34.2,
     30,      88.6,       93.8,   115.9,         97,     50.3,     52.4,      28.1
  ) |>
  pivot_longer(cols = 2:8, names_to = "GROUP", values_to = "ACT") |>
  separate(col = GROUP, sep = "_", into = c("PAR", "CONC"), remove = FALSE) |>
  mutate(CONC = as.numeric(CONC)) |>
  select(-c("GROUP", "PAR")) |>
  mutate(SOURCE = "Study examplinib_tdi")

usethis::use_data(examplinib_in_vitro_tdi, overwrite = TRUE)


## in vitro CYP3A4 mRNA induction data, fully synthetic

examplinib_in_vitro_ind <- tibble::tribble(
  ~DONOR,       ~PRECIPITANT, ~CONC,  ~OBJECT,  ~FOLD,  ~REL,                ~SOURCE,
     "A",             "test",  0.01, "CYP3A4",  0.904, -0.16, "Study examplinib_ind",
     "A",             "test",  0.03, "CYP3A4",   0.94,  -0.1, "Study examplinib_ind",
     "A",             "test",   0.1, "CYP3A4",  1.066,  0.11, "Study examplinib_ind",
     "A",             "test",   0.3, "CYP3A4",  1.182,   0.3, "Study examplinib_ind",
     "A",             "test",   0.5, "CYP3A4",  1.451,  0.74, "Study examplinib_ind",
     "A",             "test",     1, "CYP3A4",  1.943,  1.55, "Study examplinib_ind",
     "A",             "test",     3, "CYP3A4",  6.331,  8.78, "Study examplinib_ind",
     "A",             "test",     5, "CYP3A4",     NA,    NA, "Study examplinib_ind",
     "A", "positive_control",    NA, "CYP3A4", 61.699,   100, "Study examplinib_ind",
     "B",             "test",  0.01, "CYP3A4",  0.728, -2.58, "Study examplinib_ind",
     "B",             "test",  0.03, "CYP3A4",  0.881, -1.13, "Study examplinib_ind",
     "B",             "test",   0.1, "CYP3A4",  0.796, -1.93, "Study examplinib_ind",
     "B",             "test",   0.3, "CYP3A4",  0.655, -3.27, "Study examplinib_ind",
     "B",             "test",   0.5, "CYP3A4",  0.751, -2.36, "Study examplinib_ind",
     "B",             "test",     1, "CYP3A4",  0.844, -1.48, "Study examplinib_ind",
     "B",             "test",     3, "CYP3A4",  1.325,  3.08, "Study examplinib_ind",
     "B",             "test",     5, "CYP3A4",   2.31, 12.41, "Study examplinib_ind",
     "B", "positive_control",    NA, "CYP3A4", 11.552,   100, "Study examplinib_ind",
     "C",             "test",  0.01, "CYP3A4",  0.892, -1.31, "Study examplinib_ind",
     "C",             "test",  0.03, "CYP3A4",  0.796, -2.48, "Study examplinib_ind",
     "C",             "test",   0.1, "CYP3A4",  0.754, -2.99, "Study examplinib_ind",
     "C",             "test",   0.3, "CYP3A4",  0.711, -3.51, "Study examplinib_ind",
     "C",             "test",   0.5, "CYP3A4",  0.648, -4.28, "Study examplinib_ind",
     "C",             "test",     1, "CYP3A4",  0.728, -3.31, "Study examplinib_ind",
     "C",             "test",     3, "CYP3A4",    1.2,  2.43, "Study examplinib_ind",
     "C",             "test",     5, "CYP3A4",  1.585,  7.11, "Study examplinib_ind",
     "C", "positive_control",    NA, "CYP3A4",  9.223,   100, "Study examplinib_ind",
     "A",             "test",   0.3, "CYP1A2",   1.38,   5.8, "Study examplinib_ind",
     "A",             "test",     1, "CYP1A2",   2.07,   8.7, "Study examplinib_ind",
     "A",             "test",     3, "CYP1A2",   3.08, 12.94, "Study examplinib_ind",
     "A",             "test",    10, "CYP1A2",   7.37, 30.97, "Study examplinib_ind",
     "A", "positive_control",    NA, "CYP1A2",   23.8,   100, "Study examplinib_ind",
     "B",             "test",   0.3, "CYP1A2",   1.24,  4.28, "Study examplinib_ind",
     "B",             "test",     1, "CYP1A2",   1.94,  6.69, "Study examplinib_ind",
     "B",             "test",     3, "CYP1A2",   2.86,  9.86, "Study examplinib_ind",
     "B",             "test",    10, "CYP1A2",   5.09, 17.55, "Study examplinib_ind",
     "B", "positive_control",    NA, "CYP1A2",     29,   100, "Study examplinib_ind",
     "C",             "test",   0.3, "CYP1A2",   1.56,  4.47, "Study examplinib_ind",
     "C",             "test",     1, "CYP1A2",   2.25,  6.45, "Study examplinib_ind",
     "C",             "test",     3, "CYP1A2",   3.21,   9.2, "Study examplinib_ind",
     "C",             "test",    10, "CYP1A2",   4.15, 11.89, "Study examplinib_ind",
     "C", "positive_control",    NA, "CYP1A2",   34.9,   100, "Study examplinib_ind"
  )

usethis::use_data(examplinib_in_vitro_ind, overwrite = TRUE)
