## in vitro CYP3A4 activity after varying pre-incubation times
## fully synthetic data

in_vitro_сyp3a4_act_sample_data <- tibble::tribble(
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
  select(-c("GROUP", "PAR"))

usethis::use_data(in_vitro_сyp3a4_act_sample_data, overwrite = TRUE)
