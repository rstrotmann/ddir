## code to prepare `hepatic_cyp_turnover` dataset goes here
library(tidyverse)
library(usethis)

hepatic_cyp_turnover <- tribble(
  ~cyp, ~method, ~mean_hl, ~in_vivo,~reference,
  "CYP1A2", "In vitro Method 1", "51", "FALSE", "PMID: 2136526",
  "CYP1A2", "In vitro Method 2", "43", "FALSE", "DOI: 10.1007/3-540-29804-5_25",
  "CYP1A2", "In vitro Method 2", "36", "FALSE", "PMID: 10997941",
  "CYP1A2", "In vivo Method 1", "39", "TRUE", "DOI: 10.1016/j.clpt.2004.04.003",
  "CYP1A2", "In vivo Method 3", "105", "TRUE", "DOI: 10.1038/sj.clpt.6100431",
  "CYP2A6", "In vitro Method 2", "226", "FALSE", "PMID: 10997941",
  "СУР2В6", "In vitro Method 2", "32", "FALSE", "PMID: 10997941",
  "CYP2C8", "In vitro Method 2", "23", "FALSE", "PMID: 10997941",
  "СУР2С9", "In vitro Method 2", "104", "FALSE", "PMID: 10997941",
  "CYP2C19", "In vitro Method 2", "26", "FALSE", "PMID: 10997941",
  "CYP2D6", "In vitro Method 2", "70", "FALSE", "PMID: 10997941",
  "CYP2D6", "In vivo Method 2", "51", "TRUE", "DOI: 10.1097/00004714-200204000-00010",
  "CYP2E1", "In vitro Method 2", "27", "FALSE", "PMID: 10997941",
  "CYP2E1", "In vivo Method 1", "60", "TRUE", "DOI: 10.1111/j.1530-0277.1995.tb01516.x",
  "CYP2E1", "In vivo Method 2", "50", "TRUE", "PMID: 10490907.",
  "CYP3A4", "In vitro Method 1", "44", "FALSE", "PMID: 1614409",
  "CYP3A4", "In vitro Method 2", "26", "FALSE", "DOI: 10.1007/3-540-29804-5_25",
  "CYP3A4", "In vitro Method 2", "79", "FALSE", "PMID: 10997941",
  "CYP3A4", "In vivo Method 1", "72", "TRUE", "PMID: 10234596",
  "CYP3A4", "In vivo Method 3", "96", "TRUE", "",
  "CYP3A4", "In vivo Method 1", "72", "TRUE", "DOI: 10.1002/cpt1978243316",
  "CYP3A4", "In vivo Method 3", "10", "TRUE", "DOI: 10.1007/BF00685732",
  "CYP3A4", "In vivo Method 3", "94", "TRUE", "DOI: 10.1046/j.1365-2125.1999.00974.x",
  "CYP3A4", "In vivo Method 3", "70", "TRUE", "DOI: 10.1038/sj.clpt.6100431",
  "CYP3A4", "In vivo Method 3", "85", "TRUE", "DOI: 10.1128/AAC.41.5.898",
  "CYP3A4", "In vivo Method 3", "140", "TRUE", "DOI: 10.1016/S0009-9236(98)90018-2",
  "CYP35", "In vitro Method 2", "36", "FALSE", "PMID: 10997941"
) %>%
  mutate_all(str_trim) %>%
  as.data.frame() %>%
  mutate(in_vivo = str_detect(method, "In vivo")) %>%
  mutate(mean_hl=as.numeric(mean_hl)) %>%
  mutate(kdeg=round(log(2)/mean_hl, 4))

usethis::use_data(hepatic_cyp_turnover, overwrite = TRUE)
