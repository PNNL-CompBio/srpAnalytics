## Zebrafish Benchmark Database

library(data.table)
library(dplyr)
library(tidyr)
library(xlsx)

### Confirm the counts are correct----------------------------------------------

to_pivot <- fread("~/Git_Repos/srpAnalytics/zfBmd/test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv")

to_pivot %>%
  pivot_wider(
    id_cols = c(chemical.id, bottle.id, conc, plate.id, well, date), 
    names_from = endpoint, 
    values_from = value
  ) %>%
  write.csv("~/Downloads/pivot_wider.csv", quote = F, row.names = F)




