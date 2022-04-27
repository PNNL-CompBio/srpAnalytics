## Zebrafish Benchmark Database

library(data.table)
library(dplyr)
library(tidyr)

### Confirm the counts are correct----------------------------------------------

pivot_wider <- data.table::fread("~/Desktop/srpTargets/zfBmd/pivot_wider.csv")

pivot_wider2 <- data.table(pivot_wider[,1:4], do.call(cbind, lapply(5:22, function(col) {
  vals <- pivot_wider[,..col] %>% unlist()
  vals[vals > 1] <- 1
  return(vals)
})))
colnames(pivot_wider2) <- colnames(pivot_wider)

# Determine which data to filter based on NAs at 0 concentration - absolutely none should be filtered out
FilterLow <- pivot_wider2 %>%
  select(-well) %>%
  pivot_longer(colnames(pivot_wider)[5:22]) %>%
  group_by(chemical.id, conc, plate.id, name) %>%
  summarise(
    NACount = sum(is.na(value)),
    Total = length(value),
    Flag = NACount >= Total * 0.5
  )


Lies <- pivot_wider %>%
  select(-c(plate.id, well)) %>%
  pivot_longer(colnames(pivot_wider)[5:22]) %>%
  group_by(chemical.id, conc, name) %>%
  summarise(
    NumAffected = sum(value, na.rm = T),
    TotWells = length(value),
    NumEmbryos = length(value[!is.na(value)]),
    FractAffected = NumAffected / NumEmbryos
  ) %>% 
  arrange(name)

Truth <- pivot_wider2 %>%
  select(-c(plate.id, well)) %>%
  pivot_longer(colnames(pivot_wider2)[5:22]) %>%
  group_by(chemical.id, conc, name) %>%
  summarise(
    NumAffected = sum(value, na.rm = T),
    TotWells = length(value),
    NumEmbryos = length(value[!is.na(value)]),
    FractAffected = NumAffected / NumEmbryos
  ) %>% 
  arrange(name)

write.csv(RealNumAffected, "~/Desktop/srpTargets/zfBmd/RealNumAffected.csv")



