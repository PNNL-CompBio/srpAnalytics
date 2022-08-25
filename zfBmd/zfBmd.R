## Zebrafish Benchmark Database

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

### Confirm the counts are correct----------------------------------------------

to_pivot <- fread("~/Git_Repos/srpAnalytics/zfBmd/test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv")

to_pivot %>%
  pivot_wider(
    id_cols = c(chemical.id, bottle.id, conc, plate.id, well, date), 
    names_from = endpoint, 
    values_from = value
  ) %>%
  write.csv("~/Downloads/pivot_wider.csv", quote = F, row.names = F)

### Visualize LPR data 

toRM <- c('3756 0 20544 H02', '3756 0 20544 H08', '3756 0 20624 H05', 
          '3756 0 20624 H10', '3756 2.16 20544 G06', '3756 2.16 20625 G12', 
          '3756 6.8 20625 F04', '3756 14.7 20624 E07', '3756 14.7 20624 E08', 
          '3756 14.7 20625 E10', '3756 31.6 20544 D12', '3756 31.6 20624 D03', 
          '3756 31.6 20624 D11', '3756 56.2 20544 C02', '3756 56.2 20544 C04', 
          '3756 56.2 20544 C08', '3756 56.2 20625 C02', '3756 75 20544 B09', 
          '3756 75 20624 B02', '3756 75 20624 B03', '3756 75 20624 B04', 
          '3756 75 20625 B07', '3756 75 20625 B08', '3756 100 20544 A01', 
          '3756 100 20544 A03', '3756 100 20544 A07', '3756 100 20624 A09', 
          '3756 6.8 20544 F08', '3756 6.8 20625 F06', '3756 6.8 20625 F08', 
          '3756 6.8 20625 F10', '3756 14.7 20624 E03', '3756 14.7 20625 E04', 
          '3756 31.6 20544 D05', '3756 31.6 20624 D09', '3756 31.6 20625 D02', 
          '3756 31.6 20625 D03', '3756 31.6 20625 D08', '3756 56.2 20544 C09', 
          '3756 56.2 20624 C03', '3756 56.2 20624 C04', '3756 56.2 20624 C06', 
          '3756 56.2 20624 C10', '3756 56.2 20625 C10', '3756 75 20624 B01', 
          '3756 100 20624 A05', '3756 100 20625 A03')

LPR_timepoints <- fread("~/Git_Repos/srpAnalytics/zfBmd/test_files/7_PAH_zf_LPR_data_2021JAN11_3756.csv")

# Add some minimal colu
LPR_timepoints <- LPR_timepoints %>%
  mutate(
    time = (gsub("t", "", variable) %>% as.numeric()) / 10,
    TimeBin = floor(time),
    well.id = paste(chemical.id, conc, plate.id, well)
  ) %>%
  filter(
    well.id %in% toRM == FALSE
  ) %>%
  group_by(chemical.id, conc, plate.id, well, TimeBin) %>%
  summarise(
    Sum = sum(value, na.rm = T)
  ) %>%
  ungroup() 

ggplot(LPR_timepoints, aes(x = TimeBin, y = Sum)) + geom_smooth() + xlab("Minutes") + theme_bw() +
  geom_vline(xintercept = c(2, 3), color = "black") + facet_wrap(.~conc) 

# MOV1 (actually MOV2) is the movement at the start of light cycle 2
MOV_test <- LPR_timepoints %>%
  filter(TimeBin %in% c(2,3)) %>%
  group_by(chemical.id, conc, plate.id, well) %>%
  summarise(
    MOV = Sum[TimeBin == 3] - Sum[TimeBin == 2]
  ) 

test <- MOV_test %>%
  ungroup() %>%
  group_by(chemical.id, conc, plate.id) %>%
  mutate(
    Quart1 = quantile(MOV[which(conc == 0)], 0.25), 
    Quart3 = quantile(MOV[which(conc == 0)], 0.75),
    IQR_A = (Quart3 - Quart1) * 1.5
  ) %>%
  summarise(
    Hypo = sum(MOV <= 0, na.rm = T),
    Hyper = sum(MOV < (Quart1 - IQR_A) | MOV > (Quart3 + IQR_A), na.rm = T),
    Total = length(MOV)
  ) %>%
  ungroup() %>%
  group_by(chemical.id, conc) %>%
  summarise(
    Hypo = sum(Hypo, na.rm = T),
    Hyper = sum(Hyper, na.rm = T),
    Total = sum(Total),
    Response = (Hypo + Hyper) / Total
  )

