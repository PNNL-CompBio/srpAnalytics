## Zebrafish Benchmark Database

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

##########################
## PRE-PROCESS & FILTER ##
##########################

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


#############
## QC PLOT ##
#############

Summed <- LPR_timepoints %>%
  mutate(
    time = (gsub("t", "", variable) %>% as.numeric()) / 10,
    TimeBin = floor(time),
    well.id = paste(chemical.id, conc, plate.id, well),
    plate.id = as.factor(plate.id)
  ) %>%
  filter(
    well.id %in% toRM == FALSE
  ) %>%
  group_by(chemical.id, conc, plate.id, well, TimeBin) %>%
  summarise(
    Sum = sum(value, na.rm = T)
  ) %>%
  ungroup() 

ggplot(Summed, aes(x = TimeBin, y = Sum, color = plate.id)) + ylab("Movement") +
  xlab("Minutes") + theme_bw() + ggtitle("Chemical ID = 3756. Wells are summed together.") +
  annotate(geom = "rect", xmin = c(0, 6, 12, 18), xmax = c(2, 8, 14, 20), 
           ymin = rep(-5, 4), ymax = rep(380, 4),
           fill = "yellow", color = NA, alpha = 0.75) +
  annotate(geom = "rect", xmin = c(3, 9, 15, 21), xmax = c(5, 11, 17, 23), 
           ymin = rep(-5, 4), ymax = rep(380, 4),
           fill = "gray50", color = NA, alpha = 0.5) +
  geom_line() + facet_wrap(~conc)


## Read LPR done 
LPR <- read.csv("~/Downloads/LPR_done.csv") %>%
  mutate(Response = num.affected / num.nonna)

ggplot(LPR, aes(x = conc, y = Response)) + geom_line() + xlab("Concentration") + facet_wrap(.~endpoint)
