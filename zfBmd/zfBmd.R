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


###################
## CALCULATE MOV ##
###################

# MOV is the difference between the start of the dark cycle and the end of the light cycle
# Here, the light cycles end at 2, 8, 14, and 20; and the dark cycles end at 5, 11, 17, 23. 

# MOV1 (actually MOV2) is the movement at the start of light cycle 2
MOV <- Summed %>%
  filter(TimeBin %in% c(2,3)) %>%
  group_by(chemical.id, conc, plate.id, well) %>%
  summarise(
    Response = Sum[TimeBin == 3] - Sum[TimeBin == 2]
  ) 

ggplot(MOV, aes(x = as.factor(conc), y = Response)) + geom_boxplot()

################################
## CONVERT MOV TO DICHOTOMOUS ##
################################

MOV_Final <- MOV %>%
  ungroup() %>%
  group_by(chemical.id, conc, plate.id) %>%
  mutate(
    Quart1 = quantile(Response[Response > 0], 0.25), 
    Quart3 = quantile(Response[Response > 0], 0.75),
    IQR_A = (Quart3 - Quart1) * 1.5, 
    Hypo = Response <= 0,
    Hyper = Response > 0 & (Response < (Quart1 - IQR_A) | Response > (Quart3 + IQR_A))
  ) %>%
  group_by(chemical.id, conc) %>%
  summarise(
    Hypo = sum(Hypo, na.rm = T),
    Hyper = sum(Hyper, na.rm = T),
    Total = length(Response),
    Response =  Hypo / Total
  )

# Add that MOV 1 is missing
ggplot(MOV_Final, aes(x = conc, y = Response)) + geom_line() + theme_bw() + 
  ggtitle("MOV")

###################
## CALCULATE AUC ##
###################


AUC <- do.call(rbind, lapply(1:Num_Cycles, function(num) {
  
  TimeLapse <- Assumed_Cycle_Time * (num - 1)

  Light <- LPR_timepoints %>% 
    group_by(chemical.id, plate.id, conc, well) %>%
    summarise(
      Value = sum(value[time %in% (Assumed_Light_Begin + TimeLapse):(Assumed_Light_End + TimeLapse)], na.rm = T)
    )
  
  Dark <-  LPR_timepoints %>% 
    group_by(chemical.id, plate.id, conc, well) %>%
    summarise(
      Value = sum(value[time %in% (Assumed_Dark_Begin + TimeLapse):(Assumed_Dark_End + TimeLapse)], na.rm = T)
    )

  Res <- Light[,c("chemical.id", "conc", "plate.id", "well")] %>%
    ungroup() %>%
    mutate(
      endpoint = paste0("AUC", num),
      Value = (Dark$Value - Light$Value) %>% unlist()
    )
  
  return(Res)
  
}))

ggplot(AUC, aes(x = as.factor(conc), y = Value, fill = endpoint)) + geom_boxplot()

################################
## CONVERT AUC TO DICHOTOMOUS ##
################################

# Negative control is the concentration at 0. If more than 50% of the abnormal (hypo) remain,
# then remove the entire plate 
AUC_0 <- AUC %>%
  filter(conc == 0) %>%
  group_by(chemical.id, endpoint) %>%
  mutate(
    Quart1 = quantile(Value[Value > 0 & !is.na(Value)], 0.25),
    Quart3 = quantile(Value[Value > 0 & !is.na(Value)], 0.75),
    IQR_A = (Quart3 - Quart1) * 1.5
  ) %>%
  group_by(chemical.id, endpoint, plate.id) %>%
  summarise(
    Hypo = sum(Value < 0, na.rm = T),
    Hyper = sum(Value > 0 & (Value >= (Quart1 - IQR_A) | Value <= (Quart3 + IQR_A)), na.rm = T),
    Total = length(Value),
    ThresholdAt0 = Total * 0.5,
    Remove = Hyper < ThresholdAt0 
  ) %>%
  dplyr::mutate(NewID = paste(chemical.id, plate.id, endpoint))

PlatesToRm <- AUC_0 %>% filter(Remove) %>% ungroup() %>% select(NewID) %>% unlist()

AUC_Final <- AUC %>%
  mutate(NewID = paste(chemical.id, plate.id, endpoint)) %>%
  filter(NewID %in% PlatesToRm == FALSE) %>%
  group_by(chemical.id, conc, endpoint) %>%
  summarise(
    Hypo = sum(Value < 0, na.rm = T),
    Total = length(Value),
    Response = Hypo / Total
  )

ggplot(AUC_Final, aes(x = conc, y = Response, color = endpoint)) + geom_line() + theme_bw() +
  ggtitle("AUC")
