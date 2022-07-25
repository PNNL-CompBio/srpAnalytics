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

LPR_timepoints <- fread("~/Git_Repos/srpAnalytics/zfBmd/test_files/7_PAH_zf_LPR_data_2021JAN11_3756.csv")

LPR <- LPR_timepoints %>%
  mutate(
    Time = gsub("t", "", variable) %>% as.numeric() * 0.1,
    Concentration = factor(conc, levels = unique(conc)),
    Value = value,
    ID = as.factor(paste(plate.id, well))
  ) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::select(-c(variable, conc, value, plate.id, chemical.id, bottle.id, well)) %>%
  mutate(
    TimeGroup = lapply(Time, function(t) {
      if (t <= 4) {return("Remove")} else
      if (t <= 7) {return("Dark")} else
      if (t <= 10) {return("Light")} else
      if (t <= 13) {return("Dark")} else
      if (t <= 16) {return("Light")} else
      if (t <= 19) {return("Dark")} else
      if (t <= 22) {return("Light")} else {return("Remove")}}) %>% unlist()
  ) %>%
  filter(TimeGroup != "Remove") %>%
  mutate(
    TimeGroup = factor(TimeGroup, levels = c("Dark", "Light")),
  ) %>%
  group_by(Concentration, TimeGroup) %>%
  summarise(
    SumValue = sum(Value)
  )
  
# MOV1 is total dark movement - light movement 
LPR %>%
  ungroup() %>%
  group_by(Concentration) %>%
  summarise(
    `Dark - Light` = SumValue[match("Dark", TimeGroup)] - SumValue[match("Light", TimeGroup)]
  ) 




ggplot(toPlot, aes(x = TimeGroup, y = `Log Value`, fill = PlateID)) + geom_boxplot() +
  xlab("Time (min)") + theme_bw() + 
  facet_wrap(.~conc)



