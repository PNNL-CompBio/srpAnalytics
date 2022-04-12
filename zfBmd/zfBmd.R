## Zebrafish Benchmark Database

## There are several defunct bmd packages. We are using this more recent one: https://github.com/DoseResponse/bmd
## devtools::install_github("DoseResponse/drcData")
## devtools::install_github("DoseResponse/drc")
## devtools::install_github("DoseResponse/bmd")

library(data.table)
library(dplyr)

# Read in pre-zf_bmd morpohology data
Morpho <- fread("~/Desktop/Git_Repos/srpAnalytics/zfBmd/test_files/7_PAH_zf_morphology_data_2020NOV11_tall_3756.csv")

# Add a function to calculate missing categories of ANY24, ANY120, TOT_MORT, and ALL_BUT_MORT 
AddEndpoints <- function(NewDose, Endpoints, Name) {
  NewDose %>%
    subset(End_Point %in% Endpoints) %>%
    mutate(Chemical_ID, End_Point = Name, Dose, Response) %>%
    group_by(Chemical_ID, End_Point, Dose) %>%
    summarise(Response = sum(Response, na.rm = TRUE))
}

# Get all endpoints we care about
AllEndpoints <- c("AXIS", "BRN_", "CRAN", "DP24", "EDEM", "LTRK", "MO24", "MORT", 
                  "MUSC", "NC__", "SKIN", "SM24", "TCHR")

# Group by chemical id, endpoint, and concentration. Sum counts with NA completely removed.
NewDose <- Morpho %>%
  subset(endpoint %in% AllEndpoints) %>%
  group_by(chemical.id, endpoint, conc) %>%
  summarise(Response = sum(value, na.rm = T) / length(value[!is.na(value)])) %>%
  rename(Chemical_ID = chemical.id, End_Point = endpoint, Dose = conc) %>%
  rbind(
    AddEndpoints(., c("MO24", "DP24", "SM24"), "ANY24"),
    AddEndpoints(., c("AXIS", "BRN_", "CRAN", "EDEM", "LTRK", "MORT", "MUSC", "NC__", "SKIN", "TCHR", "ANY24"), "ANY120"),
    AddEndpoints(., c("MO24", "MORT"), "TOT_MORT"),
    AddEndpoints(., c("AXIS", "BRN_", "CRAN", "DP24", "EDEM", 
                      "LTRK", "MUSC", "NC__", "SKIN", "SM24", "TCHR"), "ALL_BUT_MORT")
  )

MORT <- NewDose[NewDose$End_Point == "MORT",]

library(ggplot2)
library(drc)

ggplot(MORT, aes(x = Dose, y = Response)) + geom_point()

ll <- drm(formula = Response~Dose, data = MORT, type = "binomial", fct = LL.2()) 
ll %>% summary()
bmd(ll, 0.05, def = "excess")

glm(Response~Dose, family = "binomial", data = MORT) %>%
  summary()

glm(abs(Response)~abs(Dose), family = "Gamma", data = MORT) %>%
  summary()







