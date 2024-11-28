library(trelliscope)
library(tidyverse)
library(data.table)

setwd("~/Git_Repos/srpAnalytics/zfBmd/")

# Load datasets
bmds <- rbind(fread("new_BMDS_BC.csv"),
              fread("new_BMDS_LPR.csv")) %>% 
  rename(ids = bmdrc.Endpoint.ID)
dose <- rbind(fread("new_Dose_BC.csv"),
              fread("new_Dose_LPR.csv")) %>% 
  mutate(Response = num.affected/num.nonna)
fits <- rbind(fread("new_Fits_BC.csv"),
              fread("new_Fits_LPR.csv"))

# Add an ids column to fits
fits$ids <- paste(fits$Chemical_ID, fits$End_Point)

# Make the curve plot
curve_plot <- function(ids) {
  
  ggplot() +
      geom_line(data = dplyr::filter(fits, ids == {{ids}}), mapping = aes(x = X_vals, y = Y_vals)) +
      geom_point(data = dplyr::filter(dose, ids == {{ids}}), mapping = aes(x = Dose, y = Response)) +
      geom_segment(data = dplyr::filter(dose, ids == {{ids}}), 
                   mapping = aes(x = Dose, xend = Dose, y = CI_Lo, yend = CI_Hi)) +
      geom_point(data = dplyr::filter(bmds, ids == {{ids}}), mapping = aes(x = BMD10, y = 0.1), color = "red") +
      geom_point(data = dplyr::filter(bmds, ids == {{ids}}), mapping = aes(x = BMDL, y = 0.1), color = "blue") +
      xlab("Concentration (uM)") +
      ylab("Response (Proportion Affected)") +
      ylim(c(0, 1)) +
      ggtitle(ids) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

}

# Build off BMD column
bmds %>%
  mutate(plot = panel_lazy(curve_plot)) %>%
  as_trelliscope_df(name = "Fit Dichotomous Curves", description = "Visualizing Dichotomous Curve Fits", 
                    path = "~/Downloads/bmd_trelli", force_plot = T) %>%
  view_trelliscope()
