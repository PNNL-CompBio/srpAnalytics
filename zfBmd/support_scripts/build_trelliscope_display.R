library(trelliscope)
library(tidyverse)
library(data.table)

setwd("~/Git_Repos/srpAnalytics/zfBmd/outputs/")

# Load datasets
bmds <- rbind(fread("20241204_BMDS_BC.csv"),
              fread("20241204_BMDS_LPR.csv")) %>%
  filter(Model != "")
dose <- rbind(fread("20241204_Dose_BC.csv"),
              fread("20241204_Dose_LPR.csv")) %>%
  filter(bmdrc.Endpoint.ID %in% bmds$bmdrc.Endpoint.ID) %>%
  mutate(Response = num.affected / num.nonna)
fits <- rbind(fread("20241204_Fits_BC.csv"),
              fread("20241204_Fits_LPR.csv")) %>%
  filter(bmdrc.Endpoint.ID %in% bmds$bmdrc.Endpoint.ID)

# Make the curve plot
curve_plot <- function(bmdrc.Endpoint.ID) {
  
  ggplot() +
      geom_line(data = dplyr::filter(fits, bmdrc.Endpoint.ID == {{bmdrc.Endpoint.ID}}), mapping = aes(x = X_vals, y = Y_vals)) +
      geom_point(data = dplyr::filter(dose, bmdrc.Endpoint.ID == {{bmdrc.Endpoint.ID}}), mapping = aes(x = Dose, y = Response)) +
      geom_segment(data = dplyr::filter(dose, bmdrc.Endpoint.ID == {{bmdrc.Endpoint.ID}}), 
                   mapping = aes(x = Dose, xend = Dose, y = CI_Lo, yend = CI_Hi)) +
      geom_point(data = dplyr::filter(bmds, bmdrc.Endpoint.ID == {{bmdrc.Endpoint.ID}}), mapping = aes(x = BMD10, y = 0.1), color = "red") +
      geom_point(data = dplyr::filter(bmds, bmdrc.Endpoint.ID == {{bmdrc.Endpoint.ID}}), mapping = aes(x = BMDL, y = 0.1), color = "blue") +
      xlab("Concentration (uM)") +
      ylab("Response (Proportion Affected)") +
      ylim(c(0, 1)) +
      ggtitle(bmdrc.Endpoint.ID) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

}

# Build off BMD column
bmds %>%
  mutate(plot = panel_lazy(curve_plot)) %>%
  as_trelliscope_df(name = "Fit Dichotomous Curves", description = "Visualizing Dichotomous Curve Fits", 
                    path = "~/Downloads/bmd_trelli", force_plot = T) %>%
  view_trelliscope()
