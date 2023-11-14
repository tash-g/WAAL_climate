# -------------------------------------
# Script: plot_climate_indices.R
# Author: Dr Jack Thorley; edited by Dr Natasha Gillies
# Purpose: Visualise variation in SAM, SOI, IOD
# Notes: 
# Date: 2023-10-09
# -------------------------------------


# Load functions & packages ----------------------------------------------------

# Functions

# Packages
packages <- c("ggplot2", "tidyverse", "lubridate", "patchwork") 

## Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

## Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Set ggplot theme 
climate_theme <- theme_bw() + 
  theme(axis.title = element_text(size = 14, colour = "black"), 
        plot.title = element_text(size = 22, hjust = 0.5, colour = "royalblue4"),
        axis.text = element_text(size = 12, colour = "black", angle = 90), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 14, colour = "black"))


# Southern Annular Mode =======================================================

# Load dataset
load("Data_original/SAM_monthly.RData")

# Average SAM by year
year.avg_sam <- sam %>% 
  group_by(Year) %>% 
  summarise(meanSAM = mean(SAMIndex), 
            sdSAM = sd(SAMIndex)) %>% 
  data.frame()

sam_plot <- ggplot(year.avg_sam, aes(x = Year, y = meanSAM)) +
  geom_hline(yintercept = 0, col = "darkgrey", linetype = 2) +
  geom_ribbon(aes(ymin = meanSAM - sdSAM, ymax = meanSAM + sdSAM), alpha = 0.2, fill = "dodgerblue") +
  geom_smooth(method = "lm", col = "red", se = T, fill = "red", alpha = 0.2, size = 1.1) +
  geom_path(size = 1.1) + 
  labs(x = "Year", y = "Average Index ± 1SD \n(Across years)", title = "Southern Annular Mode") + 
  scale_x_continuous(breaks = seq(1960, 2020, 5), labels = seq(1960, 2020, 5)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



# Southern Oscillation Index ===================================================

# Load dataset
load("Data_original/SOI_monthly.RData")

# Average SOI by year
year.avg_SOI <- soi %>% 
  group_by(Year) %>% 
  summarise(meanSOI = mean(SOIIndex), 
            sdSOI = sd(SOIIndex)) %>% 
  data.frame()

soi_plot <- ggplot(year.avg_SOI, aes(x = Year, y = meanSOI)) +
  geom_hline(yintercept = 0, col = "darkgrey", linetype = 2) +
  geom_ribbon(aes(ymin = meanSOI - sdSOI, ymax = meanSOI + sdSOI), alpha = 0.2, fill = "dodgerblue") +
  geom_smooth(method = "lm", col = "red", se = T, fill = "red", alpha = 0.2, size = 1.1) +
  geom_path(size = 1.1) + 
  labs(x = "Year", title = "Southern Oscillation Index") + 
  scale_x_continuous(breaks = seq(1960, 2020, 5), labels = seq(1960, 2020, 5)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 




### FIGURE 2: climatic variation through time ------------------------------------

png("Figures/FIG2_climate_index_variation.png", width = 12, height = 6, units = "in", res = 300)
ggpubr::ggarrange(sam_plot, soi_plot,
                  ncol = 2, nrow = 1,
                  widths = c(1, 0.92))
dev.off()







# Indian Ocean Dipole ==========================================================

# Load dataset
load("Data_original/IOD_monthly.RData")

# Average IOD by year
year.avg_IOD <- iod %>% 
  group_by(Year) %>% 
  summarise(meanIOD = mean(IODIndex), 
            sdIOD = sd(IODIndex)) %>% 
  data.frame()

iod_plot <- ggplot(year.avg_IOD, aes(x = Year, y = meanIOD)) +
  geom_hline(yintercept = 0, col = "darkgrey", linetype = 2) +
  geom_ribbon(aes(ymin = meanIOD - sdIOD, ymax = meanIOD + sdIOD), alpha = 0.2, fill = "dodgerblue") +
  geom_smooth(method = "lm", col = "red", se = T, fill = "red", alpha = 0.2, size = 1.1) +
  geom_path(size = 1.1) + 
  labs(x = "Year", title = "IOD") + 
  scale_x_continuous(breaks = seq(1960, 2020, 5), labels = seq(1960, 2020, 5)) +
  climate_theme +
  theme(axis.text.y = element_blank()) 



