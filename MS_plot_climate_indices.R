# -------------------------------------
# Script: plot_climate_indices.R
# Author: Dr Jack Thorley; Dr Natasha Gillies
# Purpose: Visualise variation in SAM, SOI
# Notes: 
# Date: 2023-10-09
# -------------------------------------


# Load functions & packages ----------------------------------------------------

# Functions

# Packages
packages <- c("ggplot2", "tidyverse", "lubridate", "patchwork", "forecast") 

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
year.avg_SAM <- sam %>% 
  group_by(Year) %>% 
  summarise(meanSAM = mean(SAMIndex), 
            sdSAM = sd(SAMIndex)) %>% 
  data.frame()

# Slope of effect
summary(lm(meanSAM ~ Year, data = year.avg_SAM))

sam_plot <- ggplot(year.avg_SAM, aes(x = Year, y = meanSAM)) +
  geom_hline(yintercept = 0, col = "darkgrey", linetype = 2) +
  geom_ribbon(aes(ymin = meanSAM - sdSAM, ymax = meanSAM + sdSAM), alpha = 0.2, fill = "#7ED8BF") +
  geom_smooth(method = "lm", col = "#D55E00", se = T, fill = "#D55E00", alpha = 0.2, size = 1.1) +
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

# Slope of effect
summary(lm(meanSOI ~ Year, data = year.avg_SOI))


soi_plot <- ggplot(year.avg_SOI, aes(x = Year, y = meanSOI)) +
  geom_hline(yintercept = 0, col = "darkgrey", linetype = 2) +
  geom_ribbon(aes(ymin = meanSOI - sdSOI, ymax = meanSOI + sdSOI), alpha = 0.2, fill = "#7ED8BF") +
  geom_smooth(method = "lm", col = "#D55E00", se = T, fill = "#D55E00", alpha = 0.2, size = 1.1) +
  geom_path(size = 1.1) + 
  labs(x = "Year", title = "Southern Oscillation Index") + 
  scale_x_continuous(breaks = seq(1960, 2020, 5), labels = seq(1960, 2020, 5)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 




### FIGURE 2: climatic variation through time ----------------------------------

png("Figures/FIG2_climate_index_variation.png", width = 12, height = 6, units = "in", res = 300)
ggpubr::ggarrange(sam_plot, soi_plot,
                  ncol = 2, nrow = 1,
                  widths = c(1, 0.92))
dev.off()



# Correlation between SAM and SOI ----------------------------------------------

# Simple correlation
climate_merge <- merge(year.avg_SOI, year.avg_SAM, by = "Year")
correlation <- cor(climate_merge$meanSAM, climate_merge$meanSOI)

#  Detrended correlation
climate_merge$SAM_detrended <- diff(c(NA, climate_merge$meanSAM))
climate_merge$SOI_detrended <- diff(c(NA, climate_merge$meanSOI))

## Calculate the correlation of the detrended data
cor(climate_merge$SAM_detrended, climate_merge$SOI_detrended, use = "complete.obs")

ccf_result <- ccf(climate_merge$meanSAM, climate_merge$meanSOI, na.action = na.omit)
plot(ccf_result)

# ARIMA
fit_SAM <- auto.arima(climate_merge$meanSAM, seasonal = FALSE)
fit_SOI <- auto.arima(climate_merge$meanSOI, seasonal = FALSE)

## Check residuals for independence
checkresiduals(fit_SAM)
checkresiduals(fit_SOI)

## Calculate correlation of residuals
residuals_correlation <- cor(residuals(fit_SAM), residuals(fit_SOI))



# Correlation plot
corr_plot <- ggplot() +
  geom_hline(yintercept = 0, col = "darkgrey", linetype = 2) +
  # SOI #
  geom_smooth(data = year.avg_SOI, aes(x = Year, y = meanSOI),
              method = "lm", col = "#E69F00", se = T, fill = "#E69F00", alpha = 0.2, size = 1.1) +
  geom_path(data = year.avg_SOI, aes(x = Year, y = meanSOI, colour = "#E69F00"), alpha = 0.2, size = 1.1) + 
  # SAM #
  geom_smooth(data = year.avg_SAM, aes(x = Year, y = meanSAM),
              method = "lm", col = "#049E73", se = T, fill = "#049E73", alpha = 0.2, size = 1.1) +
  geom_path(data = year.avg_SAM, aes(x = Year, y = meanSAM, colour = "#049E73"), alpha = 0.2, size = 1.1) + 
  labs(x = "Year", y = "Average Index ± 1SD \n(Across years)") + 
  scale_x_continuous(breaks = seq(1960, 2020, 5), labels = seq(1960, 2020, 5)) +
  scale_colour_manual(name = 'Climate Index', guide = 'legend',
                      values = c('#049E73'='#049E73','#E69F00'='#E69F00'), labels = c('SOI','SAM')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.92, 0.92)) 



png("Figures/FIGS1_climate_index_correlation.png", width = 8, height = 7, units = "in", res = 300)
corr_plot
dev.off()

