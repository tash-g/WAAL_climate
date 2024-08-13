# -------------------------------------
# Script: SI_distances_latitudes_by_climate_indices
# Author: Dr Natasha Gillies
# Purpose: Model foraging parameters according to broad-scale climate indices
# Date: 2024-08-12
# -------------------------------------

# Load functions & packages ----------------------------------------------------

# Packages
packages <- c("DHARMa", "tidyverse", "glmmTMB", "brms", "ggeffects", "tidybayes", 
              "patchwork", "sjPlot", "ggpubr", "magrittr")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
#  if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))


# Load data ---------------------------------------------------------------

allgps <- read.csv("Data_inputs/SI_WAAL_foraging_2010-2020_additionalMetrics.csv")

# Split by sex
allgps_F <- subset(allgps, sex == "Female")
allgps_M <- subset(allgps, sex == "Male")


## Process the data -------------------------------------------------------------

## FEMALES ##

## Set id and year to factor
allgps_F$id <- as.factor(allgps_F$id)
allgps_F$Year <- as.factor(allgps_F$Year)

## Remove days that round to 0 and NA personalities
allgps_F <- subset(allgps_F, round(tripduration.days) > 0 & 
                      !is.na(boldness))

# Isolate key variables
female_mod <- allgps_F[,c("id", "Year", "Age", "DeploymentID", "sex", "boldness",
                          "SAMIndex", "SOIIndex", 
                          "med_travel_dist", "min_lat",
                          "totalpathdistance.km", "total_landings_hmm")]

## MALES ##

## Set id and year to factors
allgps_M$id <- as.factor(allgps_M$id)
allgps_M$Year <- as.factor(allgps_M$Year)

## Remove days that round to 0 and NA personalities
allgps_M <- subset(allgps_M, round(tripduration.days) > 0 & 
                    !is.na(boldness))

## Get distance per day measure
allgps_M$distPerDay <- allgps_M$totalpathdistance.km/allgps_M$tripduration.days


# Isolate key variables
male_mod <- allgps_M[,c("id", "Year", "Age", "DeploymentID", "sex", "boldness",
                        "SAMIndex", "SOIIndex", 
                        "med_travel_dist", "min_lat",
                        "totalpathdistance.km", "total_landings_hmm")]



## Visualise relationships ----------------------------------------------------

# Travel distance with between-patch distance
dist_plot <- ggplot(allgps, aes(x = med_travel_dist, y = totalpathdistance.km, col = sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = c(0.1, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent")) +
  labs(x = "Median distance between patches (km)", y = "Total path distance (km)") +
  guides(color=guide_legend(override.aes=list(fill=NA)))

# Landings with between-patch distance
landings_plot <- ggplot(allgps, aes(y = total_landings_hmm, x = med_travel_dist, col = sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = "none") +
  labs(x = "Median distance between patches (km)", y = "Number of landings") +
  guides(color=guide_legend(override.aes=list(fill=NA)))

png("Figures/SI_distance_landings_patch_dist.png", width = 14, height = 7, units = "in", res = 300)
ggarrange(dist_plot, landings_plot, ncol = 2)
dev.off() 

# Fit models using glmmTMB ----------------------------------------------------

## Median distances ------------------------------------------------------------

# SAM #
# ~ F
f_medDist.SAM <- glmmTMB(med_travel_dist ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = female_mod)

plot(simulateResiduals(f_medDist.SAM))
summary(f_medDist.SAM)

f_medDist.SAM <- update(f_medDist.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(f_medDist.SAM)

# ~ M
m_medDist.SAM <- glmmTMB(med_travel_dist ~ 
                           SAMIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = male_mod)

plot(simulateResiduals(m_medDist.SAM))
summary(m_medDist.SAM)

m_medDist.SAM <- update(m_medDist.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(m_medDist.SAM)


# SOI #
# ~ F
f_medDist.SOI <- glmmTMB(med_travel_dist ~ 
                           SOIIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = female_mod)

plot(simulateResiduals(f_medDist.SOI))
summary(f_medDist.SOI)

f_medDist.SOI <- update(f_medDist.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(f_medDist.SOI)

# ~ M
m_medDist.SOI <- glmmTMB(med_travel_dist ~ 
                           SOIIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = male_mod)

plot(simulateResiduals(m_medDist.SOI))
summary(m_medDist.SOI)

m_medDist.SOI <- update(m_medDist.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(m_medDist.SOI)


## Latitude --------------------------------------------------------------------

# Make latitude absolute for easier plotting
female_mod$min_lat <- abs(female_mod$min_lat)
male_mod$min_lat <- abs(male_mod$min_lat)

# SAM #
# ~ F
f_minLat.SAM <- glmmTMB(min_lat ~ 
                           SAMIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = gaussian(link = "log"),
                         data = female_mod)

plot(simulateResiduals(f_minLat.SAM))
summary(f_minLat.SAM)

f_minLat.SAM <- update(f_minLat.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(f_minLat.SAM)

# ~ M
m_minLat.SAM <- glmmTMB(min_lat ~ 
                           SAMIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = gaussian(link = log),
                         data = male_mod)

plot(simulateResiduals(m_minLat.SAM))
summary(m_minLat.SAM)

m_minLat.SAM <- update(m_minLat.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(m_minLat.SAM)


# SOI #
# ~ F
f_minLat.SOI <- glmmTMB(min_lat ~ 
                           SOIIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = female_mod)

plot(simulateResiduals(f_minLat.SOI))
summary(f_minLat.SOI)

# ~ M
m_minLat.SOI <- glmmTMB(min_lat ~ 
                           SOIIndex * boldness + 
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = male_mod)

plot(simulateResiduals(m_minLat.SOI))
summary(m_minLat.SOI)

m_minLat.SOI <- update(m_minLat.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(m_minLat.SOI)



# Plot models ------------------------------------------------------------------

female_col <- "#FFC20A"
female_fill <- "#f7dc8b"
male_col <- "#0C7BDC"
male_fill <- "#91c1eb"

# Get climate index ranges
soi_range <- rbind(male_mod, female_mod) %>% 
  dplyr::select(c(SOIIndex, sex)) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

sam_range <- rbind(male_mod, female_mod) %>% 
  dplyr::select(c(SAMIndex, sex)) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

# Boldness quantiles
quantile(femalegps$boldness, c(0.1, 0.5, 0.9))
# 10%         50%         90% 
# -1.27540065  0.02783998  1.90034646 
quantile(malegps$boldness, c(0.1, 0.5, 0.9))
# 10%       50%       90% 
# -1.979724 -0.319167  1.340847 

## Median distance between patches ---------------------------------------------

## Females ##
f_medDist_sam.df <- data.frame(ggpredict(f_medDist.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

# Build the plot 
f_medDist_sam.plot <- ggplot() +
  geom_line(data = f_medDist_sam.df, aes(y = predicted, x = SAM), col = female_col, linewidth = 1) +
  geom_ribbon(data = f_medDist_sam.df, aes(y = predicted, x = SAM,
                                            ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = female_fill) +
  stat_pointinterval(data = subset(sam_range, sex == "Females"), aes(x = SAMIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Annular Mode", y = "Median distance between patches (km)") +
  scale_y_continuous(limits = c(0,80)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.tag = element_text(size = 13)) 



## Males ##
m_medDist_sam.df <- data.frame(ggpredict(m_medDist.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

# Build the plot 
m_medDist_sam.plot <- ggplot() +
  geom_line(data = m_medDist_sam.df, aes(y = predicted, x = SAM), col = male_col, linewidth = 1,
            linetype = "dashed") +
  geom_ribbon(data = m_medDist_sam.df, aes(y = predicted, x = SAM,
                                           ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = male_fill) +
  stat_pointinterval(data = subset(sam_range, sex == "Males"), aes(x = SAMIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Annular Mode", y = "Total trip distance (km)") +
  scale_y_continuous(limits = c(0,80)) + 
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        plot.tag = element_text(size = 13)) 


ggarrange(f_medDist_sam.plot, m_medDist_sam.plot, ncol = 2)



## Maximum latitude ------------------------------------------------------------

# Females #
f_minLat_soi.df <- data.frame(ggpredict(f_minLat.SOI, terms = c("SOIIndex", "boldness[-1.275, 0.0278, 1.9]"))) %>% 
  rename(SOI = x, boldness = group) %>%
  mutate(boldness = as.numeric(boldness), boldness = case_when(boldness == min(boldness) ~ "Shy",
                                                               boldness == max(boldness) ~ "Bold",
                                                               TRUE ~ "Intermediate"))

f_minLat_soi.df$boldness <- factor(f_minLat_soi.df$boldness, levels = c("Shy", "Intermediate", "Bold"))

# Build the plot 
f_minLat_soi.plot <- ggplot() +
  geom_line(data = f_minLat_soi.df, aes(y = (predicted*-1), x = SOI, col = boldness), linewidth = 1) +
  geom_ribbon(data = f_minLat_soi.df, aes(y = (predicted*-1), x = SOI,
                                          ymin = conf.low*-1, ymax = conf.high*-1, 
                                          fill = boldness), alpha = 0.2) +
  stat_pointinterval(data = subset(soi_range, sex == "Females"), aes(x = SOIIndex, y = -50), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Minimum latitude") +
  scale_colour_manual(values = c("#FFE48F", "#D9A601", "#584300"), labels = c("Shy", "Intermediate", "Bold"),
                      name = "Boldness") +
  scale_fill_manual(values = c("#FFE48F", "#FFC20A", "#A37C01"), labels = c("Shy", "Intermediate", "Bold"),
                    name = "Boldness") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) 
