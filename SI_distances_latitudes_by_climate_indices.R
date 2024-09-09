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
                          "med_travel_dist", "min_latitude", "med_latitude",
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
                        "med_travel_dist", "min_latitude", "med_latitude",
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

## Latitude --------------------------------------------------------------------

# Make latitude absolute for easier plotting
female_mod$med_latitude <- abs(female_mod$med_latitude)
male_mod$med_latitude <- abs(male_mod$med_latitude)

# SAM #
# ~ F
f_medLatitude.SAM <- glmmTMB(med_latitude ~ 
                           SAMIndex * boldness + Age +
                           (1|Year) + (1|id), 
                         family = gaussian(link = "log"),
                         data = female_mod)

plot(simulateResiduals(f_medLatitude.SAM))
summary(f_medLatitude.SAM)

f_medLatitude.SAM <- update(f_medLatitude.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(f_medLatitude.SAM)

tab_model(f_medLatitude.SAM)

# ~ M
m_medLatitude.SAM <- glmmTMB(med_latitude ~ 
                           SAMIndex * boldness + Age +
                           (1|id), 
                         family = gaussian(link = log),
                         data = male_mod)

plot(simulateResiduals(m_medLatitude.SAM))
summary(m_medLatitude.SAM)

m_medLatitude.SAM <- update(m_medLatitude.SAM, ~ SAMIndex + boldness + Age +(1|Year) + (1|id))
summary(m_medLatitude.SAM)
tab_model(m_medLatitude.SAM, show.stat = T)

# SOI #
# ~ F
f_medLatitude.SOI <- glmmTMB(med_latitude ~ 
                           SOIIndex * boldness + Age +
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = female_mod)

plot(simulateResiduals(f_medLatitude.SOI))
summary(f_medLatitude.SOI)

f_medLatitude.SOI <- update(f_medLatitude.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(f_medLatitude.SOI, show.stat = T)
tab_model(f_medLatitude.SOI, show.stat = T)

# ~ M
m_medLatitude.SOI <- glmmTMB(med_latitude ~ 
                           SOIIndex * boldness + Age +
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = male_mod)

plot(simulateResiduals(m_medLatitude.SOI))
summary(m_medLatitude.SOI)

m_medLatitude.SOI <- update(m_medLatitude.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(m_medLatitude.SOI)
tab_model(m_medLatitude.SOI, show.stat = T)

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

## Median latitude ------------------------------------------------------------

# Females #
f_medLat_sam.df <- data.frame(ggpredict(f_medLatitude.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

# Build the plot 
f_medLat_sam.plot <- ggplot() +
  geom_line(data = f_medLat_sam.df, aes(y = (predicted*-1), x = SAM), linewidth = 1) +
  geom_ribbon(data = f_medLat_sam.df, aes(y = (predicted*-1), x = SAM,
                                          ymin = conf.low*-1, ymax = conf.high*-1), alpha = 0.2) +
  stat_pointinterval(data = subset(sam_range, sex == "Females"), aes(x = SAMIndex, y = -50), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Minimum latitude") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) 
