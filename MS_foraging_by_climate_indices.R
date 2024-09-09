# -------------------------------------
# Script: foraging_by_climate_indices
# Author: Dr Natasha Gillies
# Purpose: Model foraging parameters according to broad-scale climate indices
# Date: 2022-03-11 
# -------------------------------------

# Load functions & packages ----------------------------------------------------

# Packages
packages <- c("DHARMa", "tidyverse", "glmmTMB", "brms", "ggeffects", "tidybayes", 
              "patchwork", "BNSP", "sjPlot", "ggpubr")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
#  if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

## For reference: boldness quantiles

# Sex     upperQuan lowerQuan midQuan
# <fct>       <dbl>     <dbl>   <dbl>
#   1 Females      1.80     -1.37 -0.0937
# 2 Males        1.40     -1.95 -0.376 


# Functions
## Calculate differences for different groups
predict_diffs <- function(dataset, group, group_val, val_name, val1, val2) {
  
  if (!is.na(group) && !is.na(group_val)) {
    mean_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$predicted) /
      subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$predicted * 100 
    
  } else {
    mean_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[val_name]] %in% val2)$predicted) /
      subset(dataset, dataset[[val_name]] %in% val1)$predicted * 100 
    
  }
  
  return(mean_pred)
}


## Predict changes for different groups for proportional variables
predict_diffs.prop <- function(dataset, group, group_val, val_name, val1, val2) {
  
  if (!is.na(group) && !is.na(group_val)) {
    mean_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$predicted) * 100 
    
  } else {
    mean_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[val_name]] %in% val2)$predicted) * 100 
    
  }
  
  return(mean_pred)
}



# VISUALISE DATA ===============================================================

femalegps <- read.csv("Data_inputs/WAAL_foraging_2010-2020_F.csv")
malegps <- read.csv("Data_inputs/WAAL_foraging_2010-2020_M.csv")

allgps <- rbind(femalegps, malegps)

## Path vs total distance ------------------------------------------------------

raw_dist_plot <- ggplot(allgps, aes(x = maxdistance.km, y = totalpathdistance.km, col = sex)) +
         geom_point() +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = c(0.065, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent")) +
  labs(x = "Maximum path distance (km)", y = "Total path distance (km)") 

  
log_dist_plot <- ggplot(allgps, aes(x = logmaxdistance.km, y = logtotalpathdistance.km, col = sex)) +
  geom_point() +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = "none") +
  labs(x = "Log maximum path distance (km)", y = "Log total path distance (km)") 


## Relative search time vs distance --------------------------------------------

search_dist_plot <- ggplot(allgps, aes(y = propARSvsTravel_hmm, x = totalpathdistance.km, col = sex)) +
  geom_point() +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = c(0.065, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent")) +
  labs(y = "Ratio search to travel", x = "Total path distance (km)") 

propsearch_dist_plot <- ggplot(allgps, aes(y = propSearch_hmm, x = totalpathdistance.km, col = sex)) +
  geom_point() +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = c(0.065, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent")) +
  labs(y = "Proportion of trip in search", x = "Total path distance (km)") 


## Distance against climate ----------------------------------------------------

soi_dist_plot <- ggplot(allgps, aes(y = maxdistance.km, x = SOIIndex, col = sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = c(0.065, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent")) +
  labs(y = "Max foraging distance (km)", x = "Southern Oscillation Index") 

sam_dist_plot <- ggplot(allgps, aes(y = maxdistance.km, x = SAMIndex, col = sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("dodgerblue", "orange"), name = "") +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = c(0.065, 0.95),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 12),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent")) +
  labs(y = "Max foraging distance (km)", x = "Southern Annular Mode")

png("Figures/sam_soi_distance.png", width = 8, height = 8, units = "in", res = 300)
ggpubr::ggarrange(soi_dist_plot, sam_dist_plot,
                  ncol = 1,
                  align = "hv")
dev.off()


# MODEL FITTING ================================================================

# femalegps <- read.csv("Data_inputs/WAAL_foraging_2010-2020_F.csv")
# malegps <- read.csv("Data_inputs/WAAL_foraging_2010-2020_M.csv")

## make permanent if all good
allgps <- read.csv("Data_inputs/SI_WAAL_foraging_2010-2020_additionalMetrics.csv")
femalegps <- subset(allgps, sex == "Female")
malegps <- subset(allgps, sex == "Male")

## Process the data -------------------------------------------------------------

## FEMALES ##

## Set id and year to factor
femalegps$id <- as.factor(femalegps$id)
femalegps$Year <- as.factor(femalegps$Year)

## Remove days that round to 0, NA personalities, and NA age
femalegps <- subset(femalegps, round(tripduration.days) > 0 & 
                      !is.na(boldness) & !is.na(Age))

# Isolate key variables
female_mod <- femalegps[,c("id", "Year", "Age", "DeploymentID", "sex", "boldness",
                           "totalpathdistance.km", "propARSvsTravel_hmm", "total_landings_hmm", "med_travel_dist",
                           "SAMIndex", "SOIIndex", 
                           "breeding_success", "attempted_breeding", "tripduration.days")]

## MALES ##

## Set id and year to factors
malegps$id <- as.factor(malegps$id)
malegps$Year <- as.factor(malegps$Year)

## Remove days that round to 0 and NA personalities
malegps <- subset(malegps, round(tripduration.days) > 0 & 
                    !is.na(boldness) & !is.na(Age))

# Isolate key variables
male_mod <- malegps[,c("id", "Year", "Age", "DeploymentID", "sex", "boldness",
                       "totalpathdistance.km", "propARSvsTravel_hmm", "total_landings_hmm", "med_travel_dist",
                       "SAMIndex", "SOIIndex", 
                       "breeding_success", "attempted_breeding", "tripduration.days")]


## Fit models using glmmTMB ----------------------------------------------------

#### Number of landings --------------------------------------------------------

# SAM #
# ~ F
f_landings.SAM <- glmmTMB(total_landings_hmm ~ 
                            SAMIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_landings.SAM))
summary(f_landings.SAM)

f_landings.SAM <- update(f_landings.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(f_landings.SAM)

tab_model(f_landings.SAM, show.stat = TRUE)

# ~ M
m_landings.SAM <- glmmTMB(total_landings_hmm ~ 
                            SAMIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = poisson(),
                          data = male_mod)

plot(simulateResiduals(m_landings.SAM))
summary(m_landings.SAM)

## Drop non-significant interactions
m_landings.SAM <- update(m_landings.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(m_landings.SAM)
tab_model(m_landings.SAM, show.stat = TRUE)


# SOI # 
# ~ F
f_landings.SOI <- glmmTMB(total_landings_hmm ~ 
                            SOIIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_landings.SOI))
summary(f_landings.SOI)

## Drop non-significant interactions
f_landings.SOI <- update(f_landings.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(f_landings.SOI)
tab_model(f_landings.SOI, show.stat = TRUE)

# ~ M
m_landings.SOI <- glmmTMB(total_landings_hmm ~ 
                            SOIIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = male_mod)

plot(simulateResiduals(m_landings.SOI))
summary(m_landings.SOI)
tab_model(m_landings.SOI, show.stat = TRUE)


#### Prop Search:Travel --------------------------------------------------------

# SAM #
# ~ F
f_searchTrav.SAM <- glmmTMB(propARSvsTravel_hmm ~ 
                              SAMIndex * boldness + Age +
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = female_mod)

plot(simulateResiduals(f_searchTrav.SAM))
summary(f_searchTrav.SAM)

## Drop non-significant interactions
f_searchTrav.SAM <- update(f_searchTrav.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(f_searchTrav.SAM)
tab_model(f_searchTrav.SAM, show.stat = TRUE)

# ~ M
m_searchTrav.SAM <- glmmTMB(propARSvsTravel_hmm ~ 
                              SAMIndex * boldness + Age +
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = male_mod)

plot(simulateResiduals(m_searchTrav.SAM))
summary(m_searchTrav.SAM)

## Drop non-significant interactions
m_searchTrav.SAM <- update(m_searchTrav.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(m_searchTrav.SAM)
tab_model(m_searchTrav.SAM, show.stat = TRUE)



# SOI # 
# ~ F
f_searchTrav.SOI <- glmmTMB(propARSvsTravel_hmm ~ 
                              SOIIndex * boldness + Age +
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = female_mod)

plot(simulateResiduals(f_searchTrav.SOI))
summary(f_searchTrav.SOI)

## Drop non-significant interactions
f_searchTrav.SOI <- update(f_searchTrav.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(f_searchTrav.SOI)
tab_model(f_searchTrav.SOI, show.stat = TRUE)

# ~ M
m_searchTrav.SOI <- glmmTMB(propARSvsTravel_hmm ~ 
                              SOIIndex * boldness + Age +
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = male_mod)

plot(simulateResiduals(m_searchTrav.SOI))
summary(m_searchTrav.SOI)

## Drop non-significant interactions
m_searchTrav.SOI <- update(m_searchTrav.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(m_searchTrav.SOI)
tab_model(m_searchTrav.SOI, show.stat = TRUE)


#### Path distance -------------------------------------------------------------

# SAM #
# ~ F
f_pathDist.SAM <- glmmTMB(totalpathdistance.km ~ 
                            SAMIndex * boldness + Age + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = female_mod)

plot(simulateResiduals(f_pathDist.SAM))
summary(f_pathDist.SAM)

## Drop non-significant interactions
f_pathDist.SAM <- update(f_pathDist.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(f_pathDist.SAM)

tab_model(f_pathDist.SAM, show.stat = T)

# ~ M 
m_pathDist.SAM <- glmmTMB(totalpathdistance.km ~ 
                            SAMIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = male_mod)

plot(simulateResiduals(m_pathDist.SAM))
summary(m_pathDist.SAM)

## Drop non-significant interactions
m_pathDist.SAM <- update(m_pathDist.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(m_pathDist.SAM)
tab_model(m_pathDist.SAM, show.stat = T)


# SOI # 
# ~ F
f_pathDist.SOI <- glmmTMB(totalpathdistance.km ~ 
                            SOIIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = Gamma(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_pathDist.SOI))
summary(f_pathDist.SOI)

f_pathDist.SOI <- update(f_pathDist.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(f_pathDist.SOI)

tab_model(f_pathDist.SOI, show.stat = T)


# ~ M 
m_pathDist.SOI <- glmmTMB(totalpathdistance.km ~ 
                            SOIIndex * boldness + Age +
                            (1|Year) + (1|id), 
                          family = Gamma(link = "log"),
                          data = male_mod)

plot(simulateResiduals(m_pathDist.SOI))
summary(m_pathDist.SOI)

## Drop non-significant interactions
m_pathDist.SOI <- update(m_pathDist.SOI, ~ SOIIndex + boldness + Age +(1|Year) + (1|id))
summary(m_pathDist.SOI)
tab_model(m_pathDist.SOI, show.stat = T)


### Median distance between patches --------------------------------------------

# SAM #
# ~ F
f_medDist.SAM <- glmmTMB(med_travel_dist ~ 
                           SAMIndex * boldness + Age +
                           (1|Year) + (1|id), 
                         family = Gamma(link = log),
                         data = female_mod)

plot(simulateResiduals(f_medDist.SAM))
summary(f_medDist.SAM)

f_medDist.SAM <- update(f_medDist.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(f_medDist.SAM)
tab_model(f_medDist.SAM, show.stat = T)

# ~ M
m_medDist.SAM <- glmmTMB(med_travel_dist ~ 
                           SAMIndex * boldness + Age + 
                           (1|id), 
                         family = Gamma(link = log),
                         data = male_mod)

plot(simulateResiduals(m_medDist.SAM))
summary(m_medDist.SAM)

m_medDist.SAM <- update(m_medDist.SAM, ~ SAMIndex + boldness + Age + (1|Year) + (1|id))
summary(m_medDist.SAM)
tab_model(m_medDist.SAM, show.stat = T)

# SOI #
# ~ F
f_medDist.SOI <- glmmTMB(med_travel_dist ~ 
                           SOIIndex * boldness + Age +
                           (1|Year) +(1|id), 
                         family = Gamma(link = log),
                         data = female_mod)

plot(simulateResiduals(f_medDist.SOI))
summary(f_medDist.SOI)

f_medDist.SOI <- update(f_medDist.SOI, ~ SOIIndex + boldness + Age + (1|Year) + (1|id))
summary(f_medDist.SOI)
tab_model(f_medDist.SOI, show.stat = T)

# ~ M
m_medDist.SOI <- glmmTMB(med_travel_dist ~ 
                           SOIIndex * boldness + Age +
                           (1|id), 
                         family = Gamma(link = log),
                         data = male_mod)

plot(simulateResiduals(m_medDist.SOI))
summary(m_medDist.SOI)

m_medDist.SOI <- update(m_medDist.SOI, ~ SOIIndex + boldness + Age + (1|id))
summary(m_medDist.SOI)
tab_model(m_medDist.SOI, show.stat = T)


# MODEL PLOTTING & PREDICTIONS #################################################

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

## SAM ---------------------------------------------------------------

### Landings -------------------------------------------------------------------
f_landings_sam.df <- data.frame(ggpredict(f_landings.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

# Overall reduction in landings with SOI
mean(predict_diffs(f_landings_sam.df, NA, NA, "SAM", 1, 0))

#### Build the plot
f_landings_sam.plot <- ggplot(data = f_landings_sam.df, aes(y = predicted, x = SAM)) +
  geom_point(aes(y = total_landings_hmm, x = SAMIndex), colour = "grey50", data = female_mod) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  labs(x = "Monthly Southern Annular Mode", y = "Number of landings per day") +
  ylim(0, 75) +
  theme_bw() +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


### Search : travel -----------------------------------------------------------
f_searchTrav_sam.df <- data.frame(ggpredict(f_searchTrav.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

# Bold vs shy reduction in distance with SAM
predict_diffs(f_searchTrav_sam.df,  NA, NA, "SAM", 1, 0)

# Build the plot 
f_searchTrav_sam.plot <- ggplot() +
  geom_point(aes(y = propARSvsTravel_hmm, x = SAMIndex), colour = "grey50", data = female_mod) +
  geom_line(data = f_searchTrav_sam.df, aes(y = predicted, x = SAM), linewidth = 1, col = female_col) +
  geom_ribbon(data = f_searchTrav_sam.df, aes(y = predicted, x = SAM,
              ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  labs(x = "Monthly Southern Annular Mode", y = "Ratio time spent in search : travel") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) 



### Median distance between patches -------------------------------------------------------------------
f_medDist_sam.df <- data.frame(ggpredict(f_medDist.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

# Overall reduction in landings with SOI
mean(predict_diffs(f_medDist_sam.df, NA, NA, "SAM", 0, 1))

#### Build the plot
f_medDist_sam.plot <- ggplot(data = f_medDist_sam.df, aes(y = predicted, x = SAM)) +
  geom_point(aes(y = med_travel_dist, x = SAMIndex), colour = "grey50", data = female_mod) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  ylim(25, 175) +
  labs(x = "Monthly Southern Annular Mode", y = "Distance between patches (median, km)") +
  theme_bw() +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


#### FIGURE 3: behaviours ~ SAM =================================================

png("Figures/FIG3_allBehaviour_SAM.png", width = 11, height = 10, units = "in", res = 300)
ggarrange(f_searchTrav_sam.plot, f_landings_sam.plot, 
          f_medDist_sam.plot,
          ncol = 2, nrow = 2,
          labels = "AUTO")
dev.off()

## SOI -------------------------------------------------------------------------

### Landings - F ---------------------------------------------------------------
f_landings_soi.df <- data.frame(ggpredict(f_landings.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

# Overall reduction in landings with SOI
mean(predict_diffs(f_landings_soi.df, NA, NA, "SOI", 0, 1))

#### Build the plot
f_landings_soi.plot <- ggplot(data = f_landings_soi.df, aes(y = predicted, x = SOI)) +
  geom_point(aes(y = total_landings_hmm, x = SOIIndex), colour = "grey50", data = female_mod) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  labs(x = "Monthly Southern Oscillation Index", y = "Number of landings per day") +
  ylim(0, 80) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


### Landings - M ---------------------------------------------------------------
quantile(malegps$boldness, c(0.1, 0.5, 0.9))

m_landings_soi.df <- data.frame(ggpredict(m_landings.SOI, 
                                          terms = c("SOIIndex", "boldness[-1.979724, -0.319167, 1.340847]"))) %>% 
  rename(SOI = x, boldness = group)

# Overall reduction in landings with SOI
mean(predict_diffs(m_landings_soi.df, NA, NA, "SOI", 0, 1))

# Bold vs shy reduction in landings with SOI
predict_diffs(m_landings_soi.df, "boldness", 3, "SOI", 0, 1)
predict_diffs(m_landings_soi.df, "boldness", 1, "SOI", 0, 1)

# Overall differences in landings
m_landings_soi.df$boldness <- as.character(as.numeric(m_landings_soi.df$boldness))
min_boldness.m <- min(m_landings_soi.df$boldness)
max_boldness.m <- max(m_landings_soi.df$boldness)

mean(subset(m_landings_soi.df, boldness == min_boldness.m)$predicted)
mean(subset(m_landings_soi.df, boldness == max_boldness.m)$predicted)

# Build the plot
m_landings_soi.df$boldness <- as.factor(m_landings_soi.df$boldness)

m_landings_soi.plot <- ggplot() +
  geom_point(aes(y = total_landings_hmm, x = SOIIndex), colour = "grey50", data = male_mod) +
  geom_ribbon(data = m_landings_soi.df, 
              aes(y = predicted, x = SOI, 
                  group = boldness, col = boldness, 
                  fill = boldness, ymin = conf.low, ymax = conf.high), 
              alpha = 0.15, col = NA) +
  geom_line(data = m_landings_soi.df, 
            aes(y = predicted, x = SOI, 
                group = boldness, col = boldness), linewidth = 1.5) +
  scale_colour_manual(values = c("#559de6", "#0465ba", "#042440"), labels = c("Shy", "Intermediate", "Bold"), 
                      name = "Boldness") +
  scale_fill_manual(values = c("#b9d3ed", "#76a2cf", "#4b647d"), labels = c("Shy", "Intermediate", "Bold"),
                    name = "Boldness") +
  labs(x = "Monthly Southern Oscillation Index", y = "Number of landings per day") +
  ylim(0, 80) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.85),
        axis.text.y =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        legend.key.size = unit(0.75, "cm"),
        legend.key = element_blank(), 
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent")) 



### Path distance - F ---------------------------------------------------------
f_pathDist_soi.df <- data.frame(ggpredict(f_pathDist.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

# Change in distnce with increasing SOI
predict_diffs(f_pathDist_soi.df, NA, NA, "SOI", 1, 0)

# Build the plot
f_pathDist_soi.plot <- ggplot(data = f_pathDist_soi.df, aes(y = predicted, x = SOI)) +
  geom_point(aes(y = totalpathdistance.km, x = SOIIndex), colour = "grey50", data = female_mod) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  labs(x = "Monthly Southern Oscillation Index", y = "Total trip distance (km)") +
  ylim(0, 25000) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x =  element_blank(),
        axis.title.x = element_blank())

#### Path distance - M ---------------------------------------------------------
m_pathDist_soi.df <- data.frame(ggpredict(m_pathDist.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

# Change in distnce with increasing SOI
predict_diffs(m_pathDist_soi.df, NA, NA, "SOI", 1, 0)

# Build the plot
m_pathDist_soi.plot <- ggplot(data = m_pathDist_soi.df, aes(y = predicted, x = SOI)) +
  geom_point(aes(y = totalpathdistance.km, x = SOIIndex), colour = "grey50", data = male_mod) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(linewidth = 1, col = male_col) +
  labs(x = "Monthly Southern Oscillation Index", y = "Total trip distance (km)") +
  ylim(0, 25000) + 
  theme_bw() +
  theme(axis.text.y =  element_blank(),
        axis.title.y = element_blank(),
        axis.text.x =  element_blank(),
        axis.title.x = element_blank())


##### FIGURE 4: behaviours ~ SOI ===============================================

png("Figures/FIG4_allBehaviour_SOI.png", width = 11, height = 10, units = "in", res = 300)
ggarrange(f_pathDist_soi.plot + theme(plot.margin = unit(c(0.2, 0, 0.5, 0), "cm")), 
          m_pathDist_soi.plot + theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.75), "cm")), 
          f_landings_soi.plot + theme(plot.margin = unit(c(0, 0, 0, 0.75), "cm")), 
          m_landings_soi.plot + theme(plot.margin = unit(c(0, 0.2, 0, 0.75), "cm")),
          ncol = 2, nrow = 2,
          labels = "AUTO")
dev.off()
