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

load("Data_outputs/foraging_data_climate_analyses_M.RData")
load("Data_outputs/foraging_data_climate_analyses_F.RData")

## Process the data -------------------------------------------------------------

## FEMALES ##

## Set id and year to factor
femalegps$id <- as.factor(femalegps$id)
femalegps$Year <- as.factor(femalegps$Year)

## Remove days that round to 0 and NA personalities
femalegps <- subset(femalegps, round(tripduration.days) > 0 & 
                      !is.na(boldness))

## Get distance per day measure
femalegps$distPerDay <- femalegps$totalpathdistance.km/femalegps$tripduration.days

# Isolate key variables
female_mod <- femalegps[,c("id", "Year", "Age", "DeploymentID", "sex", "totalpathdistance.km",
                           "propARSvsTravel_hmm", "total_landings_hmm", "boldness",
                           "SAMIndex", "IODIndex", "SOIIndex", "breeding_success",
                           "attempted_breeding", "distPerDay", "tripduration.days")]

## MALES ##

## Set id and year to factors
malegps$id <- as.factor(malegps$id)
malegps$Year <- as.factor(malegps$Year)

## Remove days that round to 0 and NA personalities
malegps <- subset(malegps, round(tripduration.days) > 0 & 
                    !is.na(boldness))

## Get distance per day measure
malegps$distPerDay <- malegps$totalpathdistance.km/malegps$tripduration.days


# Isolate key variables
male_mod <- malegps[,c("id", "Year", "Age", "DeploymentID", "sex", "totalpathdistance.km",
                       "propARSvsTravel_hmm", "total_landings_hmm", "boldness",
                       "SAMIndex", "IODIndex", "SOIIndex", "breeding_success",
                       "attempted_breeding", "distPerDay", "tripduration.days")]


## Fit models using glmmTMB ----------------------------------------------------

#### Path distance -------------------------------------------------------------

# SAM #
# ~ F
f_pathDist.SAM <- glmmTMB(totalpathdistance.km ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = female_mod)

plot(simulateResiduals(f_pathDist.SAM))
summary(f_pathDist.SAM)


## Drop non-significant interactions
f_pathDist.SAM <- update(f_pathDist.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(f_pathDist.SAM)

tab_model(f_pathDist.SAM, show.stat = T)

# ~ M 
m_pathDist.SAM <- glmmTMB(totalpathdistance.km ~ 
                            SAMIndex * boldness +
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = male_mod)

plot(simulateResiduals(m_pathDist.SAM))
summary(m_pathDist.SAM)

## Drop non-significant interactions
m_pathDist.SAM <- update(m_pathDist.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(m_pathDist.SAM)
tab_model(m_pathDist.SAM, show.stat = T)


# SOI # 
# ~ F
f_pathDist.SOI <- glmmTMB(totalpathdistance.km ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_pathDist.SOI))
summary(f_pathDist.SOI)
tab_model(f_pathDist.SOI, show.stat = T)


# ~ M 
m_pathDist.SOI <- glmmTMB(totalpathdistance.km ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = "log"),
                          data = male_mod)

plot(simulateResiduals(m_pathDist.SOI))
summary(m_pathDist.SOI)

## Drop non-significant interactions
m_pathDist.SOI <- update(m_pathDist.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(m_pathDist.SOI)
tab_model(m_pathDist.SOI, show.stat = T)


#### Prop Search:Travel --------------------------------------------------------

# SAM #
# ~ F
f_searchTrav.SAM <- glmmTMB(propARSvsTravel_hmm ~ 
                              SAMIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = female_mod)

plot(simulateResiduals(f_searchTrav.SAM))
summary(f_searchTrav.SAM)

## Drop non-significant interactions
f_searchTrav.SAM <- update(f_searchTrav.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(f_searchTrav.SAM)
tab_model(f_searchTrav.SAM, show.stat = TRUE)

# ~ M
m_searchTrav.SAM <- glmmTMB(propARSvsTravel_hmm ~ 
                              SAMIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = male_mod)

plot(simulateResiduals(m_searchTrav.SAM))
summary(m_searchTrav.SAM)

## Drop non-significant interactions
m_searchTrav.SAM <- update(m_searchTrav.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(m_searchTrav.SAM)
tab_model(m_searchTrav.SAM, show.stat = TRUE)



# SOI # 
# ~ F
f_searchTrav.SOI <- glmmTMB(propARSvsTravel_hmm ~ 
                              SOIIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = female_mod)

plot(simulateResiduals(f_searchTrav.SOI))
summary(f_searchTrav.SOI)

## Drop non-significant interactions
f_searchTrav.SOI <- update(f_searchTrav.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(f_searchTrav.SOI)
tab_model(f_searchTrav.SOI, show.stat = TRUE)

# ~ M
m_searchTrav.SOI <- glmmTMB(propARSvsTravel_hmm ~ 
                              SOIIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = male_mod)

plot(simulateResiduals(m_searchTrav.SOI))
summary(m_searchTrav.SOI)

## Drop non-significant interactions
m_searchTrav.SOI <- update(m_searchTrav.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(m_searchTrav.SOI)
tab_model(m_searchTrav.SOI, show.stat = TRUE)


#### Number of landings --------------------------------------------------------

# SAM #
# ~ F
f_landings.SAM <- glmmTMB(total_landings_hmm ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_landings.SAM))
summary(f_landings.SAM)

f_landings.SAM <- update(f_landings.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(f_landings.SAM)

tab_model(f_landings.SAM, show.stat = TRUE)

# ~ M
m_landings.SAM <- glmmTMB(total_landings_hmm ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(),
                          data = male_mod)

plot(simulateResiduals(m_landings.SAM))
summary(m_landings.SAM)

## Drop non-significant interactions
m_landings.SAM <- update(m_landings.SAM, ~ SAMIndex + boldness + (1|Year) + (1|id))
summary(m_landings.SAM)
tab_model(m_landings.SAM, show.stat = TRUE)



# SOI # 
# ~ F
f_landings.SOI <- glmmTMB(total_landings_hmm ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_landings.SOI))
summary(f_landings.SOI)

## Drop non-significant interactions
f_landings.SOI <- update(f_landings.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(f_landings.SOI)
tab_model(f_landings.SOI, show.stat = TRUE)

# ~ M
m_landings.SOI <- glmmTMB(total_landings_hmm ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = male_mod)

plot(simulateResiduals(m_landings.SOI))
summary(m_landings.SOI)

## Drop non-significant interactions
m_landings.SOI <- update(m_landings.SOI, ~ SOIIndex + boldness + (1|Year) + (1|id))
summary(m_landings.SOI)
tab_model(m_landings.SOI, show.stat = TRUE)


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

## Path distance ---------------------------------------------------------------

# Females #
f_pathDist_soi.df <- data.frame(ggpredict(f_pathDist.SOI, terms = c("SOIIndex", "boldness[-1.275, 0.0278, 1.9]"))) %>% 
  rename(SOI = x, boldness = group) %>%
  mutate(boldness = as.numeric(boldness), boldness = case_when(boldness == min(boldness) ~ "Shy",
                              boldness == max(boldness) ~ "Bold",
                              TRUE ~ "Intermediate"))

f_pathDist_soi.df$boldness <- factor(f_pathDist_soi.df$boldness, levels = c("Shy", "Intermediate", "Bold"))

#### Make some predictions

# Bold vs shy reduction in distance with SOI
predict_diffs(f_pathDist_soi.df, "boldness", "Bold", "SOI", 0, 1)
predict_diffs(f_pathDist_soi.df, "boldness", "Shy", "SOI", 0, 1)


#### Build the plot 
f_pathDist_soi.plot <- ggplot() +
  geom_line(data = f_pathDist_soi.df, aes(y = predicted, x = SOI, col = boldness), linewidth = 1) +
  geom_ribbon(data = f_pathDist_soi.df, aes(y = predicted, x = SOI,
              ymin = conf.low, ymax = conf.high, fill = boldness), alpha = 0.2) +
  stat_pointinterval(data = subset(soi_range, sex == "Females"), aes(x = SOIIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Total trip distance (km)") +
  scale_y_continuous(limits = c(0,15000),
                      breaks = seq(0, 15000, 2500)) +
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
        axis.title.x = element_text(size = 14)) +
  theme(plot.margin = unit(c(0, 0, 0, 0.25), "cm"))


## Males ##
m_pathDist_soi.df <- data.frame(ggpredict(m_pathDist.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

#### Make some predictions
#  Reduction in distance for 1 point increase in SOI
predict_diffs(m_pathDist_soi.df, NA, NA, "SOI", -1, 0)


#### Build the plot 
m_pathDist_soi.plot <- ggplot() +
  geom_line(data = m_pathDist_soi.df, aes(y = predicted, x = SOI), col = male_col, linewidth = 1) +
  geom_ribbon(data = m_pathDist_soi.df, aes(y = predicted, x = SOI,
                                            ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = male_fill) +
  stat_pointinterval(data = subset(soi_range, sex == "Males"), aes(x = SOIIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Total trip distance (km)") +
  scale_y_continuous(limits = c(0,15000),
                     breaks = seq(0, 15000, 2500)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        plot.tag = element_text(size = 13)) +
  theme(plot.margin = unit(c(0, 0, 0, -0.05), "cm"))

#### Search : travel -----------------------------------------------------------

# Females #
f_searchTrav_sam.df <- data.frame(ggpredict(f_searchTrav.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

#### Make some predictions
predict_diffs.prop(f_searchTrav_sam.df, NA, NA, "SAM", -0.5, 0.5)

#### Build the plot
f_searchTrav_sam.plot <- ggplot(data = f_searchTrav_sam.df, aes(y = predicted, x = SAM)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  stat_pointinterval(data = subset(sam_range, sex == "Females"), aes(x = SAMIndex, y = 0.35), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Annular Mode", y = "Ratio time spent in search : travel") +
  scale_y_continuous(limits = c(0.35, 0.75),
                     breaks = seq(0.35, 0.75, 0.05)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) 

####

f_searchTrav_soi.df <- data.frame(ggpredict(f_searchTrav.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

#### Make some predictions
predict_diffs.prop(f_searchTrav_soi.df, NA, NA, "SOI", 1, 0)

#### Build the plot
f_searchTrav_soi.plot <- ggplot(data = f_searchTrav_soi.df, aes(y = predicted, x = SOI)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  stat_pointinterval(data = subset(soi_range, sex == "Females"), aes(x = SOIIndex, y = 0.35), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Ratio time spent in search : travel") +
  scale_y_continuous(limits = c(0.35, 0.75),
                     breaks = seq(0.35, 0.75, 0.05)) +
  theme_bw()  +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        plot.tag = element_text(size = 13)) +
  theme(plot.margin = unit(c(0, 0, 0, -0.05), "cm"))


#### Landings ------------------------------------------------------------------

# Females #
# SAM #
f_landings_sam.df <- data.frame(ggpredict(f_landings.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

#### Make some predictions
# Overall reduction in landings with SOI
mean(predict_diffs(f_landings_sam.df, NA, NA, "SAM", 0, 1))

#### Build the plot
f_landings_sam.plot <- ggplot(data = f_landings_sam.df, aes(y = predicted, x = SAM)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  stat_pointinterval(data = subset(sam_range, sex == "Females"), aes(x = SAMIndex, y = 0), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Annular Mode", y = "Number of landings per day") +
  scale_y_continuous(limits = c(0, 55),
                     breaks = seq(0, 55, 10)) +
  theme_bw() +
  theme(legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


# SOI #
f_landings_soi.df <- data.frame(ggpredict(f_landings.SOI, terms = c("SOIIndex", "boldness[-1.275, 0.0278, 1.9]"))) %>% 
  rename(SOI = x, boldness = group)

#### Make some predictions
# Overall reduction in landings with SOI
mean(predict_diffs(f_landings_soi.df, NA, NA, "SOI", 0, 1))

# Bold vs shy reduction in landings with SOI
f_landings_soi.df$boldness <- as.character(as.numeric(f_landings_soi.df$boldness))
min_boldness.f <- min(f_landings_soi.df$boldness)
max_boldness.f <- max(f_landings_soi.df$boldness)

mean(subset(f_landings_soi.df, boldness == min_boldness.f)$predicted)
mean(subset(f_landings_soi.df, boldness == max_boldness.f)$predicted)

#


#### Build the plot
f_landings_soi.df$boldness <- as.factor(f_landings_soi.df$boldness)

f_landings_soi.plot <- ggplot() +
  geom_ribbon(data = f_landings_soi.df, aes(y = predicted, x = SOI, 
                                            group = boldness, col = boldness, 
                                            fill = boldness, ymin = conf.low, ymax = conf.high), 
              alpha = 0.15, col = NA) +
  geom_line(data = f_landings_soi.df, aes(y = predicted, x = SOI, 
                                          group = boldness, col = boldness),
            linewidth = 1.5) +
  scale_colour_manual(values = c("#FFC20A", "#D9A601", "#584300"), labels = c("Shy", "Intermediate", "Bold"), 
                      name = "Boldness") +
  scale_fill_manual(values = c("#FFE48F", "#FFC20A", "#A37C01"), labels = c("Shy", "Intermediate", "Bold"),
                    name = "Boldness") +
  stat_pointinterval(data = subset(soi_range, sex == "Females"), aes(x = SOIIndex, y = 0), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Number of landings per day") +
  scale_y_continuous(limits = c(0, 110),
                     breaks = seq(0, 110, 20)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        plot.tag = element_text(size = 13)) +
  theme(plot.margin = unit(c(0, 0, 0, -0.05), "cm"))


####
# SOI #
m_landings_soi.df <- data.frame(ggpredict(m_landings.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x, boldness = group) 

#### Make some predictions
predict_diffs(m_landings_soi.df, NA, NA, "SOI", 0, 1)

m_landings_soi.df$boldness <- as.factor(m_landings_soi.df$boldness)

#### Build the plot
m_landings_soi.plot <- ggplot(data = m_landings_soi.df, aes(y = predicted, x = SOI)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(linewidth = 1, col = male_col) +
  stat_pointinterval(data = subset(soi_range, sex == "Males"), aes(x = SOIIndex, y = 0), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Monthly Southern Oscillation Index", y = "Number of landings per trip") +
  scale_y_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14)) +
theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


### FIGURE 3: all behaviours ~ SAM & SOI ---------------------------------------

png("Figures/FIG3_allBehaviour_SAM_SOI.png", width = 14, height = 24, units = "in", res = 300)
ggarrange(f_pathDist_soi.plot, m_pathDist_soi.plot, 
                  f_searchTrav_sam.plot, f_searchTrav_soi.plot,
                  f_landings_sam.plot, f_landings_soi.plot,
                  m_landings_soi.plot,
                  ncol = 2, nrow = 4,
                  labels = "AUTO",
                  hjust = c(0, -3.5),
                  #widths = c(1, 0.9),
                  align = "hv")
          
dev.off()

