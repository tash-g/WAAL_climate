# -------------------------------------
# Script: foraging_by_climate_indices
# Author: Dr Jack Thorley; adapted by Dr Natasha Gillies
# Purpose: Model foraging parameters according to broad-scale climate indices
# Notes: Adapted from original script: Foraging Trip Variation_brms_new.R
# Date: 2022-03-11 (modification)
# -------------------------------------

# Load functions & packages ----------------------------------------------------

# Packages
packages <- c("DHARMa", "tidyverse", "glmmTMB", "brms", "ggeffects", "tidybayes", 
              "patchwork", "BNSP", "sjPlot")

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
## Calculate percentage changes for different groups
predict_diffs <- function(dataset, group, group_val, val_name, val1, val2) {
  
  if (!is.na(group) && !is.na(group_val)) {
    mean_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$predicted) /
      subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$predicted * 100 
    
    lower_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$conf.low -
                     subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$conf.low) /
      subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$conf.low * 100
    
    upper_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$conf.high -
                     subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$conf.high) /
      subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$conf.high * 100
  } else {
    mean_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[val_name]] %in% val2)$predicted) /
      subset(dataset, dataset[[val_name]] %in% val1)$predicted * 100 
    
    lower_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$conf.low -
                     subset(dataset, dataset[[val_name]] %in% val2)$conf.low) /
      subset(dataset, dataset[[val_name]] %in% val1)$conf.low * 100
    
    upper_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$conf.high -
                     subset(dataset, dataset[[val_name]] %in% val2)$conf.high) /
      subset(dataset, dataset[[val_name]] %in% val1)$conf.high * 100
  }
  
  df_pred <- cbind(mean_pred, lower_pred, upper_pred)
  return(df_pred)
}


predict_diffs.prop <- function(dataset, group, group_val, val_name, val1, val2) {
  
  if (!is.na(group) && !is.na(group_val)) {
    mean_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$predicted) * 100 
    
    lower_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$conf.low -
                     subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$conf.low) * 100
    
    upper_pred <- (subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val1)$conf.high -
                     subset(dataset, dataset[[group]] == group_val & dataset[[val_name]] %in% val2)$conf.high) * 100
  } else {
    mean_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$predicted -
                    subset(dataset, dataset[[val_name]] %in% val2)$predicted) * 100 
    
    lower_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$conf.low -
                     subset(dataset, dataset[[val_name]] %in% val2)$conf.low) * 100
    
    upper_pred <- (subset(dataset, dataset[[val_name]] %in% val1)$conf.high -
                     subset(dataset, dataset[[val_name]] %in% val2)$conf.high) * 100
  }
  
  df_pred <- cbind(mean_pred, lower_pred, upper_pred)
  return(df_pred)
}


# VISUALISATION ================================================================

load("Data_outputs/foraging_data_climate_analyses_M.RData")
load("Data_outputs/foraging_data_climate_analyses_F.RData")

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

png("Figures/distances.png", width = 8, height = 8, units = "in", res = 300)
ggpubr::ggarrange(raw_dist_plot, log_dist_plot,
                  ncol = 1,
                  align = "hv")
dev.off()


# ANALYSIS =====================================================================

load("Data_outputs/foraging_data_climate_analyses_M.RData")
load("Data_outputs/foraging_data_climate_analyses_F.RData")

## Process the data -------------------------------------------------------------

## FEMALES ##

# Set id and year to factor
femalegps$id <- as.factor(femalegps$id)
femalegps$Year <- as.factor(femalegps$Year)

# Remove days that round to 0 and NA personalities
femalegps <- subset(femalegps, round(tripduration.days) > 0 & 
                      !is.na(boldness))

# Isolate key variables
female_mod <- femalegps[,c("id", "Year", "Age", "DeploymentID", "sex", "totalpathdistance.km",
                           "propARSvsTravel_hmm", "total_landings_hmm", "boldness",
                           "SAMIndex", "IODIndex", "SOIIndex", "breeding_success",
                           "attempted_breeding")]

## MALES ##

# Set id and year to factors
malegps$id <- as.factor(malegps$id)
malegps$Year <- as.factor(malegps$Year)

# Remove days that round to 0 and NA personalities
malegps <- subset(malegps, round(tripduration.days) > 0 & 
                    !is.na(boldness))

male_mod <- malegps[,c("id", "Year", "Age", "DeploymentID", "sex", "totalpathdistance.km",
                       "propARSvsTravel_hmm", "total_landings_hmm", "boldness",
                       "SAMIndex", "IODIndex", "SOIIndex", "breeding_success",
                       "attempted_breeding")]

## Visualise variables ---------------------------------------------------------

# Histograms
par(mfrow = c(3,2))
hist(female_mod$total_landings_hmm, main = "Female landings"); 
hist(male_mod$total_landings_hmm, main = "Male landings");
hist(female_mod$travSearch_ratio, main = "Female travel:search"); 
hist(male_mod$travSearch_ratio, main = "Male travel:search"); 
hist(female_mod$logmaxdistance.km, main = "Female max distance"); 
hist(male_mod$logmaxdistance.km, main = "Male max distance")

# Correlations
car::scatterplotMatrix(~total_landings_hmm + travSearch_ratio + logmaxdistance.km, 
                       data = female_mod)

car::scatterplotMatrix(~total_landings_hmm + travSearch_ratio + logmaxdistance.km, 
                       data = male_mod)


## Fit models using glmmTMB ----------------------------------------------------

#### Path distance -------------------------------------------------------------

# SAM #
f_pathDist.SAM <- glmmTMB(totalpathdistance.km ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = female_mod)

plot(simulateResiduals(f_pathDist.SAM))
summary(f_pathDist.SAM)
tab_model(f_pathDist.SAM, show.stat = T)

m_pathDist.SAM <- glmmTMB(totalpathdistance.km ~ 
                            SAMIndex * boldness +
                            (1|Year) + (1|id), 
                          family = Gamma(link = log),
                          data = male_mod)

plot(simulateResiduals(m_pathDist.SAM))
summary(m_pathDist.SAM)
tab_model(m_pathDist.SAM, show.stat = T)


# SOI # 
f_pathDist.SOI <- glmmTMB(totalpathdistance.km ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_pathDist.SOI))
summary(f_pathDist.SOI)
tab_model(f_pathDist.SOI, show.stat = T)


m_pathDist.SOI <- glmmTMB(totalpathdistance.km ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = Gamma(link = "log"),
                          data = male_mod)

plot(simulateResiduals(m_pathDist.SOI))
summary(m_pathDist.SOI)
tab_model(m_pathDist.SOI, show.stat = T)


#### Prop Search:Travel --------------------------------------------------------

# SAM #
f_searchTrav.SAM <- glmmTMB(propARSvsTravel_hmm ~ 
                              SAMIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = female_mod)

plot(simulateResiduals(f_searchTrav.SAM))
summary(f_searchTrav.SAM)
tab_model(f_searchTrav.SAM, show.stat = TRUE)


m_searchTrav.SAM <- glmmTMB(propARSvsTravel_hmm ~ 
                              SAMIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = male_mod)

plot(simulateResiduals(m_searchTrav.SAM))
summary(m_searchTrav.SAM)
tab_model(m_searchTrav.SAM, show.stat = TRUE)



# SOI # 
f_searchTrav.SOI <- glmmTMB(propARSvsTravel_hmm ~ 
                              SOIIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = female_mod)

plot(simulateResiduals(f_searchTrav.SOI))
summary(f_searchTrav.SOI)
tab_model(f_searchTrav.SOI, show.stat = TRUE)


m_searchTrav.SOI <- glmmTMB(propARSvsTravel_hmm ~ 
                              SOIIndex * boldness + 
                              (1|Year) + (1|id), 
                            family = beta_family,
                            data = male_mod)

plot(simulateResiduals(m_searchTrav.SOI))
summary(m_searchTrav.SOI)
tab_model(m_searchTrav.SOI, show.stat = TRUE)



#### Number of landings --------------------------------------------------------

# SAM #
f_landings.SAM <- glmmTMB(total_landings_hmm ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_landings.SAM))
summary(f_landings.SAM)
tab_model(f_landings.SAM, show.stat = TRUE)



m_landings.SAM <- glmmTMB(total_landings_hmm ~ 
                            SAMIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(),
                          data = male_mod)

plot(simulateResiduals(m_landings.SAM))
summary(m_landings.SAM)
tab_model(m_landings.SAM, show.stat = TRUE)



# SOI # 
f_landings.SOI <- glmmTMB(total_landings_hmm ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = female_mod)

plot(simulateResiduals(f_landings.SOI))
summary(f_landings.SOI)
tab_model(f_landings.SOI, show.stat = TRUE)


m_landings.SOI <- glmmTMB(total_landings_hmm ~ 
                            SOIIndex * boldness + 
                            (1|Year) + (1|id), 
                          family = poisson(link = "log"),
                          data = male_mod)

plot(simulateResiduals(m_landings.SOI))
summary(m_landings.SOI)
tab_model(m_landings.SOI, show.stat = TRUE)


## Plot the model predictions -------------------------------------------------

female_col <- "#FFC20A"
female_fill <- "#f7dc8b"
male_col <- "#0C7BDC"
male_fill <- "#91c1eb"

# Get climate index ranges
soi_range <- rbind(male_mod, female_mod) %>% 
  dplyr::select(c(SOIIndex, sex)) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

iod_range <- rbind(male_mod, female_mod) %>% 
  dplyr::select(c(IODIndex, sex)) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

sam_range <- rbind(male_mod, female_mod) %>% 
  dplyr::select(c(SAMIndex, sex)) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))


#### Path distance -------------------------------------------------------------

# Females #
f_pathDist_soi.df <- data.frame(ggpredict(f_pathDist.SOI, terms = c("SOIIndex", "boldness[-1.37, 0.094, 1.80]"))) %>% 
  rename(SOI = x, boldness = group) %>%
  mutate(boldness = as.numeric(boldness), boldness = case_when(boldness == min(boldness) ~ "Shy",
                              boldness == max(boldness) ~ "Bold",
                              TRUE ~ "Intermediate"))

f_pathDist_soi.df$boldness <- factor(f_pathDist_soi.df$boldness, levels = c("Shy", "Intermediate", "Bold"))

#### Make some predictions

# Bold vs shy reduction in distance with SOI
predict_diffs(f_pathDist_soi.df, "boldness", "Bold", "SOI", -2, 2)
predict_diffs(f_pathDist_soi.df, "boldness", "Shy", "SOI", -2, 2)


#### Build the plot 
f_pathDist_soi.plot <- ggplot() +
  geom_line(data = f_pathDist_soi.df, aes(y = predicted, x = SOI, col = boldness), linewidth = 1) +
  geom_ribbon(data = f_pathDist_soi.df, aes(y = predicted, x = SOI,
              ymin = conf.low, ymax = conf.high, fill = boldness), alpha = 0.2) +
  stat_pointinterval(data = subset(soi_range, sex == "Females"), aes(x = SOIIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Oscillation Index (Breeding)", y = "Total trip distance (km)", title = "Females") +
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
        plot.title = element_text(face = "bold"))


## Males ##
m_pathDist_soi.df <- data.frame(ggpredict(m_pathDist.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

#### Make some predictions
# Reduction in distance for 1 point increase in SOI
predict_diffs(m_pathDist_soi.df, NA, NA, "SOI", -2, 2)

#### Build the plot 
m_pathDist_soi.plot <- ggplot(data = m_pathDist_soi.df, aes(y = predicted, x = SOI)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(linewidth = 1, col = male_col) +
  stat_pointinterval(data = subset(soi_range, sex == "Males"), aes(x = SOIIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Oscillation Index (Breeding)", title = "Males") +
  scale_y_continuous(limits = c(0,15000),
                     breaks = seq(0, 15000, 2500)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"))



### FIGURE 5: path distance ~ SOI ------------------------------------------------

png("Figures/FIG5_pathDist_by_SOI.png", width = 12, height = 6, units = "in", res = 300)
ggpubr::ggarrange(f_pathDist_soi.plot, m_pathDist_soi.plot, 
                  ncol = 2, nrow = 1,
                  widths = c(1, 0.92))
dev.off()

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
  labs(x = "Southern Annular Mode (Breeding)", y = "Ratio time spent in search : travel", title = "(a) Females") +
  scale_y_continuous(limits = c(0.35, 0.75),
                     breaks = seq(0.35, 0.75, 0.05)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

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
  labs(x = "Southern Oscillation Index (Breeding)", y = "Ratio time spent in search : travel", title = " (c) Females") +
  scale_y_continuous(limits = c(0.35, 0.75),
                     breaks = seq(0.35, 0.75, 0.05)) +
  theme_bw()  +
  theme(plot.title = element_text(face = "bold"))


# Males #
m_searchTrav_sam.df <- data.frame(ggpredict(m_searchTrav.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

m_searchTrav_sam.plot <- ggplot(data = m_searchTrav_sam.df, aes(y = predicted, x = SAM)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(linewidth = 1, col = male_col) +
  stat_pointinterval(data = subset(sam_range, sex == "Males"), aes(x = SAMIndex, y = 0.35), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Annular Mode (Breeding)", title = "(b) Males (n.s.)") +
  scale_y_continuous(limits = c(0.35, 0.75),
                     breaks = seq(0.35, 0.75, 0.05)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"))


####

m_searchTrav_soi.df <- data.frame(ggpredict(m_searchTrav.SOI, terms = "SOIIndex")) %>% 
  rename(SOI = x) 

m_searchTrav_soi.plot <- ggplot(data = m_searchTrav_soi.df, aes(y = predicted, x = SOI)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(linewidth = 1, col = male_col) +
  stat_pointinterval(data = subset(soi_range, sex == "Males"), aes(x = SOIIndex, y = 0.35), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Oscillation Index (Breeding)", title = "(d) Males (n.s.)") +
  scale_y_continuous(limits = c(0.35, 0.75),
                     breaks = seq(0.35, 0.75, 0.05)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"))


### FIGURE 6: travel : search ~ SOI, SAM, IOD ------------------------------------

png("Figures/FIG6_searchTrav_by_SOISAM.png", width = 12, height = 12, units = "in", res = 300)
ggpubr::ggarrange(f_searchTrav_sam.plot, m_searchTrav_sam.plot, 
                  f_searchTrav_soi.plot, m_searchTrav_soi.plot,
                  ncol = 2, nrow = 2,
                  widths = c(1, 0.92))
dev.off()


#### Landings ------------------------------------------------------------------

# Females #
f_landings_sam.df <- data.frame(ggpredict(f_landings.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

#### Make some predictions
predict_diffs(f_landings_sam.df, NA, NA, "SAM", 0.5, 1.5)

#### Build the plot

f_landings_sam.plot <- ggplot(data = f_landings_sam.df, aes(y = predicted, x = SAM)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(linewidth = 1, col = female_col) +
  stat_pointinterval(data = subset(sam_range, sex == "Females"), aes(x = SAMIndex, y = 0), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Annular Mode (Breeding)", y = "Number of landings per trip", title = "(a) Females") +
  scale_y_continuous(limits = c(0, 35),
                     breaks = seq(0, 35, 10)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

####

f_landings_soi.df <- data.frame(ggpredict(f_landings.SOI, terms = c("SOIIndex", "boldness[-1.37, 0.094, 1.80]"))) %>% 
  rename(SOI = x, boldness = group) 

#### Make some predictions
# Overall reduction in landings with SOI
predict_diffs(f_landings_soi.df, NA, NA, "SOI", 0, 1)

# Bold vs shy reduction in landings with SOI
f_landings_soi.df$boldness <- as.character(as.numeric(f_landings_soi.df$boldness))
min_boldness.f <- min(f_landings_soi.df$boldness)
max_boldness.f <- max(f_landings_soi.df$boldness)

predict_diffs(f_landings_soi.df, "boldness", min_boldness.f, "SOI", 0, 1)
predict_diffs(f_landings_soi.df, "boldness", max_boldness.f, "SOI", 0, 1)

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
  labs(x = "Southern Oscillation Index (Breeding)", y = "Number of landings per trip", title = "(c) Females") +
  scale_y_continuous(limits = c(0, 70),
                     breaks = seq(0, 70, 10)) +
  theme_bw() +
  theme(legend.position = c(0.85,0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(face = "bold"))


# Males #
m_landings_sam.df <- data.frame(ggpredict(m_landings.SAM, terms = "SAMIndex")) %>% 
  rename(SAM = x) 

m_landings_sam.plot <- ggplot(data = m_landings_sam.df, aes(y = predicted, x = SAM)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(linewidth = 1, col = male_col) +
  stat_pointinterval(data = subset(sam_range, sex == "Males"), aes(x = SAMIndex, y = 0), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Annular Mode (Breeding)", y = "Number of landings per trip", title = "(b) Males (N.S.)") +
  scale_y_continuous(limits = c(0, 35),
                     breaks = seq(0, 35, 10)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"))

####

m_landings_soi.df <- data.frame(ggpredict(m_landings.SOI, terms = c("SOIIndex", "boldness[01.95, -0.38, 1.40"))) %>% 
  rename(SOI = x, boldness = group) 

#### Make some predictions
# Overall reduction in landings with SOI
predict_diffs(m_landings_soi.df, NA, NA, "SOI", 0, 1)

# Bold vs shy reduction in landings with SOI
m_landings_soi.df$boldness <- as.character(as.numeric(m_landings_soi.df$boldness))
min_boldness.m <- min(m_landings_soi.df$boldness)
max_boldness.m <- max(m_landings_soi.df$boldness)

predict_diffs(m_landings_soi.df, "boldness", min_boldness.m, "SOI", 0, 1)
predict_diffs(m_landings_soi.df, "boldness", max_boldness.m, "SOI", 0, 1)

( subset(m_landings_soi.df, boldness == max(boldness) & SOI == 0)$predicted -
    subset(m_landings_soi.df, boldness == max(boldness) & SOI == 1)$predicted  ) /
  subset(m_landings_soi.df, boldness == max(boldness) & SOI == 0)$predicted  * 100

( subset(m_landings_soi.df, boldness == min(boldness) & SOI == 0)$predicted -
    subset(m_landings_soi.df, boldness == min(boldness) & SOI == 1)$predicted  ) /
  subset(m_landings_soi.df, boldness == min(boldness) & SOI == 0)$predicted  * 100

m_landings_soi.df$boldness <- as.factor(m_landings_soi.df$boldness)

#### Build the plot
m_landings_soi.plot <- ggplot() +
  geom_ribbon(data = m_landings_soi.df, aes(y = predicted, x = SOI, 
                                            group = boldness, col = boldness, fill = boldness, ymin = conf.low, ymax = conf.high), 
              alpha = 0.25, col = NA) +
  geom_line(data = m_landings_soi.df, aes(y = predicted, x = SOI, 
                                          group = boldness, col = boldness),
            linewidth = 1.5) +
  scale_colour_manual(values = c("#2797FD", "#0262BA", "#01203D"), labels = c("Shy", "Intermediate", "Bold"), 
                      name = "Boldness") +
  scale_fill_manual(values = c("#A0D1FE", "#0276DB", "#013E75"), labels = c("Shy", "Intermediate", "Bold"),
                    name = "Boldness") +
  stat_pointinterval(data = subset(soi_range, sex == "Males"), aes(x = SOIIndex, y = 0), 
                     point_size = 3.5, colour = "darkgrey") +
  labs(x = "Southern Oscillation Index (Breeding)", y = "Number of landings per trip", title = "(d) Males") +
  scale_y_continuous(limits = c(0, 70),
                     breaks = seq(0, 70, 10)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.85,0.85),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(face = "bold"))



### FIGURE 7: landings ~ SOI + SAM -----------------------------------------------

png("Figures/FIG7_landings_by_SAMSOI.png", width = 12, height = 12, units = "in", res = 300)
ggpubr::ggarrange(f_landings_sam.plot, m_landings_sam.plot, 
                  f_landings_soi.plot, m_landings_soi.plot,
                  ncol = 2, nrow = 2,
                  #align = "hv"
                  widths = c(1, 0.92))
dev.off()





# APPENDIX  ====================================================================
#####
#####
#####
#####
# BAYESIAN MULTIPLE MULTIVARIATE REGRESSION ------------------------------------

## Uses the BNSP package
## Currently: using overall variables (i.e. not separating day)

# Data pre-processing -----------------------------------------------------

# ~ FEMALES

# Set id and year to factor
femalegps$id <- as.factor(femalegps$id)
femalegps$Year <- as.factor(femalegps$Year)

# Remove days that round to 0 and NA personalities
femalegps <- subset(femalegps, round(tripduration.days) > 0 & 
                       !is.na(boldness))

# Isolate key variables
female_mod <- femalegps[,c("id", "Year", "Age", "DeploymentID", "sex", "totalpathdistance.km",
                           "propARSvsTravel_hmm", "total_landings_hmm", "boldness",
                           "SAMIndex", "IODIndex", "SOIIndex", "breeding_success",
                           "attempted_breeding")]

## Introduce weighting variable to account for multiple breeding records per individual
female_mod <- female_mod %>% 
                  group_by(id, Year) %>%
                  mutate(n_recs = n()) %>%
                  mutate(pointweight = 1/n_recs) %>%
  data.frame()

# ~ MALES

# Set id and year to factors
malegps$id <- as.factor(malegps$id)
malegps$Year <- as.factor(malegps$Year)

# Remove days that round to 0 and NA personalities
malegps <- subset(malegps, round(tripduration.days) > 0 & 
                     !is.na(boldness))

male_mod <- malegps[,c("id", "Year", "Age", "DeploymentID", "sex", "totalpathdistance.km",
                           "propARSvsTravel_hmm", "total_landings_hmm", "boldness",
                           "SAMIndex", "IODIndex", "SOIIndex", "breeding_success",
                           "attempted_breeding")]

## Introduce weighting variable to account for multiple breeding records per individual
male_mod <- male_mod %>% 
              group_by(id, Year) %>%
              mutate(n_recs = n()) %>%
               mutate(pointweight = 1/n_recs) %>%
  data.frame()


## Visualise variables

# Histograms
par(mfrow = c(3,2))
hist(female_mod$total_landings_hmm, main = "Female landings"); 
  hist(male_mod$total_landings_hmm, main = "Male landings");
hist(female_mod$travSearch_ratio, main = "Female travel:search"); 
  hist(male_mod$travSearch_ratio, main = "Male travel:search"); 
hist(female_mod$logmaxdistance.km, main = "Female max distance"); 
  hist(male_mod$logmaxdistance.km, main = "Male max distance")

# Correlations
car::scatterplotMatrix(~total_landings_hmm + travSearch_ratio + logmaxdistance.km, 
                       data = female_mod)

car::scatterplotMatrix(~total_landings_hmm + travSearch_ratio + logmaxdistance.km, 
                       data = male_mod)


## Visualise breeding success

# Get one row per individual
all_data <- rbind(male_mod, female_mod)

all_data <- all_data %>% 
              group_by(id, Year) %>%
              filter(row_number()==1) %>%
            data.frame()

ggplot(all_data %>% 
         filter(!is.na(breeding_success)) %>% 
         group_by(Age) %>% 
         summarise(mean = mean(breeding_success)), 
       aes(x = Age, y = mean)) + 
  geom_point()

# Distribution of ages

ggplot(all_data %>% 
         filter(!is.na(breeding_success)), aes(x = Age)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ sex)



# Identify model structures using glmmTMB ---------------------------------

## Path distance
m_pathDist <- glmmTMB(totalpathdistance.km ~ SOIIndex * boldness + IODIndex * boldness +
                      SAMIndex * boldness + (1|id) + (1|Year), 
              family = gaussian("log"),
              data = male_mod)

plot(simulateResiduals(m_pathDist))

# Gamma looks terrible; gaussian better but not great; poisson inappropriate 
# Log gaussian is a little underdispersed but this will usually bias p-values 
# to the conservative side. This is probably because the model is a little
# over-parameterised 

## Prop Search:Travel
m_travSearch <- glmmTMB(propARSvsTravel_hmm ~ SOIIndex * boldness + IODIndex * boldness +
              SAMIndex * boldness + (1|id) + (1|Year), 
             family = beta_family,
             data = male_mod)

plot(simulateResiduals(m_travSearch))

# Gaussian looks surprisingly good but need to account for proportional response.
# Beta family looks lovely! Some quantile deviations but really minor


## Number of landings
m_landings <- glmmTMB(total_landings_hmm ~ SOIIndex * boldness + IODIndex * boldness +
               SAMIndex * boldness + (1|id) + (1|Year),
               family = poisson(link = "log"),
               data = male_mod)

plot(simulateResiduals(m_landings))

# No distribution seems perfect, though poisson probably makes most sense. Residual
# vs predicted plot probably supports a non-linear relationship but can't really
# see this in raw data and no a priori predictions to this effect


## Fitness 

# Just look at age first

# Scale age
male_mod$age_s <- as.numeric(scale(male_mod$Age, center = TRUE, scale = TRUE))
all_data$age_s <- as.numeric(scale(all_data$Age, center = TRUE, scale = TRUE))
rsdata$age_s <- as.numeric(scale(rsdata$Age, center = TRUE, scale = TRUE))

# One record per individual
male_mod2 <- male_mod %>%
              group_by(id, Year) %>%
              filter(row_number()==1) %>%
            data.frame()

male_mod2 <- subset(male_mod, !is.na(Age))
male_mod2 <- male_mod[order(male_mod$id),]

m_fitness <- glmmTMB(breeding_success ~ age_s + I(age_s^2) + (1|Year), 
        data = filter(rsdata, !is.na(breeding_success)), 
        family = "binomial")

plot(simulateResiduals(m_fitness))

summary(m_fitness)

m_m1_age <- ggpredict(m_fitness, terms = "age_s [all]")
p_m_m1_age <- plot(m_m1_age, show.title = FALSE)  +
  scale_y_continuous(limits = c(0,1)) + 
  theme_jt + 
  labs(x = "Age (years, scaled)", y = "P| Breeding")


m_fitnessForage <- glmmTMB(breeding_success ~ SOIIndex * boldness + 
                  IODIndex * boldness + SAMIndex * boldness + (1|id) + (1|Year), 
                  family = binomial,
                  data = male_mod)

plot(simulateResiduals(m_fitnessForage))

# QQplot mostly fine (a little underdispersed); some within-group deviartions but
# sounds like not necessarily an issue


## Fitness ~ Age
male_mod.sub <- subset(male_mod, !is.na(Age))

m_fitnessAge <- glmmTMB(breeding_success ~ Age + Age^2 + boldness + Age:boldness + (1|id) + (1|Year), 
                        family = binomial(link = "logit"),
                        data = male_mod.sub)

plot(simulateResiduals(m_fitnessAge))
summary(m_fitnessAge)

m_fitnessAge_plot <- ggpredict(m_fitnessAge, terms = "Age [all]")
p_m_m1_age <- plot(m_fitnessAge_plot, show.title = FALSE)

# This doesn't really work in glmmTMB framework, but residuals look fine and have
# already tried different distributions in brms and figured this out


# Construct the Bayesian models ---------------------------------------------------------


## Specify models
brms_path <- bf(totalpathdistance.km ~ SOIIndex * boldness + IODIndex * boldness +
              SAMIndex * boldness + 
              (1|id) + (1|Year), 
              family = gaussian(link = "log"))

brms_ratio <- bf(propARSvsTravel_hmm ~ SOIIndex * boldness + IODIndex * boldness +
              SAMIndex * boldness +
              (1|id) + (1|Year), 
              family = Beta())

brms_landings <- bf(total_landings_hmm ~ SOIIndex * boldness + IODIndex * boldness +
              SAMIndex * boldness +
              (1|id) + (1|Year), 
              family = poisson(link = "log"))

brms_rs <- bf(breeding_success|weights(pointweight) ~  poly(Age,2)*boldness +
                SOIIndex * boldness + IODIndex * boldness + SAMIndex * boldness +
                (1|id) + (1|Year), 
                family = bernoulli(link = "logit"))

### MALES 

## Remove unknown ages ( n = 19 )
male_mod.sub <- subset(male_mod, !is.na(Age))
nrow(male_mod) - nrow(male_mod.sub)

## Testing models in turn
mod_test <- brm(brms_rs,
    set_rescor(FALSE),
    data = dat.sub, 
    cores = 3, chains = 3, 
    warmup = 1000, iter = 3000,
    control = list(adapt_delta = 0.9), 
    seed = 12)

summary(mod_test)

pp_check(mod_test, resp = "totalpathdistancekm", ndraws = 50) # fine 
pp_check(mod_test, resp = "totallandingshmm", ndraws = 50)  # fine
pp_check(mod_test, resp = "propARSvsTravelhmm", ndraws = 50)  # fine
pp_check(mod_test, resp = "Age", ndraws = 50)

conditional_effects(male_model_brms, "boldness", resp = "travSearchratio")




male_model_brms <-  brm(brms_path + brms_ratio + brms_landings + brms_rs,
                          set_rescor(FALSE),
                          data = male_mod.sub, 
                          cores = 3, chains = 3, 
                          warmup = 1000, iter = 3000,
                          control = list(adapt_delta = 0.9), 
                          seed = 12)

save(male_model_brms, file = "Data_outputs/climate_brm_M.RData")
print(summary(male_model_brms), digits = 2)  

pp_check(male_model_brms, resp = "totalpathdistancekm", ndraws = 50) # fine 
pp_check(male_model_brms, resp = "totallandingshmm", ndraws = 50)  # fine
pp_check(male_model_brms, resp = "propARSvsTravelhmm", ndraws = 50)  # fine

bayes_R2(male_model_brms)

conditional_effects(male_model_brms, "boldness", resp = "travSearchratio")



### FEMALES 
female_model_brms <-  brm(brms1 + brms2 + brms3 +
                          set_rescor(FALSE),
                        data = female_mod, 
                        cores = 3, chains = 3, 
                        warmup = 1000, iter = 3000,
                        control = list(adapt_delta = 0.9), 
                        seed = 12)

pp_check(female_model_brms, resp = "totalpathdistancekm", ndraws = 50) # fine 
pp_check(female_model_brms, resp = "totallandingshmm", ndraws = 50)  # fine
pp_check(female_model_brms, resp = "travSearchratio", ndraws = 50)  # fine

save(female_model_brms, file = "Data_outputs/climate_brm_F.RData")
print(summary(female_model_brms), digits = 2) 



# ~ #### VISUALISE RESULTS #### ------------------------------------------------

## (A) Estimate effects of climate indices on foraging metrics -----------------

sexes <- c("female", "male")
sex_lab <- c("F", "M")
indices <- c("IOD", "SAM", "SOI")

for (i in 1:length(indices)) {
  
  for (s in 1:length(sexes)) {
    
    # Find predictive range
    myindex <- paste0(indices[i], "Index")
    pred_range <- as.vector(quantile( get(paste0(sexes[s],"gps"))[, get("myindex")], 
                                     probs = c(0.025, 0.975))) # -0.2220  0.4535 
    predict_over <- seq(pred_range[1], pred_range[2], length.out = 14)
    
    # Use mean boldness and other climate variables 
    boldness_mean <- mean(na.omit(get(paste0(sexes[s],"gps"))$boldness))
    IOD_mean <- mean(na.omit(get(paste0(sexes[s],"gps"))$IODIndex))
    SOI_mean <- mean(na.omit(get(paste0(sexes[s],"gps"))$SOIIndex))
    SAM_mean <- mean(na.omit(get(paste0(sexes[s],"gps"))$SAMIndex))
    
    # Make a dataset to predict from
    newdat <- data.frame(predIndex = predict_over,
                           boldness = boldness_mean,
                           SOIIndex = SOI_mean,
                           SAMIndex = SAM_mean,
                           IODIndex = IOD_mean)
    
    newdat[,myindex] <- NULL
    colnames(newdat)[1] <- myindex
    
    # Set model
    brms_mod <- get(paste0(sexes[s], "_model_brms"))
    
    ## Predicted landings ------------------------------------------------------
    landings_pred <- fitted(brms_mod, newdata = newdat,
                            resp = "totallandingshmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
    
    landings_pred <- as.data.frame(landings_pred)
    colnames(landings_pred) <- predict_over
    
    landings_df <- gather(landings_pred, index, landings, factor_key = TRUE)
    colnames(landings_df)[1] <- myindex
    
    ## Predicted distance ------------------------------------------------------
    distance_pred <- fitted(brms_mod, newdata = newdat,
                            resp = "totalpathdistancekm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
    
    distance_pred <- as.data.frame(distance_pred)
    colnames(distance_pred) <- predict_over
    
    distance_df <- gather(distance_pred, index, distance, factor_key = TRUE)
    colnames(distance_df)[1] <- myindex
    
    ## Predicted ratio travel : search -----------------------------------------
    ratio_pred <- fitted(brms_mod, newdata = newdat,
                         resp = "travSearchratio", ndraws = 1000, 
                         summary = FALSE, re_formula = NA)
    
    ratio_pred <- as.data.frame(ratio_pred)
    colnames(ratio_pred) <- predict_over
    
    ratio_df <- gather(ratio_pred, index, ratio, factor_key = TRUE)
    colnames(ratio_df)[1] <- myindex
    
    ## Bind up the data --------------------------------------------------------
    
    climate_predictions <- cbind(landings_df, distance_df[,2], ratio_df[,2])
    colnames(climate_predictions)[3:4] <- c("distance", "travSearch")
    climate_predictions[,1] <- as.numeric(as.character(climate_predictions[,1]))
    
    save(climate_predictions, file = paste0("Data_outputs/", indices[i], "_predictions_", 
                                            sex_lab[s], ".RData"))
    
    assign(paste0(sexes[s], "_climate_predictions"), climate_predictions)
    
    # Summarise the predictions
    climate_prediction_summ <- climate_predictions %>%
      group_by_at(1) %>%
      median_qi(.width = c(.66, .95)) %>% 
      mutate(sex = sexes[s]) %>%
      data.frame()
    
    climate_prediction_summ[,1] <- as.numeric(as.character(climate_prediction_summ[,1]))
    
    save(climate_prediction_summ, 
         file = paste0("Data_outputs/", indices[i], "_predictions_", 
                       sex_lab[s],"_summ.RData"))
    
    assign(paste0(sexes[s], "_climate_prediction_summ"), climate_prediction_summ)
    
  # }
  # 
  # ## (B) Make the figures ----------------------------------------------------------
  # 
  # # Prepare the data
  # twosex_prediction <- rbind(female_climate_predictions %>% mutate(sex = "Females"), 
  #                            male_climate_predictions %>% mutate(sex = "Males"))
  # 
  # twosex_prediction_summ <- rbind(female_climate_prediction_summ %>% mutate(sex = "Females"), 
  #                                male_climate_prediction_summ %>% mutate(sex = "Males"))
  # Find climate ranges
  climate_range <- rbind(malegps, femalegps) %>% 
    select(!!sym(myindex), sex) %>% 
    mutate(sex = if_else(sex == "Male", "Males", "Females"))
  
  # Set sex and relevant colour (individual plots)
  twosex_prediction_summ <- climate_prediction_summ
  twosex_prediction_summ$sex <- ifelse(s == 2, "Males", "Females")
    
  plot_col <- ifelse(s == 2, "orange", "dodgerblue")
  
  ## [1] Build the landings plot -------------------------------------------------
  landings_plot <- ggplot(data = twosex_prediction_summ, 
                          aes(x = twosex_prediction_summ[,myindex], 
                              y = landings)) + # , group = sex)
    stat_pointinterval(data = climate_range, aes(x = climate_range[,myindex], y = 4), 
                       point_size = 3.5, colour = "darkgrey") +
    geom_pointinterval(data = twosex_prediction_summ, 
                       aes(ymin = landings.lower, ymax = landings.upper), 
                       point_size = 2.5, colour = "black") +
    geom_point(data = twosex_prediction_summ, aes(colour = plot_col), size = 2) +
    theme_bw() + 
    theme(axis.title = element_text(colour = "black", size = 14), 
          axis.text = element_text(colour = "black", size = 12), 
          strip.text = element_text(colour = "black", size = 14), 
          legend.position = "none") +
    labs(x = paste0(indices[i], " Index"), y = "Total number of landings") + 
    #facet_wrap(~sex) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(4, 35)) + 
    #scale_colour_manual(values = c("dodgerblue", "orange"))
    scale_colour_manual(values = plot_col)
  
  # png(paste0("Figures/", indices[i], "_landings.png"), 
  #     width = 8, height = 6, units = "in", res = 300)
  # print(landings_plot)
  # dev.off()
  
  png(paste0("Figures/", sex_lab[s], "_", indices[i], "_landings.png"), 
           width = 5, height = 6, units = "in", res = 300)
  print(landings_plot)
  dev.off()
  
  ## [2] Build the ratios plot -------------------------------------------------
  ratio_plot <- ggplot(data = twosex_prediction_summ, 
                          aes(x = twosex_prediction_summ[,myindex], 
                              y = travSearch)) + # , group = sex
    stat_pointinterval(data = climate_range, aes(x = climate_range[,myindex], y = 0.625), 
                       point_size = 3.5, colour = "darkgrey") +
    geom_pointinterval(data = twosex_prediction_summ, 
                       aes(ymin = travSearch.lower, ymax = travSearch.upper), 
                       point_size = 2.5, colour = "black") +
    geom_point(data = twosex_prediction_summ, aes(colour = plot_col), size = 2) +  
    theme_bw() + 
    theme(axis.title = element_text(colour = "black", size = 14), 
          axis.text = element_text(colour = "black", size = 12), 
          strip.text = element_text(colour = "black", size = 14), 
          legend.position = "none") +
    labs(x = paste0(indices[i], " Index"), y = "Ratio travel : search time") + 
    #facet_wrap(~sex) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0.6, 1.4)) + 
    #scale_colour_manual(values = c("dodgerblue", "orange")) +
    scale_colour_manual(values = plot_col) 
  
  # png(paste0("Figures/", indices[i], "_ratioTravSearch.png"), 
  #     width = 8, height = 6, units = "in", res = 300)
  # print(ratio_plot)
  # dev.off()
  
  png(paste0("Figures/", sex_lab[s], "_", indices[i], "_ratioTravSearch.png"), 
      width = 5, height = 6, units = "in", res = 300)
  print(ratio_plot)
  dev.off()
  
  ## [3] Build the distance plot -----------------------------------------------
  distance_plot <- ggplot(data = twosex_prediction_summ, 
                       aes(x = twosex_prediction_summ[,myindex], 
                           y = distance)) + # , group = sex
    stat_pointinterval(data = climate_range, aes(x = climate_range[,myindex], y = 3000), 
                       point_size = 3.5, colour = "darkgrey") +
    geom_pointinterval(data = twosex_prediction_summ, 
                       aes(ymin = distance.lower, ymax = distance.upper), 
                       point_size = 2.5, colour = "black") +
    geom_point(data = twosex_prediction_summ, aes(colour = plot_col), size = 2) + 
    theme_bw() + 
    theme(axis.title = element_text(colour = "black", size = 14), 
          axis.text = element_text(colour = "black", size = 12), 
          strip.text = element_text(colour = "black", size = 14), 
          legend.position = "none") +
    labs(x = paste0(indices[i], " Index"), y = "Total path distance (km)") + 
    #facet_wrap(~sex) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(3000, 6000)) + 
    #scale_colour_manual(values = c("dodgerblue", "orange")) 
    scale_colour_manual(values = plot_col)
  
  # png(paste0("Figures/", indices[i], "_distance.png"), 
  #     width = 8, height = 6, units = "in", res = 300)
  # print(distance_plot)
  # dev.off()
  
  png(paste0("Figures/", sex_lab[s], "_", indices[i], "_distance.png"), 
      width = 5, height = 6, units = "in", res = 300)
  print(distance_plot)
  dev.off()
  }
}


# END --------------------------------------------------------------------------


## Isolate model variables
female_mod <- femalegps[,c("tripduration.days", "propSearch_hmm", "logtotalpathdistance.km",
         "SAMIndex", "SOIIndex", "IODIndex", "Year", "id", "boldness"),]
female_mod$id <- as.factor(female_mod$id)

# Remove days that round to 0 and NA personalities
female_mod <- subset(female_mod, round(tripduration.days) > 0 & 
                                 !is.na(boldness))

male_mod <- malegps[,c("tripduration.days", "propSearch_hmm", "logtotalpathdistance.km",
                           "SAMIndex", "SOIIndex", "IODIndex", "Year", "id", "boldness"),]
male_mod$id <- as.factor(male_mod$id)

# Remove days that round to 0 and NA personalities
male_mod <- subset(male_mod, round(tripduration.days) > 0 & 
                       !is.na(boldness))





####OLD ANALYSIS BELOW###### ------------------------------------------------------



### Check the data -------------------------------------------------------------

## Look at correlations in climate data

# Load the data
load("Data_outputs/climate_indices_1960-2021.RData")

# Scatterplots
plot(climate_data[,c(3:5)])
car::scatterplotMatrix(~ SAMIndex + SOIIndex + IODIndex, data = climate_data)

### Visualise the presumed relationships

# Direct: trip duration ~ SAM + IOD + SOI
# Mediator A: prop search ~ SAM + IOD + SOI
# Mediator B: path dist ~ prop search + SAM + IOD + SOI
# Indirect: trip duration ~ prop search + path dist

sex <- "male"

get(paste0(sex, "_mod")) %>% 
  select(tripduration.days, propSearch_hmm, logtotalpathdistance.km,
         SAMIndex, SOIIndex, IODIndex) %>% 
  psych::pairs.panels()

##### Run models separately for each climate index - set index here
climInd <- "SOIIndex" # "SAMIndex" "IODIndex"

## (A) Mediation analysis via causal models  -----------------------------------

## IV on mediator 1: prop. search
 # Males
m1_M <- glmmTMB(propSearch_hmm ~ get(climInd) * boldness + (1|Year:id), data = male_mod)
summary(m1_M)

plot(simulateResiduals(m1_M)) # good

# Females
m1_F <- glmmTMB(propSearch_hmm ~ get(climInd) * boldness + (1|Year:id), data = female_mod)
summary(m1_F)

plot(simulateResiduals(m1_F)) 
testDispersion(m1_F) # minor deviation but histogram looks good - carrying on

## Result: No effect of any indices on proportion search


## IV on mediator 2: log( path distance )
 # Males
m2_M <- glmmTMB(logtotalpathdistance.km ~ propSearch_hmm + get(climInd) * boldness + 
                  (1|Year:id), data = male_mod)
summary(m2_M)

plot(simulateResiduals(m2_M)) # fine
testDispersion(m2_M)

plot(ggpredict(m2_M, terms = c("propSearch_hmm"))) # More search = shorter trips
plot(ggpredict(m2_M, terms = c("SOIIndex"))) # Shorter trips with positive SOI
plot(ggpredict(m2_M, terms = c("propSearch_hmm", "SOIIndex"))) 

 # Females
m2_F <- glmmTMB(logtotalpathdistance.km ~ propSearch_hmm +get(climInd) * boldness + 
                  (1|Year:id), data = femalegps)
summary(m2_F)

plot(simulateResiduals(m2_F)) # okay - histogram good
testDispersion(m2_F)

plot(ggpredict(m2_F, terms = c("propSearch_hmm"))) # More search = shorter trips
plot(ggpredict(m2_F, terms = c("SOIIndex"))) # Shorter trips with positive SOI
plot(ggpredict(m2_F, terms = c("propSearch_hmm", "SOIIndex"))) 


## Mediator on DV : trip duration ~ search + distance
 # Males

## zero truncated neg binom: ' Zero-truncated negative binomial regression is 
## used to model count data for which the value zero cannot occur and for 
## which over dispersion exists'
m3_M <-
  glmmTMB(round(tripduration.days) ~ propSearch_hmm + logtotalpathdistance.km +
      (1 | Year:id),
    family = truncated_nbinom2(link = "log"),
    data = male_mod, 
    control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
    na.action = na.omit)

summary(m3_M) 

plot(simulateResiduals(m3_M))

plotResiduals(m3_M, malegps$logtotalpathdistance.km[!is.na(malegps$propSearch_hmm)])
plotResiduals(m3_M, na.omit(malegps$propSearch_hmm))

plot(ggpredict(m3_M, terms = c("propSearch_hmm", "logpathdistance.km"))) # more search + longer distances associated with longer trips

# Females
m3_F <-
  glmmTMB(round(tripduration.days) ~ propSearch_hmm + logtotalpathdistance.km +
            (1 | Year:id),
          family = truncated_nbinom2(link = "log"),
          data = female_mod, 
          control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")),
          na.action = na.omit)

summary(m3_F) 
plot(simulateResiduals(m3_F))

plot(ggpredict(m3_F, terms = c("propSearch_hmm", "logpathdistance.km"))) # more search + longer distances associated with longer trips



## (B) SEM-style mediation model in brms -------------------------------------------------------------------------

## Look at personality * climate interactions for each index 
## NB Have to manually change climate index

# Define the individual models
brms1_M <- bf(propSearch_hmm ~ IODIndex * boldness + (1|Year:id))
brms2_M <- bf(logtotalpathdistance.km ~ propSearch_hmm + IODIndex * boldness + (1|Year:id))
# brms3_M <- bf(round(tripduration.days) | trunc(lb = 1) ~ 
#                 propSearch_hmm + logtotalpathdistance.km + 
#                 SAMIndex * boldness + (1|Year:id),
#               family = negbinomial(link = "log"))

brms3_M <- bf(round(tripduration.days) | trunc(lb = 1) ~ 
                propSearch_hmm + logtotalpathdistance.km + 
                IODIndex * boldness + (1|Year:id),
              family = "poisson")

male_model_brms <-  brm(brms1_M + brms2_M + brms3_M + set_rescor(FALSE),
                        data = male_mod, 
                        cores = 3, chains = 3, 
                        warmup = 1000, iter = 3000,
                        control = list(adapt_delta = 0.9), 
                        seed = 12)

save(male_model_brms, file = paste0("Data_outputs/", climInd, "_brm_M.RData"))
plot(male_model_brms)  
print(summary(male_model_brms), digits = 2)  

pp_check(male_model_brms, resp = "PropSearch", ndraws = 50)
pp_check(male_model_brms, resp = "logtotalpathdistancekm", ndraws = 50)
pp_check(male_model_brms, resp = "rounddur", ndraws = 50)

## Females

# Define the individual models
brms1_F <- bf(propSearch_hmm ~ IODIndex * boldness + (1|Year:id))
brms2_F <- bf(logtotalpathdistance.km ~ propSearch_hmm + IODIndex * boldness + (1|Year:id))
brms3_F <- bf(round(tripduration.days) | trunc(lb = 1) ~ 
                propSearch_hmm + logtotalpathdistance.km + 
                IODIndex * boldness + (1|Year:id),
              family = "poisson")

female_model_brms <-  brm(brms1_F + brms2_F + brms3_F + set_rescor(FALSE),
                        data = female_mod, 
                        cores = 3, chains = 3, 
                        warmup = 1000, iter = 3000,
                        control = list(adapt_delta = 0.9), 
                        seed = 12)

save(female_model_brms, file = paste0("Data_outputs/", climInd, "_brm_F.RData"))
plot(female_model_brms)  
print(summary(female_model_brms), digits = 2)  

#pp_check(female_model_brms, resp = "propSearch_hmm", ndraws = 50)  
#pp_check(female_model_brms, resp = "logpathdistancekm", ndraws = 50)
#pp_check(female_model_brms, resp = "tripduration.days", ndraws = 50)


# ~ #### ~RESULTS~ #### --------------------------------------------------

allgps <- rbind(femalegps, malegps)
allgps$sex <- ifelse(allgps$sex == "Female", "Females", "Males")

# (A) Estimate the effect of IOD on path distance ------------------------------

# Estimate the effect of IOD on path distance by carrying forward the posterior prediction through the mediation model

## Males -----------------------------------------------------------------------

# Load the model
load("Data_outputs/IODIndex_brm_M.RData")

mean(malegps$IODIndex)
quantile(malegps$IODIndex, probs = c(0.025, 0.975)) 

# Predict from -0.22 to 0.50 in increments of 0.05
predict_over <- seq(-0.22, 0.49, by = 0.05)

# Use mean boldness, as not important in models
boldness_mean <- mean(na.omit(malegps$boldness))

male_list <- list()

for (i in 1:length(predict_over)) {
  
  print(i)
  
  newdat_m <- data.frame(IODIndex = predict_over[i],
                         boldness = boldness_mean)
  
  propsearch_pred <- fitted(male_model_brms, newdata = newdat_m,
                            resp = "propSearchhmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_m2 <- expand.grid(IODIndex = predict_over[i], 
                           boldness = boldness_mean,
                           propSearch_hmm = as.vector(propsearch_pred))
  
  pathl_pred <- fitted(male_model_brms, newdata = newdat_m2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(propSearch_hmm = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  pathl_pred$boldness <- boldness_mean
  pathl_pred$IODIndex <- predict_over[i]
  pathl_pred$pathdist <- exp(pathl_pred$logtotalpathdistance.km) 
 
  male_list[[i]] <- data.frame(IODIndex = predict_over[i], 
                               pathdist = pathl_pred$pathdist)
  
}

male_prediction <- do.call(rbind, male_list)

male_prediction_summ <- male_prediction %>%
  group_by(IODIndex) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Males")


save(male_prediction, file = "Data_outputs/IOD_pathDist_M.RData")
save(male_prediction_summ, file = "Data_outputs/IOD_pathDist_M_summ.RData")

## Females ---------------------------------------------------------------------

# Load the model
load("Data_outputs/IODIndex_brm_F.RData")

mean(femalegps$IODIndex)
quantile(femalegps$IODIndex, probs = c(0.025, 0.975)) 

# Predict from -0.22 to 0.45 in increments of 0.05
predict_over <- seq(-0.22, 0.45, by = 0.05)

# Use mean boldness, as not important in models
boldness_mean <- mean(na.omit(femalegps$boldness))

female_list <- list()

for (i in 1:length(predict_over)) {
  
  print(i)
  
  newdat_f <- data.frame(IODIndex = predict_over[i],
                         boldness = boldness_mean)
  
  propsearch_pred <- fitted(female_model_brms, newdata = newdat_f,
                            resp = "propSearchhmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_f2 <- expand.grid(IODIndex = predict_over[i], 
                           boldness = boldness_mean,
                           propSearch_hmm = as.vector(propsearch_pred))
  
  pathl_pred <- fitted(female_model_brms, newdata = newdat_f2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(propSearch_hmm = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  pathl_pred$boldness <- boldness_mean
  pathl_pred$IODIndex <- predict_over[i]
  pathl_pred$pathdist <- exp(pathl_pred$logtotalpathdistance.km) 
  
  female_list[[i]] <- data.frame(IODIndex = predict_over[i], 
                               pathdist = pathl_pred$pathdist)
  
}

female_prediction <- do.call(rbind, female_list)

female_prediction_summ <- female_prediction %>%
  group_by(IODIndex) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Females")


save(female_prediction, file = "Data_outputs/IOD_pathDist_F.RData")
save(female_prediction_summ, file = "Data_outputs/IOD_pathDist_F_summ.RData")


## ~ FIGURE ~ Path distance ~ IOD + Sex ----------------------------------------

twosex_prediction <- rbind(male_prediction %>% mutate(sex = "Males"), 
                           female_prediction %>% mutate(sex = "Females"))
twosex_prediction_summ <- rbind(male_prediction_summ, female_prediction_summ)

## Find IOD ranges
iodranges <- rbind(malegps, femalegps) %>% 
  select(IODIndex, sex) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

## Build the plot
iod_plot <- ggplot(data = twosex_prediction, 
                   aes(x = IODIndex, y = pathdist, group = sex)) + 
  geom_point(data = allgps, aes(x = IODIndex, y = totalpathdistance.km, 
                                col = sex), alpha = 0.5,
             position = position_jitter()) +
  stat_pointinterval(data = iodranges, aes(x = IODIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  geom_pointinterval(data = twosex_prediction_summ, aes(ymin = .lower, ymax = .upper), 
                     point_size = 2.5, colour = "black") +
  geom_point(data = twosex_prediction_summ, aes(colour = sex), size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = "none") +
  labs(x = "IOD Index", y = "Total path distance, km") + 
  facet_wrap(~sex) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  scale_colour_manual(values = c("dodgerblue", "orange"))

png("Figures/iod_effects_raw_points.png", width = 8, height = 6, units = "in", res = 300)
iod_plot
dev.off()




# (B) Estimate the effect of SOI on path distance ------------------------------

# Estimate the effect of SOI on path distance by carrying forward the posterior prediction through the mediation model

## Males -----------------------------------------------------------------------

# Load the model
load("Data_outputs/SOIIndex_brm_M.RData")

quantile(malegps$SOIIndex, probs = c(0.025, 0.975)) 

# Predict from -3.6 to 4.5 in increments of 0.5
predict_over <- seq(-3.6, 4.5, by = 0.5)

# Use mean boldness, as not important in models
boldness_mean <- mean(na.omit(malegps$boldness))

male_list <- list()

for (i in 1:length(predict_over)) {
  
  print(i)
  
  newdat_m <- data.frame(SOIIndex = predict_over[i],
                         boldness = boldness_mean)
  
  propsearch_pred <- fitted(male_model_brms, newdata = newdat_m,
                            resp = "propSearchhmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_m2 <- expand.grid(SOIIndex = predict_over[i], 
                           boldness = boldness_mean,
                           propSearch_hmm = as.vector(propsearch_pred))
  
  pathl_pred <- fitted(male_model_brms, newdata = newdat_m2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(propSearch_hmm = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  pathl_pred$pathdist <- exp(pathl_pred$logtotalpathdistance.km) 
  
  male_list[[i]] <- data.frame(SOIIndex = predict_over[i], 
                               pathdist = pathl_pred$pathdist)
  
}

male_prediction <- do.call(rbind, male_list)

male_prediction_summ <- male_prediction %>%
  group_by(SOIIndex) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Males")


save(male_prediction, file = "Data_outputs/SOI_pathDist_M.RData")
save(male_prediction_summ, file = "Data_outputs/SOI_pathDist_M_summ.RData")

## Females ---------------------------------------------------------------------

# Load the model
load("Data_outputs/SOIIndex_brm_F.RData")

quantile(femalegps$SOIIndex, probs = c(0.025, 0.975)) 

# Predict from -3.6 to 4.5 in increments of 0.5
predict_over <- seq(-3.6, 4.5, by = 0.5)

# Use mean boldness, as not important in models
boldness_mean <- mean(na.omit(femalegps$boldness))

female_list <- list()

for (i in 1:length(predict_over)) {
  
  print(i)
  
  newdat_f <- data.frame(SOIIndex = predict_over[i],
                         boldness = boldness_mean)
  
  propsearch_pred <- fitted(female_model_brms, newdata = newdat_f,
                            resp = "propSearchhmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_f2 <- expand.grid(SOIIndex = predict_over[i], 
                           boldness = boldness_mean,
                           propSearch_hmm = as.vector(propsearch_pred))
  
  pathl_pred <- fitted(female_model_brms, newdata = newdat_f2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(propSearch_hmm = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  pathl_pred$pathdist <- exp(pathl_pred$logtotalpathdistance.km) 
  
  female_list[[i]] <- data.frame(SOIIndex = predict_over[i], 
                                 pathdist = pathl_pred$pathdist)
  
}

female_prediction <- do.call(rbind, female_list)

female_prediction_summ <- female_prediction %>%
  group_by(SOIIndex) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Females")


save(female_prediction, file = "Data_outputs/SOI_pathDist_F.RData")
save(female_prediction_summ, file = "Data_outputs/SOI_pathDist_F_summ.RData")


## ~ FIGURE ~ Path distance ~ SOI + Sex ----------------------------------------

twosex_prediction <- rbind(male_prediction %>% mutate(sex = "Males"), 
                           female_prediction %>% mutate(sex = "Females"))
twosex_prediction_summ <- rbind(male_prediction_summ, female_prediction_summ)

## Find SOI ranges
soiranges <- rbind(malegps, femalegps) %>% 
  select(SOIIndex, sex) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

## Build the plot
soi_plot <- ggplot(data = twosex_prediction, 
                   aes(x = SOIIndex, y = pathdist, group = sex)) +
  geom_point(data = allgps, aes(x = SOIIndex, y = totalpathdistance.km, 
                                col = sex), alpha = 0.5,
             position = position_jitter()) +
  stat_pointinterval(data = soiranges, aes(x = SOIIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  geom_pointinterval(data = twosex_prediction_summ, aes(ymin = .lower, ymax = .upper), 
                     point_size = 2.5, colour = "black") +
  geom_point(data = twosex_prediction_summ, aes(colour = sex), size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = "none") +
  labs(x = "SOI Index", y = "Total path distance, km") + 
  facet_wrap(~sex) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  scale_colour_manual(values = c("dodgerblue", "orange"))

png("Figures/soi_effects_raw_points.png", width = 8, height = 6, units = "in", res = 300)
soi_plot
dev.off()




# (C) Estimate the effect of SAM on path distance ------------------------------

# Estimate the effect of SAM on path distance by carrying forward the posterior prediction through the mediation model

## Males -----------------------------------------------------------------------

# Load the model
load("Data_outputs/SAMIndex_brm_M.RData")

quantile(malegps$SAMIndex, probs = c(0.025, 0.975)) 

# Predict from -1.87 to 4.36 in increments of 0.5
predict_over <- seq(-1.87, 4.36, by = 0.5)

# Use mean boldness, as not important in models
boldness_mean <- mean(na.omit(malegps$boldness))

male_list <- list()

for (i in 1:length(predict_over)) {
  
  print(i)
  
  newdat_m <- data.frame(SAMIndex = predict_over[i],
                         boldness = boldness_mean)
  
  propsearch_pred <- fitted(male_model_brms, newdata = newdat_m,
                            resp = "propSearchhmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_m2 <- expand.grid(SAMIndex = predict_over[i], 
                           boldness = boldness_mean,
                           propSearch_hmm = as.vector(propsearch_pred))
  
  pathl_pred <- fitted(male_model_brms, newdata = newdat_m2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(propSearch_hmm = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  pathl_pred$boldness <- boldness_mean
  pathl_pred$SAMIndex <- predict_over[i]
  pathl_pred$pathdist <- exp(pathl_pred$logtotalpathdistance.km) 
  
  male_list[[i]] <- data.frame(SAMIndex = predict_over[i], 
                               pathdist = pathl_pred$pathdist)
  
}

male_prediction <- do.call(rbind, male_list)

male_prediction_summ <- male_prediction %>%
  group_by(SAMIndex) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Males")


save(male_prediction, file = "Data_outputs/SAM_pathDist_M.RData")
save(male_prediction_summ, file = "Data_outputs/SAM_pathDist_M_summ.RData")

## Females ---------------------------------------------------------------------

# Load the model
load("Data_outputs/SAMIndex_brm_F.RData")

quantile(femalegps$SAMIndex, probs = c(0.025, 0.975)) 

# Predict from -1.87 to 4.36 in increments of 0.5
predict_over <- seq(-1.87, 4.36, by = 0.5)

# Use mean boldness, as not important in models
boldness_mean <- mean(na.omit(femalegps$boldness))

female_list <- list()

for (i in 1:length(predict_over)) {
  
  print(i)
  
  newdat_f <- data.frame(SAMIndex = predict_over[i],
                         boldness = boldness_mean)
  
  propsearch_pred <- fitted(female_model_brms, newdata = newdat_f,
                            resp = "propSearchhmm", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_f2 <- expand.grid(SAMIndex = predict_over[i], 
                           boldness = boldness_mean,
                           propSearch_hmm = as.vector(propsearch_pred))
  
  pathl_pred <- fitted(female_model_brms, newdata = newdat_f2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(propSearch_hmm = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  pathl_pred$boldness <- boldness_mean
  pathl_pred$SAMIndex <- predict_over[i]
  pathl_pred$pathdist <- exp(pathl_pred$logtotalpathdistance.km) 
  
  female_list[[i]] <- data.frame(SAMIndex = predict_over[i], 
                                 pathdist = pathl_pred$pathdist)
  
}

female_prediction <- do.call(rbind, female_list)

female_prediction_summ <- female_prediction %>%
  group_by(SAMIndex) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Females")


save(female_prediction, file = "Data_outputs/SAM_pathDist_F.RData")
save(female_prediction_summ, file = "Data_outputs/SAM_pathDist_F_summ.RData")


## ~ FIGURE ~ Path distance ~ SAM + Sex ----------------------------------------

twosex_prediction <- rbind(male_prediction %>% mutate(sex = "Males"), 
                           female_prediction %>% mutate(sex = "Females"))
twosex_prediction_summ <- rbind(male_prediction_summ, female_prediction_summ)

## Find SAM ranges
SAMranges <- rbind(malegps, femalegps) %>% 
  select(SAMIndex, sex) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

## Build the plot
SAM_plot <- ggplot(data = twosex_prediction, 
                   aes(x = SAMIndex, y = pathdist, group = sex)) + 
  geom_point(data = allgps, aes(x = SAMIndex, y = totalpathdistance.km, 
                                col = sex), alpha = 0.5,
             position = position_jitter()) +
  stat_pointinterval(data = SAMranges, aes(x = SAMIndex, y = 4), 
                     point_size = 3.5, colour = "darkgrey") +
  geom_pointinterval(data = twosex_prediction_summ, aes(ymin = .lower, ymax = .upper), 
                     point_size = 2.5, colour = "black") +
  geom_point(data = twosex_prediction_summ, aes(colour = sex), size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = "none") +
  labs(x = "SAM Index", y = "Total path distance, km") + 
  facet_wrap(~sex) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  scale_colour_manual(values = c("dodgerblue", "orange"))

png("Figures/SAM_effects_raw_points.png", width = 8, height = 6, units = "in", res = 300)
SAM_plot
dev.off()








#-----------------

# UNEDITED SCRIPT BELOW ---------------------------------------------------

# Setting up casual models for the effect of wind and personality on foraging performance # 

lapply(c("DHARMa", "tidyverse", "glmmTMB", "brms", "ggeffects", "tidybayes", "patchwork"), FUN = "library", character.only = TRUE)

setwd("C:/Users/jackt/OneDrive/Documents/Liverpool/WAAL GPS ANALYSES/ForagingTripVariation")

femalegps <- read.csv("femaleincubationtrips.csv", header = TRUE)
malegps <- read.csv("maleincubationtrips.csv", header = TRUE) # stick to Jan and Feb to be sure...

setwd("C:/Users/jackt/OneDrive/Desktop/Albatross paper pipeline")

# Get the aggregated data for each trip 
femalegps$year <- as.factor(femalegps$year) ; femalegps$id <- as.factor(femalegps$id)
femalegps$logmaxdistance.km <- log(femalegps$maxdistance.km)
femalegps$logtotalpathdistance.km <- log(femalegps$totalpathdistance.km)
femalegps$meanwindspeed.kph <- femalegps$meanwindspeed*3.6

malegps$year <- as.factor(malegps$year) ; malegps$id <- as.factor(malegps$id)
malegps$logmaxdistance.km <- log(malegps$maxdistance.km)
malegps$meanwindspeed.kph <- malegps$meanwindspeed*3.6
malegps$logtotalpathdistance.km <- log(malegps$totalpathdistance.km)
malegps <- subset(malegps, meanwindspeed.kph > 20) # outliers 

# merge them together and give a new Trackid that links to the interpolated tracks
gpsdata <- rbind(malegps, femalegps) %>% 
  mutate(trackid = NA, 
         start_time = as.POSIXct(start_time, format = "%d/%m/%Y %H:%M")) %>% 
  select(-trackid)

# Get the interpolated 'fine-scale' data. Interpolated for every 15 minutes. 
interpol <- read.csv("C:/Users/jackt/OneDrive/Documents/Liverpool/WAAL GPS ANALYSES/WAAL_GPS_interpolated15mins.csv", header = TRUE)
interpol$time <- as.POSIXct(interpol$time, format = "%d/%m/%Y %H:%M")

start_times <- interpol %>% 
  group_by(DeploymentID, trackid) %>% 
  summarise(start_time = min(time))

gpsdata <- left_join(gpsdata, start_times, by = c("DeploymentID", "start_time")) %>% 
  filter(!is.na(trackid))
gpsdata$year <- droplevels(gpsdata$year)

# now incorporate the other information from hmms
hmm <- read.csv("C:/Users/jackt/OneDrive/Documents/Liverpool/WAAL GPS ANALYSES/hmm_datout.csv", header = TRUE)
hmm_summ <- hmm %>% 
  rename(trackid = ID) %>% 
  group_by(trackid, State) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(trackid, names_from = State, values_from = n) %>% 
  mutate(PropRest = Rest/(Rest + Search + Travel), 
         PropSearch = Search/(Rest + Search + Travel), 
         PropTravel = Travel/(Rest + Search + Travel), 
         RatioTravelSearch = Travel/(Search + Travel))

gpsdata <- left_join(gpsdata, hmm_summ)

malegps <- filter(gpsdata, sex == "Male")
malegps$meanwindspeed_s <- as.numeric(scale(malegps$meanwindspeed.kph)) # to help the non-brms models converge

femalegps <- filter(gpsdata, sex == "Female")
femalegps$meanwindspeed_s <- as.numeric(scale(femalegps$meanwindspeed.kph))  # to help the non-brms models converge

#--------- Fit the models -------------------

# Prop Search ~ wind speed
m1 <- glmmTMB(PropSearch ~ meanwindspeed_s + I(meanwindspeed_s^2) + (1|Year) + (1|id), data = malegps)
summary(m1) # both linear and quadratic effect of wind speed
plot(ggpredict(m1, terms = c("meanwindspeed_s [all]")))

# path length

head(malegps)
# Max distance ~ PropSearch + wind speed
m2 <- glmmTMB(logtotalpathdistance.km ~ PropSearch + meanwindspeed_s + I(meanwindspeed_s^2) + (1|year) + (1|id), data = malegps)
summary(m2)
plot(ggpredict(m2, terms = c("PropSearch")))
plot(ggpredict(m2, terms = c("meanwindspeed_s [all]")))

# Trip duration ~ Prop Search + max distance
round(tripduration.days) ~ logmaxdistance.km + PropSearch
m3 <- glmmTMB(round(tripduration.days) ~ PropSearch + logtotalpathdistance.km +  (1|year) + (1|id), family = "poisson", data = malegps)
summary(m3) 
#plot(simulateResiduals(m3))
plot(ggpredict(m3, terms = c("PropSearch", "logtotalpathdistance.km")))

#----------------------  THE FEMALE MODEL ACCORDING TO CORNIOLEY ---------

# Prop Search ~ wind speed
f1 <- glmmTMB(PropSearch ~ meanwindspeed_s + I(meanwindspeed_s^2) + (1|year) + (1|id), data = femalegps)
summary(f1) # both linear and quadratic effect of wind speed
plot(ggpredict(f1, terms = c("meanwindspeed_s [all]")))

# Max distance ~ PropSearch + wind speed
f2 <- glmmTMB(logtotalpathdistance.km ~ PropSearch  +  meanwindspeed_s + I(meanwindspeed_s^2) + (1|year) + (1|id), data = femalegps)
summary(f2)
plot(ggpredict(m2, terms = c("PropSearch")))
plot(ggpredict(m2, terms = c("meanwindspeed_s [all]")))

# Trip duration ~ Prop Search + max distance
f3 <- glmmTMB(round(tripduration.days) ~ PropSearch + logtotalpathdistance.km +  (1|year) + (1|id), family = "poisson", data = femalegps)
summary(f3) 
#plot(simulateResiduals(m3))
plot(ggpredict(f3, terms = c("logtotalpathdistance.km", "PropSearch")))

# ------------ SET UP THE MODELS AS MULTIVARIATE [MEDIATION] MODELS IN BRMS  ---------

## NB Add interaction with personality to each of climate variables

# Male mediation model

model_1 <- bf(PropSearch ~  meanwindspeed.kph + I(meanwindspeed.kph^2) + (1|Year) + (1|id))
model_2 <- bf(logtotalpathdistance.km ~ PropSearch + meanwindspeed.kph +  I(meanwindspeed.kph^2) + (1|Year) + (1|id))
model_3 <- bf(round(tripduration.days) ~ PropSearch + logtotalpathdistance.km + (1|Year) + (1|id), 
              family = poisson)

male_model_brms <-  brm(model_1 + model_2 + model_3 + set_rescor(FALSE),
                        data = malegps, 
                        cores = 3, chains = 3, 
                        warmup = 1000, iter = 3000,
                        control = list(adapt_delta = 0.9), 
                        seed = 12)
plot(male_model_brms)  
print(summary(male_model_brms), digits = 2)  
#pp_check(model_brms, resp = "meanspeedkph", ndraws = 50)  
#pp_check(model_brms, resp = "logmaxdistancekm", ndraws = 50)
#pp_check(model_brms, resp = "tripdurationdays", ndraws = 50)

# Female mediation model
female_model_brms <-  brm(model_1 + model_2 + model_3 + set_rescor(FALSE),
                          data = femalegps, 
                          cores = 3, chains = 3, 
                          warmup = 1000, iter = 3000,
                          control = list(adapt_delta = 0.9), 
                          seed = 12)
#plot(female_model_brms)
print(summary(female_model_brms), digits = 2)


# ---------- Estimate the effect of SOI on path distance -----------------------

# Estimate the effect of wind speed on trip duration by carrying forward the posterior prediction through the mediation model
# Wind speed: Males
mean(malegps$meanwindspeed.kph)
quantile(malegps$meanwindspeed.kph, probs = c(0.025, 0.975)) # predict from 25 to 43

male_list <- list()

for (i in 25:43) {
  
  print(i)
  
  newdat_m <- data.frame(meanwindspeed.kph = i)
  propsearch_pred <- fitted(male_model_brms, newdata = newdat_m,
                            resp = "PropSearch", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_m2 <- expand.grid(meanwindspeed.kph = i, PropSearch = as.vector(propsearch_pred))
  pathl_pred <- fitted(male_model_brms, newdata = newdat_m2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(PropSearch = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  tripd_pred <- fitted(male_model_brms, newdata = pathl_pred,
                       resp = "roundtripdurationdays", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  tripd_pred <- as.matrix(diag(tripd_pred)) %>% 
    data.frame() %>% 
    bind_cols(pathl_pred) %>% 
    rename(tripduration.days = names(.)[1])
  
  male_list[[i - 24]] <- data.frame(meanwindspeed = i, tripduration = tripd_pred$tripduration.days)
  
}

male_prediction <- do.call(rbind, male_list)
male_prediction_summ <- male_prediction %>%
  group_by(meanwindspeed) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Males")




#-------------------------------------------------------------------------------------

# Estimate the effect of wind speed on trip duration by carrying forward the posterior prediction through the mediation model
# Wind speed: Males
mean(malegps$meanwindspeed.kph)
quantile(malegps$meanwindspeed.kph, probs = c(0.025, 0.975)) # predict from 25 to 43

male_list <- list()

for (i in 25:43) {
  
  print(i)
  
  newdat_m <- data.frame(meanwindspeed.kph = i)
  propsearch_pred <- fitted(male_model_brms, newdata = newdat_m,
                            resp = "PropSearch", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_m2 <- expand.grid(meanwindspeed.kph = i, PropSearch = as.vector(propsearch_pred))
  pathl_pred <- fitted(male_model_brms, newdata = newdat_m2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(PropSearch = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  tripd_pred <- fitted(male_model_brms, newdata = pathl_pred,
                       resp = "roundtripdurationdays", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  tripd_pred <- as.matrix(diag(tripd_pred)) %>% 
    data.frame() %>% 
    bind_cols(pathl_pred) %>% 
    rename(tripduration.days = names(.)[1])
  
  male_list[[i - 24]] <- data.frame(meanwindspeed = i, tripduration = tripd_pred$tripduration.days)
  
}

male_prediction <- do.call(rbind, male_list)
male_prediction_summ <- male_prediction %>%
  group_by(meanwindspeed) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Males")

#------------------------------

# Wind speed: Females
mean(femalegps$meanwindspeed.kph)
quantile(femalegps$meanwindspeed.kph, probs = c(0.025, 0.975)) # predict from 22 to 40

female_list <- list()

for (i in 22:40) {
  
  print(i)
  newdat_f <- data.frame(meanwindspeed.kph = i)
  propsearch_pred <- fitted(female_model_brms, newdata = newdat_f,
                            resp = "PropSearch", ndraws = 1000, 
                            summary = FALSE, re_formula = NA)
  
  newdat_f2 <- expand.grid(meanwindspeed.kph = i, PropSearch = as.vector(propsearch_pred))
  pathl_pred <- fitted(female_model_brms, newdata = newdat_f2,
                       resp = "logtotalpathdistancekm", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  pathl_pred <- as.matrix(diag(pathl_pred)) %>% 
    data.frame() %>% 
    mutate(PropSearch = propsearch_pred) %>% 
    rename(logtotalpathdistance.km = names(.)[1])
  
  tripd_pred <- fitted(female_model_brms, newdata = pathl_pred,
                       resp = "roundtripdurationdays", ndraws = 1000, 
                       summary = FALSE, re_formula = NA)
  
  tripd_pred <- as.matrix(diag(tripd_pred)) %>% 
    data.frame() %>% 
    bind_cols(pathl_pred) %>% 
    rename(tripduration.days = names(.)[1])
  
  
  female_list[[i - 21]] <- data.frame(meanwindspeed = i, tripduration = tripd_pred$tripduration.days)
  
}

female_prediction <- do.call(rbind, female_list)
female_prediction_summ <- female_prediction %>%
  group_by(meanwindspeed) %>%
  median_qi(.width = c(.66, .95)) %>% 
  mutate(sex = "Females")

twosex_prediction <- rbind(male_prediction %>% mutate(sex = "Males"), female_prediction %>% mutate(sex = "Females"))
twosex_prediction_summ <- rbind(male_prediction_summ, female_prediction_summ)

windranges <- rbind(malegps, femalegps) %>% 
  select(meanwindspeed.kph, sex) %>% 
  mutate(sex = if_else(sex == "Male", "Males", "Females"))

wind_plot <- ggplot(data = twosex_prediction, aes(x = meanwindspeed, y = tripduration, group = sex)) + 
  #geom_jitter(aes(colour = sex), alpha = 0.1, width = 0.2) + 
  stat_pointinterval(data = windranges, aes(x = meanwindspeed.kph, y = 4), point_size = 3.5, colour = "darkgrey") +
  geom_pointinterval(data = twosex_prediction_summ, aes(ymin = .lower, ymax = .upper), point_size = 2.5, colour = "black") +
  geom_point(data = twosex_prediction_summ, aes(colour = sex), size = 2) +
  theme_bw() + 
  theme(axis.title = element_text(colour = "black", size = 14), 
        axis.text = element_text(colour = "black", size = 12), 
        strip.text = element_text(colour = "black", size = 14), 
        legend.position = "none") +
  labs(x = "Mean wind speed, kph", y = "Trip duration, days") + 
  facet_wrap(~sex) + 
  scale_y_continuous(breaks = seq(4,14,2), labels = seq(4,14,2)) + 
  scale_colour_manual(values = c("dodgerblue", "orange"))


#----------------- GENERATE OTHER RELEVANT PLOTS FROM THE DATA --------------------------------

# I will make conditional effects plots for all of the effects (without personality)
theme_cond <-  theme(axis.title = element_text(colour = "black", size = 13), 
                     axis.text = element_text(colour = "black", size = 12.5), 
                     strip.text = element_text(size = 13))

# first sub-model (Prop search ~ wind speed)
p1 <- conditional_effects(female_model_brms, "meanwindspeed.kph", resp = "PropSearch")[[1]] %>% 
  mutate(sex = "Females") 
p1 <- ggplot(p1, aes(x = meanwindspeed.kph, y = estimate__)) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "dodgerblue", alpha = 0.3) +
  geom_line(colour = "dodgerblue", lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  labs(x = "Mean wind speed, kph", y = "Proportion of time\n in search behaviour") + 
  facet_wrap(~sex)  +
  ylim(c(0.2, 0.7))


p2 <- conditional_effects(male_model_brms, "meanwindspeed.kph", resp = "PropSearch")[[1]] %>% 
  mutate(sex = "Males") 
p2 <- ggplot(p2, aes(x = meanwindspeed.kph, y = estimate__)) +  
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "orange", alpha = 0.3) +
  geom_line(colour = "orange", lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  labs(x = "Mean wind speed, kph", y = "Proportion of time\n in search behaviour") + 
  facet_wrap(~sex) +
  scale_x_continuous(breaks = seq(20, 50, 5)) +
  ylim(c(0.2, 0.7))

p1 | p2

# Second sub-model (path length ~ mean speed)
p3 <- conditional_effects(female_model_brms, "PropSearch", resp = "logtotalpathdistancekm")[[1]] %>% 
  mutate(sex = "Females") 
p3 <- ggplot(p3, aes(x = PropSearch, y = estimate__)) +  
  geom_ribbon(aes(ymin =  lower__, ymax =  upper__), fill = "dodgerblue", alpha = 0.3) +
  geom_line(colour = "dodgerblue", lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  labs(x = "Proportion of time\n in search behaviour", y = "ln(Path length), km") + 
  facet_wrap(~sex) + 
  scale_x_continuous(breaks = seq(0.1, 0.8, 0.1)) +
  ylim(c(3.5, 10))


p4 <- conditional_effects(male_model_brms, "PropSearch", resp = "logtotalpathdistancekm")[[1]] %>% 
  mutate(sex = "Males") 
p4 <- ggplot(p4, aes(x = PropSearch, y = estimate__)) +  
  geom_ribbon(aes(ymin =  lower__, ymax =  upper__), fill = "orange", alpha = 0.3) +
  geom_line(colour = "orange", lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  labs(x =  "Proportion of time\n in search behaviour", y = "ln(Path length), km") + 
  facet_wrap(~sex) + 
  scale_x_continuous(breaks = seq(0.1, 0.8, 0.1)) +
  ylim(c(3.5, 10))

p3 | p4

# Second sub-model (trip distance ~ wind speed)
p5 <- conditional_effects(female_model_brms, "meanwindspeed.kph", resp = "logtotalpathdistancekm")[[1]] %>% 
  mutate(sex = "Females") 
p5 <- ggplot(p5, aes(x = meanwindspeed.kph, y = estimate__)) +  
  geom_ribbon(aes(ymin =  lower__, ymax =  upper__), fill = "dodgerblue", alpha = 0.3) +
  geom_line(colour = "dodgerblue", lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  labs(x = "Mean wind speed, kph", y = "ln(Path length), km") + 
  facet_wrap(~sex) + 
  ylim(c(4.5, 10))

p6 <- conditional_effects(male_model_brms, "meanwindspeed.kph", resp = "logtotalpathdistancekm")[[1]] %>% 
  mutate(sex = "Males") 
p6 <- ggplot(p6, aes(x = meanwindspeed.kph, y = estimate__)) +  
  geom_ribbon(aes(ymin =  lower__, ymax =  upper__), fill = "orange", alpha = 0.3) +
  geom_line(colour = "orange", lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  labs(x = "Mean wind speed, kph", y = "ln(Path length), km") + 
  facet_wrap(~sex) + 
  ylim(c(4.5, 10))

p5 | p6

# Third sub-model (trip duration ~ prop search + foraging distance)

p7 <- conditional_effects(female_model_brms, effects = c("logtotalpathdistance.km"), conditions = data.frame(PropSearch = c(0.2, 0.4, 0.6)), resp = "roundtripdurationdays")[[1]] %>% 
  mutate(sex = "Females")
p7 <- ggplot(p7, aes(x = logtotalpathdistance.km, y = estimate__, group = PropSearch , colour = as.factor(PropSearch), 
                     fill = as.factor(PropSearch))) +  
  geom_ribbon(aes(ymin =  lower__, ymax =  upper__), alpha = 0.3, colour = NA) +
  geom_line(lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  guides(fill = "none") +
  theme(legend.position = c(0.4, 0.7)) +
  labs(x = "ln(Path length), km", y = "Trip durations, days", colour = "Prop. search") + 
  facet_wrap(~sex) + 
  scale_colour_manual(values = c("aquamarine", "dodgerblue", "blue")) +
  scale_fill_manual(values = c("aquamarine", "dodgerblue", "blue"))

p8 <- conditional_effects(male_model_brms, effects = c("logtotalpathdistance.km"), conditions = data.frame(PropSearch = c(0.2, 0.4, 0.6)), resp = "roundtripdurationdays")[[1]] %>% 
  mutate(sex = "Males")
p8 <- ggplot(p8, aes(x = logtotalpathdistance.km, y = estimate__, group = PropSearch , colour = as.factor(PropSearch), 
                     fill = as.factor(PropSearch))) +    
  geom_ribbon(aes(ymin =  lower__, ymax =  upper__), alpha = 0.3, colour = NA) +
  geom_line(lwd = 1.5) +
  theme_bw() + 
  theme_cond + 
  guides(fill = "none") +
  theme(legend.position = c(0.4, 0.7)) +
  labs(x = "ln(Path length), km", y = "Trip durations, days", colour = "Prop. search") + 
  facet_wrap(~sex) + 
  scale_colour_manual(values = c("gold", "orange", "red")) +
  scale_fill_manual(values = c("gold", "orange", "red")) 

p7 | p8

plota <- p1 | p2
plotb <- p3 | p4
plotc <- p5 | p6
plotd <- p7 | p8

p1+ p2 + p3 + p4 + p5 + p6 + p7 + p8 + 
  plot_layout(ncol = 4, nrow = 2)


##---------------- GET STANDARDISED VARIABLES FOR THE MODEL COMPONENTS ------------------------

# FOR THIS NEED TO REFIT THE MODELS WITH WIND AS AN ORTHOGONAL POLYNOMIAL 
# THIS IS PREFERABLE OVER FITTING ORIGINALLY AS ORTHOGONAL POLYNOMIAL AND BACK-TRANSFORMING COEFFICIENTS, WHICH IS A PAIN FOR PREDICTION


# For the models without personality info

malegps <- data.frame(malegps, poly(malegps$meanwindspeed.kph, 2))
names(malegps)[51:52] <- c("windpoly", "wind2poly")

femalegps <- data.frame(femalegps, poly(femalegps$meanwindspeed.kph, 2))
names(femalegps)[51:52] <- c("windpoly", "wind2poly")


model_7 <- bf(PropSearch ~ windpoly + wind2poly + (1|year) + (1|id))
model_8 <- bf(logtotalpathdistance.km ~ PropSearch + windpoly + wind2poly + (1|year) + (1|id))

male_model_poly_brms <-  brm(model_7 + model_8 + model_3 + set_rescor(FALSE),
                             data = malegps, 
                             cores = 3, chains = 3, 
                             warmup = 1000, iter = 3000,
                             control = list(adapt_delta = 0.9), 
                             seed = 12)
print(summary(male_model_poly_brms), digits = 3)
bayes_R2(male_model_poly_brms) # get the bayes R2

female_model_poly_brms <-  brm(model_7 + model_8 + model_3 + set_rescor(FALSE),
                               data = femalegps, 
                               cores = 3, chains = 3, 
                               warmup = 1000, iter = 3000,
                               control = list(adapt_delta = 0.9), 
                               seed = 12)
print(summary(female_model_poly_brms), digits = 3)
bayes_R2(female_model_poly_brms) # get the bayes R2

# Get the standardised variables for males

get_variables(male_model_poly_brms)[1:11]
b_male_windspeed_search <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_PropSearch_windpoly")))*sd(malegps$windpoly)/sd(malegps$PropSearch) # ground speed ~ wind speed
mean(b_male_windspeed_search)  ; quantile(b_male_windspeed_search, probs = c(0.025, 0.975)) # very similar to ETI used by brms

b_male_windspeed2_search <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_PropSearch_wind2poly")))*sd(malegps$wind2poly)/sd(malegps$PropSearch) # ground speed ~ wind speed2
mean(b_male_windspeed2_search)  ; quantile(b_male_windspeed2_search, probs = c(0.025, 0.975))

b_male_propsearch_pathl <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_logtotalpathdistancekm_PropSearch")))*sd(malegps$PropSearch)/sd(malegps$logtotalpathdistance.km) # foraging distance ~ ground speed
mean(b_male_propsearch_pathl)  ; quantile(b_male_propsearch_pathl, probs = c(0.025, 0.975))

b_male_windspeed_pathl <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_logtotalpathdistancekm_windpoly")))*sd(malegps$windpoly)/sd(malegps$logtotalpathdistance.km) # foraging distance ~ wind speed
mean(b_male_windspeed_pathl)  ; quantile(b_male_windspeed_pathl, probs = c(0.025, 0.975))

b_male_windspeed2_pathl <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_logtotalpathdistancekm_wind2poly")))*sd(malegps$wind2poly)/sd(malegps$logtotalpathdistance.km) # foraging distance ~ wind speed
mean(b_male_windspeed2_pathl)  ; quantile(b_male_windspeed2_pathl, probs = c(0.025, 0.975))

# For the poisson model standardisation see https://jslefche.github.io/sem_book/coefficients.html
R2.bayes.male <- cor(round(malegps$tripduration.days), fitted(male_model_poly_brms, type = "response", resp = "roundtripdurationdays")[,1])^2
sd.yhat.bayes.male <- sqrt(var(fitted(male_model_poly_brms, scale = "linear", resp = "roundtripdurationdays")[,1])/R2.bayes.male)

b_male_propsearch_tripd <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_roundtripdurationdays_PropSearch"))) * sd(malegps$PropSearch)/sd.yhat.bayes.male # trip duration ~ ground speed
mean(b_male_propsearch_tripd)  ; quantile(b_male_propsearch_tripd, probs = c(0.025, 0.975))

b_male_pathl_tripd <- as.vector(unlist(as_draws_list(male_model_poly_brms, variable = "b_roundtripdurationdays_logtotalpathdistance.km"))) * sd(malegps$logtotalpathdistance.km)/sd.yhat.bayes.male # trip duration ~ ground speed # trip duration ~ ground speed
mean(b_male_pathl_tripd)  ; quantile(b_male_pathl_tripd, probs = c(0.025, 0.975))

# Get the standardised variables for females
get_variables(female_model_poly_brms)[1:11]
b_female_windspeed_search <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_PropSearch_windpoly")))*sd(femalegps$windpoly)/sd(femalegps$PropSearch) # ground speed ~ wind speed
mean(b_female_windspeed_search)  ; quantile(b_female_windspeed_search, probs = c(0.025, 0.975)) # very similar to ETI used by brms

b_female_windspeed2_search <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_PropSearch_wind2poly")))*sd(femalegps$wind2poly)/sd(femalegps$PropSearch) # ground speed ~ wind speed2
mean(b_female_windspeed2_search)  ; quantile(b_female_windspeed2_search, probs = c(0.025, 0.975))

b_female_propsearch_pathl <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_logtotalpathdistancekm_PropSearch")))*sd(femalegps$PropSearch)/sd(femalegps$logtotalpathdistance.km) # foraging distance ~ ground speed
mean(b_female_propsearch_pathl)  ; quantile(b_female_propsearch_pathl, probs = c(0.025, 0.975))

b_female_windspeed_pathl <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_logtotalpathdistancekm_windpoly")))*sd(femalegps$windpoly)/sd(femalegps$logtotalpathdistance.km) # foraging distance ~ wind speed
mean(b_female_windspeed_pathl)  ; quantile(b_female_windspeed_pathl, probs = c(0.025, 0.975))

b_female_windspeed2_pathl <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_logtotalpathdistancekm_wind2poly")))*sd(femalegps$wind2poly)/sd(femalegps$logtotalpathdistance.km) # foraging distance ~ wind speed
mean(b_female_windspeed2_pathl)  ; quantile(b_female_windspeed2_pathl, probs = c(0.025, 0.975))

# For the poisson model standardisation see https://jslefche.github.io/sem_book/coefficients.html
R2.bayes.female <- cor(round(femalegps$tripduration.days), fitted(female_model_poly_brms, type = "response", resp = "roundtripdurationdays")[,1])^2
sd.yhat.bayes.female <- sqrt(var(fitted(female_model_poly_brms, scale = "linear", resp = "roundtripdurationdays")[,1])/R2.bayes.female)

b_female_propsearch_tripd <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_roundtripdurationdays_PropSearch"))) * sd(femalegps$PropSearch)/sd.yhat.bayes.female # trip duration ~ ground speed
mean(b_female_propsearch_tripd)  ; quantile(b_female_propsearch_tripd, probs = c(0.025, 0.975))

b_female_pathl_tripd <- as.vector(unlist(as_draws_list(female_model_poly_brms, variable = "b_roundtripdurationdays_logtotalpathdistance.km"))) * sd(femalegps$logtotalpathdistance.km)/sd.yhat.bayes.female # trip duration ~ ground speed # trip duration ~ ground speed
mean(b_female_pathl_tripd)  ; quantile(b_female_pathl_tripd, probs = c(0.025, 0.975))



################## END ##################################################################
# save.image("Foraging Trip Variation_PropSearch.RData")

