# -------------------------------------
# Script: RS_by_climate.R
# Author: Dr Natasha Gillies
# Purpose: Lookg at variation in reproductive success with climate indices
# Notes: Based on work by Dr Jack Thorley
# Date: 2023-09-26
# -------------------------------------

# Load functions & packages ----------------------------------------------------

# Functions

# Packages
packages <- c("ggplot2", "dplyr", "glmmTMB", "magrittr", "ggeffects", "tidyr",
              "lme4", "lmerTest", "sjPlot")

## Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())
# 
# if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
# }

## Load packages
invisible(lapply(packages, library, character.only = TRUE))

## For reference: boldness quantiles

# Sex     upperQuan lowerQuan midQuan
# <fct>       <dbl>     <dbl>   <dbl>
#   1 Females      1.80     -1.37 -0.0937
# 2 Males        1.40     -1.95 -0.376 

## Set colours
female_col <- "#FFC20A"
female_fill <- "#f7dc8b"
male_col <- "#0C7BDC"
male_fill <- "#91c1eb"


# Function
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




# PREPARE THE DATASETS =========================================================

# Load the fitness data
waal_rs <- read.csv("Data_original/WAAL_RS_glmmTMBdata.csv", header = TRUE)

## Remove individuals reproducing before 7 (probably erroneous) or those with no age info
waal_rs %<>% 
  filter(AFR >=6) %>%
  filter(StatutBaguage == "P", !is.na(Sex))

## Set attempted breeding and breeding success to binomial variables
waal_rs %<>% 
  mutate(attempted_breeding = if_else(StateCode != 3, 1, 0), 
         breeding_success = case_when(StateCode == 0 ~ 0, 
                                              StateCode %in% c(1,2) ~ 1, 
                                              TRUE ~ as.numeric(NA))) 

## Add a variable for the previous year's reproductive attempt
waal_rs %<>%  
  group_by(Metalring) %>% 
  mutate(prevyear = case_when(lag(breeding_success) == "0" ~ "failedrep", 
                              lag(breeding_success) == "1" ~ "successfulrep", 
                              TRUE ~ "no attempt")) %>% 
  slice(-1) %>% # Remove their first ever record as this is one year before their first attempt
  data.frame()


# Change formatting of data and extract relevant variables
waal_rs$year <- as.factor(waal_rs$year)
waal_rs$prevyear <- factor(waal_rs$prevyear, levels = c("no attempt", "failedrep", "successfulrep"))


## Process the climate data ----------------------------------------------------

load("Data_original/IOD_monthly.RData")
load("Data_original/SOI_monthly.RData")
load("Data_original/SAM_monthly.RData")

# Average climate indices in pre-breeding (Sep-Nov) and breeding (Dec - Feb) periods

## Pre-breeding ## 

climate_prebreeding <- soi %>%
  group_by(Year) %>% 
  summarise(avgSOI_prebreeding = mean(SOIIndex[month %in% 9:11])) %>%
  mutate(Year = Year + 1) %>% 
  left_join(sam %>% 
              group_by(Year) %>% 
              summarise(avgSAM_prebreeding  = mean(SAMIndex[month %in% 9:11])) %>% 
              mutate(Year = Year + 1)) %>% 
  left_join(iod %>% 
              group_by(Year) %>% 
              summarise(avgIOD_prebreeding = mean(IODIndex[month %in% 9:11])) %>%
              mutate(Year = Year + 1)) %>% 
 mutate(Year = as.factor(Year)) %>%
 rename(year = Year)

## Breeding ##

climate_breeding <- soi %>% 
  #mutate(Year = ifelse(month == "12", Year + 1, Year)) %>%    
  group_by(Year) %>% 
  summarise(avgSOI_breeding = mean(SOIIndex[month %in% 1:4])) %>% 
  left_join(sam %>% 
              #mutate(year = ifelse(month == "12", year + 1, year)) %>% 
              group_by(Year) %>% 
              summarise(avgSAM_breeding = mean(SAMIndex[month %in% 1:4]))) %>% 
  left_join(iod %>% 
              #mutate(year = ifelse(month == "12", year + 1, year)) %>% 
              group_by(Year) %>% 
              summarise(avgIOD_breeding = mean(IODIndex[month %in% 1:4]))) %>% 
  mutate(Year = as.factor(Year)) %>%
  rename(year = Year)


## Merge into full dataset
waal_rs <- climate_breeding %>% 
  left_join(climate_prebreeding) %>% 
  right_join(waal_rs) %>% 
  data.frame()


# VISUALISE THE DATA ===========================================================

### Breeding success ~ Age
ggplot(waal_rs %>% 
         filter(!is.na(breeding_success)) %>% 
         group_by(Age) %>% 
         summarise(mean = mean(breeding_success)), 
       aes(x = as.factor(Age), y = mean)) + 
  geom_point() +
  geom_smooth() # single threshold should be good and best fitted to the data


### Distribution of ages
ggplot(waal_rs %>% 
         filter(!is.na(breeding_success)), aes(x = Age)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ Sex)


### Distribution of AFR
ggplot(waal_rs %>% 
         group_by(Metalring) %>% 
         slice(1), aes(x = AFR)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ Sex)

### Breeding success ~ climate
ggplot(waal_rs, 
       aes(x = avgIOD_breeding, y = breeding_success)) + 
  geom_point() +
  geom_smooth() 

# FIT THE MODELS ===============================================================

# Individual level effects of climate on reproductive success ------------------

### Isolate variables and rename
waal_rs %<>% dplyr::select(c("Metalring", "Sex", "year", "ReproCode", "Age", "AFR", 
                             "boldness_BLUP_mean", "attempted_breeding", "breeding_success",
                             "prevyear",  "avgSOI_breeding", "avgSAM_breeding",    
                             "avgIOD_breeding", "avgSOI_prebreeding", "avgSAM_prebreeding",
                             "avgIOD_prebreeding")) %>%
  rename(ring = Metalring) 


# Separate males and females
female_rs <- waal_rs %>% filter(Sex == "F" & !is.na(boldness_BLUP_mean))
male_rs <- waal_rs %>% filter(Sex == "M" & !is.na(boldness_BLUP_mean))

### Standardise the variables
female_rs %<>% 
  mutate(age_s = scale(Age),
         AFR_s = scale(AFR),
         boldness_s = scale(boldness_BLUP_mean),
         avgSAM_breeding_s = scale(avgSAM_breeding),
         avgSAM_prebreeding_s = scale(avgSAM_prebreeding),
         avgSOI_breeding_s = scale(avgSOI_breeding),
         avgSOI_prebreeding_s = scale(avgSOI_prebreeding),
         avgIOD_breeding_s = scale(avgIOD_breeding),
         avgIOD_prebreeding_s = scale(avgIOD_prebreeding)) 

male_rs %<>% 
  mutate(age_s = scale(Age),
         AFR_s = scale(AFR),
         boldness_s = scale(boldness_BLUP_mean),
         avgSAM_breeding_s = scale(avgSAM_breeding),
         avgSAM_prebreeding_s = scale(avgSAM_prebreeding),
         avgSOI_breeding_s = scale(avgSOI_breeding),
         avgSOI_prebreeding_s = scale(avgSOI_prebreeding),
         avgIOD_breeding_s = scale(avgIOD_breeding),
         avgIOD_prebreeding_s = scale(avgIOD_prebreeding)) 


## How does the probability of breeding success change with climate ------------

### SAM ========================================================================

#### FEMALES -------------------------------------------------------------------
f_SAM_glmm <- glmmTMB(breeding_success ~ age_s*avgSAM_prebreeding_s*boldness_s +   
                   I(age_s^2)*avgSAM_prebreeding_s*boldness_s +
                   age_s*avgSAM_breeding_s*boldness_s + 
                   I(age_s^2)*avgSAM_breeding_s*boldness_s +
                   prevyear +
                   (1|year/ring), 
                 data =  female_rs, 
                 family = "binomial")

summary(f_SAM_glmm)

## Remove non-significant random effects
f_SAM_glmm <- update(f_SAM_glmm, ~ age_s + I(age_s^2) + avgSAM_prebreeding_s +
                           avgSAM_breeding_s + boldness_s + prevyear +
                           I(age_s^2):avgSAM_breeding_s +
                           (1|year/ring))

summary(f_SAM_glmm)

f_SAM_glmm <- update(f_SAM_glmm, ~ age_s + I(age_s^2) + avgSAM_prebreeding_s +
                       avgSAM_breeding_s + boldness_s + prevyear + (1|year/ring))

tab_model(f_SAM_glmm, show.stat = T)

f_SAM_pred <- data.frame(ggpredict(f_SAM_glmm, terms = c("age_s [all]")))
f_SAM_pred$age <- (f_SAM_pred$x*sd(female_rs$Age)) + mean(female_rs$Age)

# Build the plot
p_f_age <- ggplot(f_SAM_pred, aes(x = age, y = predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = female_fill) +
  geom_line(size = 1.1, col = female_col) +
  labs(x = "Age, years", y = "P|Successful", title = "Females") + 
  ylim(c(0,1)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.15),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())




#### MALES ---------------------------------------------------------------------
m_SAM_glmm <- glmmTMB(breeding_success ~ age_s*avgSAM_prebreeding_s*boldness_s +   
                   I(age_s^2)*avgSAM_prebreeding_s*boldness_s +
                   age_s*avgSAM_breeding_s*boldness_s + 
                   I(age_s^2)*avgSAM_breeding_s*boldness_s +
                   prevyear + 
                   (1|year/ring), 
                 data =  male_rs, 
                 family = "binomial")

summary(m_SAM_glmm)

## Drop non-significant interactions
m_SAM_glmm <- update(m_SAM_glmm, ~ age_s + I(age_s^2) + avgSAM_prebreeding_s +
                       avgSAM_breeding_s + boldness_s + prevyear + (1|year/ring))

summary(m_SAM_glmm)
tab_model(m_SAM_glmm, show.stat = T)


# Get the plot data
m_SAM_pred <- data.frame(ggpredict(m_SAM_glmm, terms = c("age_s [all]")))
m_SAM_pred$age <- (m_SAM_pred$x*sd(male_rs$Age)) + mean(male_rs$Age)

## Build the plot
p_m_age <- ggplot(m_SAM_pred, aes(x = age, y = predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = male_fill) +
  geom_line(size = 1.1, col = male_col) +
  labs(x = "Age, years", y = "P|Successful", title = "Males") + 
  ylim(c(0,1)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.15),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())



### SOI ========================================================================

#### FEMALES -------------------------------------------------------------------
f_SOI_glmm <- glmmTMB(breeding_success ~ age_s*avgSOI_prebreeding_s*boldness_s +   
                   I(age_s^2)*avgSOI_prebreeding_s*boldness_s +
                   age_s*avgSOI_breeding_s*boldness_s + 
                   I(age_s^2)*avgSOI_breeding_s*boldness_s +
                   prevyear +
                   (1|year/ring), 
                 data =  female_rs, 
                 family = "binomial")

summary(f_SOI_glmm)

## Drop non-significant interactions
f_SOI_glmm <- update(f_SOI_glmm, ~ age_s +  I(age_s^2) + avgSOI_prebreeding_s + 
                          avgSOI_breeding_s + boldness_s + prevyear + (1|year/ring))

summary(f_SOI_glmm)
tab_model(f_SOI_glmm, show.stat = T)


#### MALES -------------------------------------------------------------------
m_SOI_glmm <- glmmTMB(breeding_success ~ age_s*avgSOI_prebreeding_s*boldness_s +   
                   I(age_s^2)*avgSOI_prebreeding_s*boldness_s +
                   age_s*avgSOI_breeding_s*boldness_s + 
                   I(age_s^2)*avgSOI_breeding_s*boldness_s +
                   prevyear +
                   (1|year/ring), 
                 data =  male_rs, 
                 family = "binomial")

summary(m_SOI_glmm)

## Drop non-significant interactions
m_SOI_glmm <- update(m_SOI_glmm, ~ age_s +  I(age_s^2) + avgSOI_prebreeding_s + 
                       avgSOI_breeding_s + boldness_s + prevyear + (1|year/ring))


summary(m_SOI_glmm)
tab_model(m_SOI_glmm, show.stat = T)


# FIGURE SX: can do an age effect plot here if useful ==========================

png("Figures/FIG4_individual_RS_by_climate.png", width = 12, height = 12, units = "in", res = 300)
ggpubr::ggarrange(p_f_SAM, p_m_SAM,
                  p_f_SOI, p_m_SOI,
                  ncol = 2, nrow = 2,
                  widths = c(1, 0.92))
dev.off()




### Population-level breeding success with climate -----------------------------

# Calculate mean RS per year
annual_rs <- waal_rs %>%
  group_by(year, avgSAM_breeding, avgSOI_breeding) %>%
  summarise(n_birds = n_distinct(ring),
            mean_rs = mean(breeding_success, na.rm = T)) %>%
  filter(n_birds > 3)

# Get weighted overall mean
annual_means <- waal_rs %>%
  group_by(year) %>%
  summarize(annual_rs = mean(breeding_success, na.rm = T),
            n_inds = n())

annual_means$weights <- annual_means$n_inds/sum(annual_means$n_inds)
xm <- weighted.mean(annual_means$annual_rs, annual_means$weights)
#[1] 0.7809899

weighted_var <- sum(annual_means$weights * (annual_means$annual_rs - xm)^2)
sqrt(weighted_var)

# Breeding success as a function of climate variables
RS_climate_glmer <- glmmTMB(breeding_success ~ avgSAM_breeding + avgSAM_prebreeding +
                              avgSOI_breeding + avgSOI_prebreeding + (1|year/ring), 
                            data = waal_rs, family = "binomial")

tab_model(RS_climate_glmer, show.stat = T)

## Plot results

climate_pred <- data.frame(ggpredict(RS_climate_glmer, terms = c("avgSAM_breeding [all]"))) %>% 
  rename(SAM = x) 

ggpredict(RS_climate_glmer, terms = c("avgSAM_breeding [0,1]"))
# avgSAM_breeding | Predicted |       95% CI
# ------------------------------------------
#   0 |      0.77 | [0.76, 0.79]
#   1 |      0.79 | [0.77, 0.80]

# SAM : 
RS_sam_plot <- ggplot() + 
  geom_point(data = annual_rs, aes(x = avgSAM_breeding, y = mean_rs, 
                                   size = n_birds)) +
  geom_ribbon(data = climate_pred,
              aes(x = SAM, ymin = conf.low, ymax = conf.high), alpha = 0.5, fill = "coral") +
  geom_line(data = climate_pred, aes(x = SAM, y = predicted), size = 1.1, col = "red") +
  ggrepel::geom_label_repel(data = annual_rs, aes(x = avgSAM_breeding, 
                                                  y = mean_rs, label = year),
                            segment.color = "darkorange",
                            min.segment.length = 0.1) +
  labs(x = "Arithmetic mean Southern Annular Mode (January to April)",
       y = "Population-level breeding success") +
  scale_y_continuous(limit = c(0.65, 1), breaks = seq(0.7, 1, by = 0.1)) +
  theme_bw() + 
  theme(legend.position = "none")


# SOI : 
RS_soi_plot <- ggplot() +
  geom_point(data = annual_rs, aes(x = avgSOI_breeding, y = mean_rs, 
                                   size = n_birds)) +
  ggrepel::geom_label_repel(data = annual_rs, aes(x = avgSOI_breeding, 
                                                  y = mean_rs, label = year),
                            segment.color = "darkorange",
                            min.segment.length = 0.1) +
  labs(x = "Arithmetic mean Southern Oscillation Index (January to April)",
       y = "Population-level breeding success") +
  scale_y_continuous(limit = c(0.65, 1), breaks = seq(0.7, 1, by = 0.1)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# FIGURE 6: breeding success ~ climate =========================================

png("Figures/FIG6_RS_by_climate.png", width = 12, height = 6, units = "in", res = 300)
ggpubr::ggarrange(RS_sam_plot, RS_soi_plot,
                  ncol = 2, nrow = 1,
                  widths = c(1, 0.92))
dev.off()






# APPENDIX ----------------------------------------------------------------




## ALTERNATIVE - SAM focused
f_SAM_pred <- data.frame(ggpredict(f_SAM_glmm, terms = c("age_s [all]", "avgSAM_breeding_s [all]"))) %>% 
  rename(SAM = group) 
f_SAM_pred$age <- (f_SAM_pred$x*sd(female_rs$Age)) + mean(female_rs$Age)
f_SAM_pred$SAM <- as.numeric(as.character(f_SAM_pred$SAM))
f_SAM_pred$SAM <- (f_SAM_pred$SAM*sd(female_rs$avgSAM_breeding)) + mean(female_rs$avgSAM_breeding)
f_age_levels <- unique(f_SAM_pred$age)[as.numeric(round(quantile(1:length(unique(f_SAM_pred$age)), 
                                                                 probs = c(0.10, 0.5, 0.90))))]

f_SAM_pred %<>% 
  filter(age %in% f_age_levels) %>% 
  mutate(age = as.numeric(as.character(age))) %>% 
  mutate(`Age` = case_when(age == min(age) ~ "Young (10 years)", 
                           age == max(age) ~ "Old (44 years)", 
                           TRUE ~ "Mid (27 years)"))

f_SAM_pred$Age <- factor(f_SAM_pred$Age, levels = c("Young (10 years)", "Mid (27 years)", "Old (44 years)"))


#### Make some predictions
maxSAM <- max(f_SAM_pred$SAM)
minSAM <- min(f_SAM_pred$SAM)

predict_diffs.prop(f_SAM_pred, "Age", "Young (10 years)", "SAM", maxSAM, minSAM)
predict_diffs.prop(f_SAM_pred, "Age", "Mid (27 years)", "SAM", maxSAM, minSAM)
predict_diffs.prop(f_SAM_pred, "Age", "Old (44 years)", "SAM", maxSAM, minSAM)

#### Build the plot
p_f_SAM <- ggplot(f_SAM_pred, aes(x = SAM, y = predicted, colour = Age, 
                                  group = Age, fill = Age)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, colour = NA) +
  geom_line(size = 1.1) +
  labs(x = "Southern Annular Mode (Breeding)", y = "P|Successful", title = "(a) Females") + 
  scale_colour_viridis_d() +
  ylim(c(0,1)) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.15),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())


### What is the probability of breeding success given have tried to breed? -----

## Females : 
f_RS_glmm <- glmmTMB(breeding_success ~ age_s*boldness_s + I(age_s^2)*boldness_s + 
                       AFR_s*boldness_s + I(AFR_s^2)*boldness_BLUP_mean + prevyear*boldness_s +                
                       (1|ring) + (1|year), 
                     data =  female_rs,
                     family = "binomial")
summary(f_RS_glmm)

f_RS_boldness <- ggpredict(f_RS_glmm, terms = c("boldness_s"))
p_f_RS_boldness <- plot(f_RS_boldness, show.title = FALSE) + 
  scale_y_continuous(limits = c(0,1))  + 
  theme_bw + 
  labs(x = "Boldness BLUP", y = "P| Successful")

f_RS_age <- ggpredict(f_RS_glmm, terms = c("age_s [all]", "boldness_s"))
p_f_RS_age <- plot(f_RS_age, show.title = FALSE, show.legend = FALSE) + 
  scale_y_continuous(limits = c(0,1))  + 
  theme_bw() +
  labs(x = "Age (years, scaled)", y = "P| Successful")

f_RS_prevyear <- ggpredict(f_RS_glmm, terms = c("prevyear [all]", "boldness_s"))
p_f_RS_prevyear <- plot(f_RS_prevyear, show.title = FALSE) + 
  scale_y_continuous(limits = c(0,1))  + 
  scale_x_continuous(breaks = c(1,2,3), labels = c("DNB", "FB", "SB")) +
  theme_bw() +
  labs(x = "Reproduction (t - 1)", y = "P| Successful")

f_RS_AFR <- ggpredict(f_RS_glmm, terms = c("AFR_s", "boldness_s"))
p_f_RS_AFR <- plot(f_RS_AFR, show.title = FALSE, show.legend = FALSE) + 
  scale_y_continuous(limits = c(0,1))  + 
  theme_jt + 
  labs(x = "AFR (years, scaled", y = "P| Successful")



