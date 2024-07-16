# -------------------------------------
# Script: RS_by_climate.R
# Author: Dr Natasha Gillies, Dr Jack Thorley
# Purpose: Looking at variation in reproductive success with climate indices
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



# PREPARE THE DATASETS =========================================================

# Load the fitness data
waal_rs <- read.csv("Data_inputs/WAAL_breedingSuccess_1965-2020.csv", header = TRUE)

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
  group_by(id) %>% 
  mutate(prevyear = case_when(lag(breeding_success) == "0" ~ "failedrep", 
                              lag(breeding_success) == "1" ~ "successfulrep", 
                              TRUE ~ "no attempt")) %>% 
  slice(-1) %>% # Remove their first ever record as this is one year before their first attempt
  data.frame()


# Change formatting of data and extract relevant variables
waal_rs$year <- as.factor(waal_rs$year)
waal_rs$prevyear <- factor(waal_rs$prevyear, levels = c("no attempt", "failedrep", "successfulrep"))

## Remove problematic years - 1982, 1983, 1984
## NB. Protocol was not established at this point and these values seem to be unusually high
waal_rs %<>% filter(!year %in% c(1982, 1983, 1984))


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
  geom_smooth() 


### Distribution of ages
ggplot(waal_rs %>% 
         filter(!is.na(breeding_success)), aes(x = Age)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ Sex)


### Distribution of AFR
ggplot(waal_rs %>% 
         group_by(id) %>% 
         slice(1), aes(x = AFR)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ Sex)


# FIT THE MODELS ===============================================================

# Individual level effects of climate on reproductive success ------------------

### Isolate variables and rename
waal_rs %<>% dplyr::select(c("id", "Sex", "year", "Age", "AFR",
                             "boldness_BLUP_mean", "attempted_breeding", "breeding_success",
                             "prevyear",  "avgSOI_breeding", "avgSAM_breeding",
                             "avgIOD_breeding", "avgSOI_prebreeding", "avgSAM_prebreeding",
                             "avgIOD_prebreeding")) %>%
  rename(ring = id)


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
         avgIOD_prebreeding_s = scale(avgIOD_prebreeding),
         year = as.numeric(as.character(year))) %>%
  filter(!is.na(breeding_success))

male_rs %<>%
  mutate(age_s = scale(Age),
         AFR_s = scale(AFR),
         boldness_s = scale(boldness_BLUP_mean),
         avgSAM_breeding_s = scale(avgSAM_breeding),
         avgSAM_prebreeding_s = scale(avgSAM_prebreeding),
         avgSOI_breeding_s = scale(avgSOI_breeding),
         avgSOI_prebreeding_s = scale(avgSOI_prebreeding),
         avgIOD_breeding_s = scale(avgIOD_breeding),
         avgIOD_prebreeding_s = scale(avgIOD_prebreeding),
         year = as.numeric(as.character(year))) %>%
  filter(!is.na(breeding_success))



## How does the probability of breeding success change with climate ------------

### SAM ========================================================================

#### FEMALES -------------------------------------------------------------------

f_SAM_glmm <- glmmTMB(breeding_success ~ age_s*avgSAM_prebreeding_s*boldness_s +
                   I(age_s^2)*avgSAM_prebreeding_s*boldness_s +
                   age_s*avgSAM_breeding_s*boldness_s +
                   I(age_s^2)*avgSAM_breeding_s*boldness_s +
                   (1|year) + (1|ring),
                 data =  female_rs,
                 family = "binomial")

summary(f_SAM_glmm)

## Remove non-significant random effects
f_SAM_glmm <- update(f_SAM_glmm, ~ age_s + I(age_s^2) + avgSAM_prebreeding_s +
                           avgSAM_breeding_s + boldness_s + 
                           I(age_s^2):avgSAM_breeding_s + (1|year) + (1|ring))
summary(f_SAM_glmm)

f_SAM_glmm <- update(f_SAM_glmm, ~ age_s + I(age_s^2) + avgSAM_prebreeding_s +
                       avgSAM_breeding_s + boldness_s + (1|year) + (1|ring))

summary(f_SAM_glmm)
tab_model(f_SAM_glmm, show.stat = T)


## Get residuals
female_rs$residuals <- resid(f_SAM_glmm)

f_SAM_pred <- data.frame(ggpredict(f_SAM_glmm, terms = c("age_s [all]")))
f_SAM_pred$age <- (f_SAM_pred$x*sd(female_rs$Age)) + mean(female_rs$Age)



#### MALES ---------------------------------------------------------------------
m_SAM_glmm <- glmmTMB(breeding_success ~ age_s*avgSAM_prebreeding_s*boldness_s +
                        I(age_s^2)*avgSAM_prebreeding_s*boldness_s +
                        age_s*avgSAM_breeding_s*boldness_s +
                        I(age_s^2)*avgSAM_breeding_s*boldness_s +
                        (1|year) + (1|ring),
                      data =  male_rs,
                      family = "binomial")

summary(m_SAM_glmm)

## Drop non-significant interactions
m_SAM_glmm <- update(m_SAM_glmm, ~ age_s + I(age_s^2) + avgSAM_prebreeding_s +
                       avgSAM_breeding_s + boldness_s + (1|year) + (1|ring))

summary(m_SAM_glmm)
tab_model(m_SAM_glmm, show.stat = T)

# Get the plot data
m_SAM_pred <- data.frame(ggpredict(m_SAM_glmm, terms = c("age_s [all]")))
m_SAM_pred$age <- (m_SAM_pred$x*sd(male_rs$Age)) + mean(male_rs$Age)



### SOI ========================================================================

#### FEMALES -------------------------------------------------------------------
f_SOI_glmm <- glmmTMB(breeding_success ~ age_s*avgSOI_prebreeding_s*boldness_s +   
                   I(age_s^2)*avgSOI_prebreeding_s*boldness_s +
                   age_s*avgSOI_breeding_s*boldness_s + 
                   I(age_s^2)*avgSOI_breeding_s*boldness_s +
                   (1|year) + (1|ring), 
                 data =  female_rs, 
                 family = "binomial")

summary(f_SOI_glmm)

## Drop non-significant interactions
f_SOI_glmm <- update(f_SOI_glmm, ~ age_s +  I(age_s^2) + avgSOI_prebreeding_s + 
                          avgSOI_breeding_s + boldness_s + (1|year) + (1|ring))

summary(f_SOI_glmm)
tab_model(f_SOI_glmm, show.stat = T)


#### MALES -------------------------------------------------------------------
m_SOI_glmm <- glmmTMB(breeding_success ~ age_s*avgSOI_prebreeding_s*boldness_s +   
                   I(age_s^2)*avgSOI_prebreeding_s*boldness_s +
                   age_s*avgSOI_breeding_s*boldness_s + 
                   I(age_s^2)*avgSOI_breeding_s*boldness_s +
                  (1|year) + (1|ring), 
                 data =  male_rs, 
                 family = "binomial")

summary(m_SOI_glmm)

## Drop non-significant interactions
m_SOI_glmm <- update(m_SOI_glmm, ~ age_s +  I(age_s^2) + avgSOI_prebreeding_s + 
                       avgSOI_breeding_s + boldness_s + (1|year) + (1|ring))

summary(m_SOI_glmm)
tab_model(m_SOI_glmm, show.stat = T)



# VISUALISE RESULTS ============================================================

### Summarise population-level breeding success with climate -------------------

# Calculate mean RS per year
annual_rs <- waal_rs %>%
  group_by(year, avgSAM_breeding, avgSAM_prebreeding, avgSOI_prebreeding, avgSOI_breeding) %>%
  summarise(n_birds = n_distinct(ring),
            mean_rs = mean(breeding_success, na.rm = T)) %>%
  filter(n_birds > 3) %>%
  mutate(year = as.numeric(as.character(year)))

# Get weighted overall mean
annual_means <- waal_rs %>%
  group_by(year) %>%
  summarize(annual_rs = mean(breeding_success, na.rm = T),
            n_inds = n())

annual_means$weights <- annual_means$n_inds/sum(annual_means$n_inds)
xm <- weighted.mean(annual_means$annual_rs, annual_means$weights)
#[1] 0.7796432

weighted_var <- sum(annual_means$weights * (annual_means$annual_rs - xm)^2)
sqrt(weighted_var)
# [1] 0.03933454


### Make the plots -------------------------------------------------------------

# SAM : 
RS_sam_breeding.plot <- ggplot() + 
  geom_point(data = annual_rs, aes(x = avgSAM_breeding, y = mean_rs, 
                                   size = n_birds, col = year)) +
  scale_colour_gradient(high = "#006125", low = "#ebf0ed", name = "Year") +
  labs(x = "Mean Southern Annular Mode (January to April)",
       y = "Population-level mean breeding success") +
  scale_y_continuous(limit = c(0.65, 0.85), breaks = seq(0.7, 1, by = 0.1)) +
  scale_size_continuous(name = "n birds") +
  theme_bw() + 
  theme(legend.position = c(0.85, 0.15),
        legend.box = "horizontal",
        legend.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


RS_sam_prebreeding.plot <- ggplot() + 
  geom_point(data = annual_rs, aes(x = avgSAM_prebreeding, y = mean_rs, 
                                   size = n_birds, col = year)) +
  scale_colour_gradient(high = "#006125", low = "#ebf0ed") +
  labs(x = "Mean Southern Annular Mode (September to November)",
       y = "Population-level mean breeding success") +
  scale_y_continuous(limit = c(0.65, 0.85), breaks = seq(0.7, 1, by = 0.1)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


# SOI : 
RS_soi_breeding.plot <- ggplot() +
  geom_point(data = annual_rs, aes(x = avgSOI_breeding, y = mean_rs, 
                                   size = n_birds, col = year)) +
  scale_colour_gradient(high = "#006125", low = "#ebf0ed") +
  labs(x = "Mean Southern Oscillation Index (January to April)",
       y = "Population-level breeding success") +
  scale_y_continuous(limit = c(0.65, 0.85), breaks = seq(0.7, 1, by = 0.1)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


RS_soi_prebreeding.plot <- ggplot() +
  geom_point(data = annual_rs, aes(x = avgSOI_prebreeding, y = mean_rs, 
                                   size = n_birds, col = year)) +
  scale_colour_gradient(high = "#006125", low = "#ebf0ed") +
  labs(x = "Mean Southern Oscillation Index (September to November)",
       y = "Population-level breeding success") +
  scale_y_continuous(limit = c(0.65, 0.85), breaks = seq(0.7, 1, by = 0.1)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))


# FIGURE 6: breeding success ~ climate =========================================

png("Figures/FIG6_RS_by_climate.png", width = 12, height = 12, units = "in", res = 300)
ggpubr::ggarrange(RS_sam_breeding.plot, RS_sam_prebreeding.plot,
                  RS_soi_breeding.plot, RS_soi_prebreeding.plot,
                  ncol = 2, nrow = 2,
                  widths = c(1, 0.92))
dev.off()

