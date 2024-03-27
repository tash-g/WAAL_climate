# Plastic behaviour buffers climate variability in the wandering albatross _Diomedea exulans_
Natasha Gillies, Jack Thorley, Henri Weimerskirch, Stéphanie Jenouvrier, Samantha C. Patrick
## Overview
This repository contains scripts and data to recreate the main results and figures of this paper (currently in submission).

## Scripts
A short description of each script is given below.

- **MS_RS_by_climate.R** Contains analyses to examine the effect of climate variation on reproductive success. 
- **MS_foraging_by_climate_indices.R** Contains analyses to examine the effect of climate variation on foraging behaviour.
- **MS_kernel_map.R** Calculates distribution kernels and plots them. Note this code depends on sp packages, which are not supported by recent R distributions.
- **MS_plot_climate_indices.R** Plots changes in climate indices over time. 

 
## Data inputs 
These datasets are used in the above scripts. Note that individual IDs have been recoded and so cannot be linked to existing datasets. Please contact the authors if you would like to make sure of these, as we may be able to offer addiitonal information, data, or advice. 

- **SAM_monthly.RData** This contains data on the monthly Southern Annular Mode indices for each year of the study. The file _SOI_monthly.RData_ contains the same information for the Southern Oscillation Index. Each row contains a monthly value. Data columns are as follows:
  - _Year_: Numeric year of record
  - _Month_: Three-letter month code of record
  - _SAM/SOIIndex_: Numeric value of index
  - _month_: Numeric value for month
  - _counter_: Numeric counter for dataset
  - _cat_: Factor indicating whether the value is negative (neg) or positive (pos)
  - _moving_avg_: Numeric moving average for index, rolling window = 6
    
- **WAAL_breedingSuccess_1965-2020.csv** Dataset used in reproductive success analyses. Each record represents an individual breeding attempt and outcomes for that year. Data columns are as follows:
  - _id_: Factor encoding unique ID of bird; anonymised from original 
  - _Sex_: Factor encoding whether bird is female (F) or male (M)
  - _StatusBaguage_: Factor indicating whether age data available (P) or not (A)
  - _year_: Factor encoding year of record
  - _Age_: Numeric variable encoding bird age
  - _StateCode_: Factor encoding bird's breeding state, where 0 = failed to breed, 1/2 = successfully bred, 3 = did not attempt to breed
  - _AFR_: Numeric age at first reproductive attempt
  - _boldness_BLUP_mean_: Numeric boldness value for individual
 
- **WAAL_foraging_2010-2020_F/M.csv** This is the main dataset for the statistical analyses. Each row corresponds to an individual foraging trip. 'F/M' indicates female and male dataset respectively. Data columns are as follows:
  - _id_: Factor encoding unique ID of bird; anonymised from original data
  - _Year_: Factor encoding year of record
  - _Age_: Numeric age of bird
  - _trackid_: Factor encoding unique ID of foraging trip
  - _start_time_: Datetime variable indicating start time of foraging trip
  - _DeploymentID_: Factor encoding unique ID of each deployment (per bird and year)
  - _sex_: Factor encoding sex of bird; female or male
  - _tripduration.days_: Numeric duration of foraging trip duration in days
  - _maxdistance.km_: Numeric maximum distance travelled from colony, km
  - _totalpathdistance.km_: Numeric total distance covered in foraging trip, km
  - _total_landings_hmm_: Numeric number of landings per day, determined from HMM
  - _propRest_hmm_: Numeric proportion of trip spent in rest behaviour, calculated from HMM-determined behavioural categorisation
  - _propSearch_hmm_: Numeric proportion of trip spent in search behaviour, calculated from HMM-determined behavioural categorisation 
  - _propTravel_hmm_: Numeric proportion of trip spent in travel behaviour, calculated from HMM-determined behavioural categorisation
  - _propARSvsTravel_hmm_: Numeric proportion of trip spent in search/ARS behaviour relatively to travel behaviour, calculated from HMM-determined behavioural categorisation
  - _logmaxdistance.km_: Numeric logged value of maximum distance travelled from colony, km
  - _logtotalpathdistance.km_: Numeric logged value of total distance covered in foraging trip, km
  - _southLat_: Numeric value indicating most southerly latitude reached during trip, degrees
- _northlat_: Numeric value indicating most northerly latitude reached during trip, degrees
  - _latRange_: Numeric value indicating total latitudinal range
  - _SAMIndex_: Numeric value indicating SAM value for the month in which foraging trip took place
  - _SOIIndex_: Numeric value indicating SAM value for the month in which foraging trip took place
  - _boldness_: Numeric value indicating boldness of bird tracked
  - _attempted_breeding_: Factor indicating whether bird attempted breeding (1) or not (0)
  - _breeding_success_: Factor indicating whether bird was succcessful (1) or not (0) in breeding attempt
  - _prevyear_: Factor indicating whether bird's breeding success in the previous year; failedrep = attempted but failed to reproduce, no attempt = did not attempt to breed, successful rep = successfully reproduced
 
 - **WAAL_gpsLocations_1989-2020.csv** Dataset providing GPS fixes for birds tracked between 1989 and 2020. Each row gives a GPS fix for a different indiviudal. Data columns are as follows:
   - _ring_: Factor encoding unique ID of bird; anonymised from original data
   - _DateTime_: Datetime variable giving date and time GPS fix was taken
   - _Sex_: Factor indicating sex of bird, F (female) or M (male)
   - _Year_: Factor giving year of GPS fix
   - _Longitude_: Numeric value indicating longitude of fix, degrees
   - _Latitude_: Numeric value indicating latitude of fix, degrees
   
