# Flexible behaviour buffers climate variability in the wandering albatross _Diomedea exulans_
Natasha Gillies, Jack Thorley, Henri Weimerskirch, Stéphanie Jenouvrier, Samantha C. Patrick
## Overview
This repository contains scripts and data to recreate the main results and figures of this paper (currently in prep).

## Scripts
A short description of each script is given below.

- **MS_foraging_by_climate_indcices.R** This script
- **MS_kernel_map.R** This script
- **MS_plot_climate_indices.R** This script
- **MS_RS_by_climate.R** This script
- 
## Data inputs 
These datasets are used in the above scripts. Note that individual IDs have been recoded and so cannot be linked to existing datasets. Please contact the authors if you would like to make sure of these, as we may be able to offer addiitonal information, data, or advice. 

- **foraging_data_climate_analyses_F.RData** This is the main dataset for the statistical analyses. Each row corresponds to an individual foraging trip. A second dataset is available for males (_foraging_data_climate_analyses_M.RData_) and has an identical structure. The variables are as follows:
	- _iD_: Encodes unique ID of bird; anonymised from original data (factor)
- **SOI_monthly.RData** This contains data on the monthly Southern Oscillation Index indices for each year of the study. The file _SAM_monthly.RData_ contains the same information for the Southern Annular Mode.
	- _BirdID_: Encodes unique ID of bird; anonymised from original data (factor)
