## Project Summary
This repository holds the code and data for designating Cases and Controls for specific enteric pathogens for the manuscript "HLA Class I and II Associations with Common Enteric Pathogens in the First Year of Life" (McCowin et al, 2021). This work is published in <i>eBioMedicine</i> and is available by clicking on the following [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8093888/). 

The goal of this analysis was to identify associations with specific HLA haplotypes and susceptibility to 12 specific enteric pathogens in the PROVIDE cohort (see Kirkpatrick et al, 2015 for more information on this population). To that aim, we first designated individuals as cases or controls for each pathogen before analyzing associations with HLA data. The process of designating individuals as cases or controls is outlined here.

## Overview of Files Included

|                 File                  |                                                    Description                                                            |
| :-----------------------------------: | :-----------------------------------------------------------------------------------------------------------------------: |
| Determination of Cases and Controls.R | Outlines the process of designating cases and controls based on qPCR data for specific pathogens from diarrheal episodes. |
| lab_bv_hla_f10.csv                    | Contains HLA data and SIDs for all PROVIDE participants tested.                                                           |
| lab_bv_tacd_f10_Supp_r2.csv           | Contains TAC card data for the 12 pathogens analyzed in this study for all diarrheal samples.                             |
| mgmt_rto_rotatrial_outcome_f10.csv    | Contains whether each individual in PROVIDE was dropped within the first year of life or not.                             |

## Dependancies
This code requires only the tidyverse package (Wickham et al, 2019) in R to run.

## Contact Information

G. Brett Moreau - gbm5pn@virginia.edu
