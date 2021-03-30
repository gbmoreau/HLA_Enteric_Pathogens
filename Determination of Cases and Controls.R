#### Grouping PROVIDE Individuals as Cases or Controls Based on Diarrheal TAC Data #######################################################################################
# NAME: G. Brett Moreau
# DATE: October 30, 2020

#### PACKAGES ############################################################################################################################################################
#install.packages("tidyverse")
library(tidyverse)


#### INTRODUCTION ########################################################################################################################################################
# We are interested in identifying associations between HLA alleles/haplotypes and incidence of enteric pathogens in the PROVIDE 
# cohort There are 12 pathogens in particular that we are interested in (EAEC, Adenovirus 40/41, LT-ETEC, Rotavirus, ST-ETEC, 
# Typical EPEC, C.jejuni/coli, Sapovirus, Norovirus GII, Astrovirus, EIEC&Shigella, and Cryptosporidium). I'll be processing the 
# data, organizing each of the individuals for whom we have HLA data into cases or controls for each of these pathogens based on 
# if they ever had diarrheal samples that were PCR positive during the first year of life.


#### COLLECT SIDS WITH HLA DATA ###########################################################################################################################################
# I'm only interested in individuals with HLA data, so I'll first collect the SIDs for individuals with HLA data. 

HLA <- read.csv("lab_bv_hla_f10.csv") # This is the full HLA data set. 
HLA.SIDs <- select(HLA, sid) # This filters for only the SIDs of individuals with HLA data.

length(unique(HLA.SIDs$sid)) # There are 640 individuals with HLA data.


# There are 640 individuals with HLA data, but not all of them were followed for the entire first year of life. I'll be pulling 
# whether they were dropped within the first 365 days from the PROVIDE file "mgmt_rto_rotatrial_outcome_f10", which lists this 
# information.


all <- read.csv("mgmt_rto_rotatrial_outcome_f10.csv") # This should include all 700 individuals in the PROVIDE cohort.

followed <- filter(all, drop365 == 0) 
length(unique(followed$sid)) # This should include only those individuals who were followed for the entire first year 
# of life. The total number of individuals followed for the first year of life is 611.


### SIDs ONLY ###
all.SIDs <- select(all, sid)
followed.SIDs <- select(followed, sid)


# Finally, I'll combine the SIDs for individuals with HLA data and the SIDs for individuals who were followed for the entire 
# first year of life, keeping only those SIDs that appear in both data sets.


HLA.followed <- inner_join(HLA.SIDs, followed.SIDs, by = "sid") 
length(unique(HLA.followed$sid)) # There are 601 individuals who were followed for the entire first year of life and also have 
# HLA info. 



#### FILTERING TAC DATA FOR EPISODES WITHIN THE FIRST YEAR OF LIFE ###############################################################
# Cases and controls are designated based on the diarrheal TAC data, which has been edited by Mami. I'll be using that data 
# set for analysis.

TAC.original <- read.csv("lab_bv_tacd_f10_Supp_r2.csv")

# Mami's TAC data set had spaces in between each episode, which introduced many "NA"s into the data set. To remove them, I'll 
# merge this data set with the file containing all PROVIDE SIDs, keeping only the information with PROVIDE SIDs.
TAC <- left_join(all.SIDs, TAC.original, by = "sid")


# This also introduces individuals who are not found in the TAC data set as SIDs with no episode information. I can collect 
# these SIDs to identify individuals without any TAC information, who will be designated as controls for all pathogens.

TAC.nodata <- filter(TAC, is.na(specdt)) # Filters only for individuals with no specimen collection data (i.e., individuals 
# without TAC data).
TAC.nodata <- select(TAC.nodata, sid) # Collects only the SIDs.


# Then I'll remove these SIDs, so I have a TAC data set that includes only individuals with TAC data.

TAC <- anti_join(TAC, TAC.nodata, by = "sid") # 130 observations were removed, corresponding to the 130 individuals without 
# TAC data.


# These data sets still include individuals without HLA data or who were not followed for the full first year of life. So 
# I'll filter out individuals who didn't meet these criteria. There are 601 individuals who meet these criteria, so the 
# total for these two data sets should add up to 601.

TAC <- inner_join(TAC, HLA.followed, by = "sid") # Filters for only individuals with no TAC data AND who were followed for 
# 1 year and have HLA data.
length(unique(TAC$sid)) # There are 539 individuals who were followed out to 1 year with HLA and TAC data.


TAC.nodata <- inner_join(TAC.nodata, HLA.followed, by = "sid") # Filters for only individuals with no TAC data AND who 
# were followed for 1 year and have HLA data.
length(unique(TAC.nodata$sid)) # There are 62 individuals who were followed out to 1 year with HLA and no TAC data. 
# This adds up to a total of 601 individuals.


# We are only interested in diarrheal episodes that occurred during the first year of life. Therefore, I'll filter to 
# remove any diarrheal episodes that occurred after the first year of life. This data is indicated by the "episode.1yr" column,
# which marks whether an specific diarrheal episode occurred within the first year of life ("1") or not ("0").

TAC.1year <- filter(TAC, episode.1yr == 1) 
length(unique(TAC.1year$sid)) # Overall, 208 diarrheal episodes were excluded for being outside the first year 
# of life, giving a total of 1917 diarrheal episodes for 529 unique individuals.


# There are 10 individuals with TAC data, but not in the first year of life. These SIDs will be collected and 
# added to the TAC.nodata data set, since they will also count as controls for our purposes.


TAC.nodata.1year <- anti_join(TAC, TAC.1year, by = "sid")
TAC.nodata.1year <- TAC.nodata.1year %>%
  group_by(sid) %>%
  arrange(specdt) %>%
  slice(1L) %>%
  select(sid) # This will only collect the unique SIDs

TAC.nodata.1year <- rbind(TAC.nodata, TAC.nodata.1year) # This will combine the data set of individuals with 
# no TAC data and the data set of individuals with TAC data only in the second year of life.


#### DETERMINING WHETHER AN INDIVIDUAL IS POSITIVE OR NEGATIVE FOR A GIVEN PATHOGEN ####################################################################################
# Finally, I need to determine whether each individual is positive for each pathogen. Individuals who ever 
# had a positive instance with a specific pathogen will be designated as cases, and those who did not have 
# a positive instance with a specific pathogen will be designated as controls. Positive instances are 
# defined as a TAC reading of less than 35, and controls as a TAC reading of greater than or equal to 35.

# There are 12 pathogens in which we're interested. Of the 12, 10 of them are based on TAC data from a 
# single column. The other two (EAEC and ST-ETEC) are based on two columns, so that a positive TAC instance 
# in either column counts designates an individual as a case. 

# The TAC data used for each pathogen is specified under each header.


#### ADENOVIRUS ########################################################################################################################################################
# Adenovirus cases are defined as a positive TAC result from the "adeno40" column of the TAC data set, 
# which measures Adenovirus F.


# First, I'll make a new column which defines individuals with a TAC value of >35 as cases and <=35 as controls.
TAC.1year$adeno40.positive[TAC.1year$adeno40 < 35 & TAC.1year$adeno40 > 0] <- 1 
TAC.1year$adeno40.positive[TAC.1year$adeno40 >= 35] <- 0 

# Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$adeno40.positive[TAC.1year$adeno40 == -9] <- 0 

table(TAC.1year$adeno40.positive) # All episodes are accounted for. There are 1282 negative and 635 positive episodes.


# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "adeno40.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative)
# and the "freq" column counts how many times that outcome was observed per individual. 


adeno40 <- plyr::count(TAC.1year, c("sid", "adeno40.positive")) 
# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need to 
# get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the positive 
# column in descending order, then take only one slice per individual. By doing this, I'll take only the positive instance 
# for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the only number 
# for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

adeno40 <- adeno40 %>%
  group_by(sid) %>%
  arrange(desc(adeno40.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, adeno40.positive)

length(unique(adeno40$sid)) # There are 529 individuals in the data set, as expected.


# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.adeno40 <- TAC.nodata.1year
TAC.nodata.adeno40$adeno40.positive <- 0 # All are designated as control.

adeno40 <- rbind(adeno40, TAC.nodata.adeno40)

length(unique(adeno40$sid)) # The total number of individuals is 601, as expected.


table(adeno40$adeno40.positive) # For Adeno40, there are 251 individuals designated as controls and 350 individuals designated 
# as cases.

write.table(adeno40, file = "adeno40.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### EAEC ######################################################################################################################################################
# EAEC cases are defined as a positive TAC result from EITHER the "eaecic" column OR the "eaecat" column. I'll first start 
# by counting "eaecic" positive instances, then counting "eaecat" positive instances, and then combining this data to 
# identify controls and cases.

### EAECIC ###
#I'll make a new column which defines individuals with a TAC value of >35 as cases and <=35 as controls. Individuals with 
# "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$eaecic.positive[TAC.1year$eaecic < 35 & TAC.1year$eaecic > 0] <- 1 
TAC.1year$eaecic.positive[TAC.1year$eaecic >= 35] <- 0 
TAC.1year$eaecic.positive[TAC.1year$eaecic == -9] <- 0 

table(TAC.1year$eaecic.positive) 

# All 1917 episodes are accounted for. There are 1269 negative and 648 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "eaecic.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

eaecic <- plyr::count(TAC.1year, c("sid", "eaecic.positive")) 

eaecic <- eaecic %>%
  group_by(sid) %>%
  arrange(desc(eaecic.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, eaecic.positive)

length(unique(eaecic$sid)) 

# There are 529 individuals in the data set, as expected.

### EAECAT ###
#I'll make a new column which defines individuals with a TAC value of >35 as cases and <=35 as controls. Individuals with 
# "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$eaecat.positive[TAC.1year$eaecat < 35 & TAC.1year$eaecat > 0] <- 1 
TAC.1year$eaecat.positive[TAC.1year$eaecat >= 35] <- 0 
TAC.1year$eaecat.positive[TAC.1year$eaecat == -9] <- 0 

table(TAC.1year$eaecat.positive) 

# All 1917 episodes are accounted for. There are 901 negative and 1016 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "eaecat.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

eaecat <- plyr::count(TAC.1year, c("sid", "eaecat.positive")) 

eaecat <- eaecat %>%
  group_by(sid) %>%
  arrange(desc(eaecat.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, eaecat.positive)

length(unique(eaecat$sid)) 

# There are 529 individuals in the data set, as expected.


# Now, to get an actual count of cases and controls, I'll add together the counts of whether an individual was positive for 
# eaecic and the counts of whether an individual was positive for eaecat. If an individual is marked as negative for both, 
# they will be designated as a control, whereas if they're marked as positive for either TAC marker or both they'll be 
# designated as a case.

EAEC <- full_join(eaecic, eaecat, by = "sid") # Combine eaecic and eaecat positive counts.
EAEC$EAEC.positive <- NA # Create column for total positive or negative count.
EAEC$EAEC.positive[EAEC$eaecic.positive == 1] <- 1 # Indicate individuals who were eaecic positive as cases.
EAEC$EAEC.positive[EAEC$eaecat.positive == 1] <- 1 # Indicate individuals who were eaecat positive as cases.
EAEC$EAEC.positive[EAEC$eaecic.positive == 0 & EAEC$eaecat.positive == 0] <- 0 # Indicate individuals who were negative for 
# both eaecic and eaecat as controls.

table(EAEC$EAEC.positive)

# Overall, of the individuals with TAC data, there were 73 controls and 456 cases. This adds up to 529 individuals, so 
# everyone is accounted for.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.EAEC <- TAC.nodata.1year
TAC.nodata.EAEC$EAEC.positive <- 0 # All are designated as control.

EAEC <- rbind(EAEC, TAC.nodata.EAEC)

length(unique(EAEC$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number 
# of individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(EAEC$EAEC.positive) 

# Overall for EAEC, there are 145 individuals designated as controls and 456 individuals designated as cases.

write.table(EAEC, file = "EAEC.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### ST-ETEC ###############################################################################################################################################################
# ST-ETEC cases are defined as a positive TAC result from EITHER the "etecsth" column OR the "etecstp" column. I'll first 
# start by counting "etecsth" positive instances, then counting "etecstp" positive instances, and then combining this data 
# to identify controls and cases.

### ETECSTH ###
#I'll make a new column which defines individuals with a TAC value of >35 as cases and <=35 as controls. Individuals with 
# "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$etecsth.positive[TAC.1year$etecsth < 35 & TAC.1year$etecsth > 0] <- 1 
TAC.1year$etecsth.positive[TAC.1year$etecsth >= 35] <- 0 
TAC.1year$etecsth.positive[TAC.1year$etecsth == -9] <- 0 

table(TAC.1year$etecsth.positive) 

# All 1917 episodes are accounted for. There are 1446 negative and 471 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "etecsth.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

etecsth <- plyr::count(TAC.1year, c("sid", "etecsth.positive")) 

etecsth <- etecsth %>%
  group_by(sid) %>%
  arrange(desc(etecsth.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, etecsth.positive)

length(unique(etecsth$sid)) 

# There are 529 individuals in the data set, as expected.

### ETECSTP ###
#I'll make a new column which defines individuals with a TAC value of >35 as cases and <=35 as controls. Individuals with 
# "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$etecstp.positive[TAC.1year$etecstp < 35 & TAC.1year$etecstp > 0] <- 1 
TAC.1year$etecstp.positive[TAC.1year$etecstp >= 35] <- 0 
TAC.1year$etecstp.positive[TAC.1year$etecstp == -9] <- 0 

table(TAC.1year$etecstp.positive) 

# All 1917 episodes are accounted for. There are 1771 negative and 146 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "etecstp.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

etecstp <- plyr::count(TAC.1year, c("sid", "etecstp.positive")) 

etecstp <- etecstp %>%
  group_by(sid) %>%
  arrange(desc(etecstp.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, etecstp.positive)

length(unique(etecstp$sid)) 

# There are 529 individuals in the data set, as expected.

# Now, to get an actual count of cases and controls, I'll add together the counts of whether an individual was positive for 
# etecsth and the counts of whether an individual was positive for etecstp. If an individual is marked as negative for both, 
# they will be designated as a control, whereas if they're marked as positive for either TAC marker or both they'll be 
# designated as a case.

etecst <- full_join(etecsth, etecstp, by = "sid") # Combine etecsth and etecstp positive counts.
etecst$etecst.positive <- NA # Create column for total positive or negative count.
etecst$etecst.positive[etecst$etecsth.positive == 1] <- 1 # Indicate individuals who were etecsth positive as cases.
etecst$etecst.positive[etecst$etecstp.positive == 1] <- 1 # Indicate individuals who were etecstp positive as cases.
etecst$etecst.positive[etecst$etecsth.positive == 0 & etecst$etecstp.positive == 0] <- 0 # Indicate individuals who were 
# negative for both etecsth and etecstp as controls.

table(etecst$etecst.positive)

# Overall, of the individuals with TAC data, there were 211 controls and 318 cases. This adds up to 529 individuals, so 
# everyone is accounted for.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.etecst <- TAC.nodata.1year
TAC.nodata.etecst$etecst.positive <- 0 # All are designated as control.

etecst <- rbind(etecst, TAC.nodata.etecst)

length(unique(etecst$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(etecst$etecst.positive) 

# Overall for etecst, there are 283 individuals designated as controls and 318 individuals designated as cases.

write.table(etecst, file = "etecst.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### C.JEJUNI/COLI ######################################################################################################################################################
# C. jejuni/coli cases are defined as a positive TAC result from the "cjeju" column. 

# First, I'll start with "cjeju" TAC column. I'll make a new column which defines individuals with a TAC value of 
# >35 as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered 
# controls.
TAC.1year$cjeju.positive[TAC.1year$cjeju < 35 & TAC.1year$cjeju > 0] <- 1 
TAC.1year$cjeju.positive[TAC.1year$cjeju >= 35] <- 0 
TAC.1year$cjeju.positive[TAC.1year$cjeju == -9] <- 0 

table(TAC.1year$cjeju.positive) 

# All 1917 episodes are accounted for. There are 1450 negative and 467 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "cjeju.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

cjeju <- plyr::count(TAC.1year, c("sid", "cjeju.positive")) 

cjeju <- cjeju %>%
  group_by(sid) %>%
  arrange(desc(cjeju.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, cjeju.positive)

length(unique(cjeju$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.cjeju <- TAC.nodata.1year
TAC.nodata.cjeju$cjeju.positive <- 0 # All are designated as control.

cjeju <- rbind(cjeju, TAC.nodata.cjeju)

length(unique(cjeju$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(cjeju$cjeju.positive) 

# For cjejuvirus, there are 335 individuals designated as controls and 266 individuals designated as cases.

write.table(cjeju, file = "cjeju.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### NOROVIRUS G2 ##########################################################################################################################################################
# Norovirus G2 cases are defined as a positive TAC result from the "norog2" column. 

# First, I'll start with "norog2" TAC column. I'll make a new column which defines individuals with a TAC value of 
# >35 as cases and <=35 as controls.
TAC.1year$norog2.positive[TAC.1year$norog2 < 35 & TAC.1year$norog2 > 0] <- 1 
TAC.1year$norog2.positive[TAC.1year$norog2 >= 35] <- 0 

# Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$norog2.positive[TAC.1year$norog2 == -9] <- 0 

table(TAC.1year$norog2.positive) # All 1917 episodes are accounted for. There are 1517 negative and 400 positive episodes.


# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "norog2.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative)
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

norog2 <- plyr::count(TAC.1year, c("sid", "norog2.positive"))

norog2 <- norog2 %>%
  group_by(sid) %>%
  arrange(desc(norog2.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, norog2.positive)

length(unique(norog2$sid)) # There are 529 individuals in the data set, as expected.


# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.norog2 <- TAC.nodata.1year
TAC.nodata.norog2$norog2.positive <- 0 # All are designated as control.

norog2 <- rbind(norog2, TAC.nodata.norog2)

length(unique(norog2$sid)) # The total number of individuals is 601, as expected.


table(norog2$norog2.positive) # For Norovirus G2, there are 337 individuals designated as controls and 264 individuals 
# designated as cases.

write.table(norog2, file = "norog2.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### SAPOVIRUS #########################################################################################################################################################
# Sapovirus cases are defined as a positive TAC result from the "sapo" column. 

# First, I'll start with "sapo" TAC column. I'll make a new column which defines individuals with a TAC value of >35 
# as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.

TAC.1year$sapo.positive[TAC.1year$sapo < 35 & TAC.1year$sapo > 0] <- 1 
TAC.1year$sapo.positive[TAC.1year$sapo >= 35] <- 0 
TAC.1year$sapo.positive[TAC.1year$sapo == -9] <- 0 

table(TAC.1year$sapo.positive) 

# All 1917 episodes are accounted for. There are 1526 negative and 391 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "sapo.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

sapo <- plyr::count(TAC.1year, c("sid", "sapo.positive")) 

sapo <- sapo %>%
  group_by(sid) %>%
  arrange(desc(sapo.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, sapo.positive)

length(unique(sapo$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.sapo <- TAC.nodata.1year
TAC.nodata.sapo$sapo.positive <- 0 # All are designated as control.

sapo <- rbind(sapo, TAC.nodata.sapo)

length(unique(sapo$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(sapo$sapo.positive) 

# For Sapovirus, there are 344 individuals designated as controls and 257 individuals designated as cases.

write.table(sapo, file = "sapo.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### ROTAVIRUS ######################################################################################################################################################
# Rotavirus cases are defined as a positive TAC result from the "rota" column. 

# First, I'll start with "rota" TAC column. I'll make a new column which defines individuals with a TAC value of >35 
# as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$rota.positive[TAC.1year$rota < 35 & TAC.1year$rota > 0] <- 1 
TAC.1year$rota.positive[TAC.1year$rota >= 35] <- 0 
TAC.1year$rota.positive[TAC.1year$rota == -9] <- 0 

table(TAC.1year$rota.positive) 

# All 1917 episodes are accounted for. There are 1558 negative and 359 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
#  individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "rota.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

rota <- plyr::count(TAC.1year, c("sid", "rota.positive")) 

rota <- rota %>%
  group_by(sid) %>%
  arrange(desc(rota.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, rota.positive)

length(unique(rota$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.rota <- TAC.nodata.1year
TAC.nodata.rota$rota.positive <- 0 # All are designated as control.

rota <- rbind(rota, TAC.nodata.rota)

length(unique(rota$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(rota$rota.positive) 

# For rotavirus, there are 329 individuals designated as controls and 272 individuals designated as cases.

write.table(rota, file = "rota.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### ASTROVIRUS ######################################################################################################################################################
# Astrovirus cases are defined as a positive TAC result from the "astro" column. 

# First, I'll start with "astro" TAC column. I'll make a new column which defines individuals with a TAC value of >35 
# as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$astro.positive[TAC.1year$astro < 35 & TAC.1year$astro > 0] <- 1 
TAC.1year$astro.positive[TAC.1year$astro >= 35] <- 0 
TAC.1year$astro.positive[TAC.1year$astro == -9] <- 0 

table(TAC.1year$astro.positive) 

# All 1917 episodes are accounted for. There are 1644 negative and 273 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "astro.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need to 
# get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the positive 
# column in descending order, then take only one slice per individual. By doing this, I'll take only the positive instance 
# for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the only number 
# for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

astro <- plyr::count(TAC.1year, c("sid", "astro.positive")) 

astro <- astro %>%
  group_by(sid) %>%
  arrange(desc(astro.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, astro.positive)

length(unique(astro$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.astro <- TAC.nodata.1year
TAC.nodata.astro$astro.positive <- 0 # All are designated as control.

astro <- rbind(astro, TAC.nodata.astro)

length(unique(astro$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(astro$astro.positive) 

# For astrovirus, there are 394 individuals designated as controls and 207 individuals designated as cases.

write.table(astro, file = "astro.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### LT-ETEC ######################################################################################################################################################
# LT-ETEC cases are defined as a positive TAC result from the "eteclt" column. 

# First, I'll start with "eteclt" TAC column. I'll make a new column which defines individuals with a TAC value of >35 
# as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$eteclt.positive[TAC.1year$eteclt < 35 & TAC.1year$eteclt > 0] <- 1 
TAC.1year$eteclt.positive[TAC.1year$eteclt >= 35] <- 0 
TAC.1year$eteclt.positive[TAC.1year$eteclt == -9] <- 0 

table(TAC.1year$eteclt.positive) 

# All 1917 episodes are accounted for. There are 1255 negative and 662 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "eteclt.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need to 
# get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the positive 
# column in descending order, then take only one slice per individual. By doing this, I'll take only the positive instance 
# for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the only number 
# for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

eteclt <- plyr::count(TAC.1year, c("sid", "eteclt.positive")) 

eteclt <- eteclt %>%
  group_by(sid) %>%
  arrange(desc(eteclt.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, eteclt.positive)

length(unique(eteclt$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.eteclt <- TAC.nodata.1year
TAC.nodata.eteclt$eteclt.positive <- 0 # All are designated as control.

eteclt <- rbind(eteclt, TAC.nodata.eteclt)

length(unique(eteclt$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(eteclt$eteclt.positive) 

# For etecltvirus, there are 261 individuals designated as controls and 340 individuals designated as cases.

write.table(eteclt, file = "eteclt.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#### SHIGELLA AND EITC ######################################################################################################################################################
# Shigella and EITC cases are defined as a positive TAC result from the "eiecsh" column. 

# First, I'll start with "eiecsh" TAC column. I'll make a new column which defines individuals with a TAC value of >35 
# as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$eiecsh.positive[TAC.1year$eiecsh < 35 & TAC.1year$eiecsh > 0] <- 1 
TAC.1year$eiecsh.positive[TAC.1year$eiecsh >= 35] <- 0 
TAC.1year$eiecsh.positive[TAC.1year$eiecsh == -9] <- 0 

table(TAC.1year$eiecsh.positive) 

# All 1917 episodes are accounted for. There are 1723 negative and 194 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two different 
# rows. The "eiecsh.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) and the 
# "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

eiecsh <- plyr::count(TAC.1year, c("sid", "eiecsh.positive")) 

eiecsh <- eiecsh %>%
  group_by(sid) %>%
  arrange(desc(eiecsh.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, eiecsh.positive)

length(unique(eiecsh$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.eiecsh <- TAC.nodata.1year
TAC.nodata.eiecsh$eiecsh.positive <- 0 # All are designated as control.

eiecsh <- rbind(eiecsh, TAC.nodata.eiecsh)

length(unique(eiecsh$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(eiecsh$eiecsh.positive) 

# For Shigella and EITC, there are 458 individuals designated as controls and 143 individuals designated as cases.

write.table(eiecsh, file = "eiecsh.txt", sep = "\t", row.names = FALSE, quote = FALSE)




#### CRYPTOSPORIDIUM ######################################################################################################################################################
# Cryptosporidium cases are defined as a positive TAC result from the "cryp" column. 

# First, I'll start with "cryp" TAC column. I'll make a new column which defines individuals with a TAC value of >35 
# as cases and <=35 as controls. Individuals with "Missing Values" (-9) for their TAC data will also be considered controls.
TAC.1year$cryp.positive[TAC.1year$cryp < 35 & TAC.1year$cryp > 0] <- 1 
TAC.1year$cryp.positive[TAC.1year$cryp >= 35] <- 0 
TAC.1year$cryp.positive[TAC.1year$cryp == -9] <- 0 

table(TAC.1year$cryp.positive) 

# All 1917 episodes are accounted for. There are 1760 negative and 157 positive episodes.

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "cryp.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual. 

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the 
# positive column in descending order, then take only one slice per individual. By doing this, I'll take only the positive 
# instance for cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the 
# only number for these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.

cryp <- plyr::count(TAC.1year, c("sid", "cryp.positive")) 

cryp <- cryp %>%
  group_by(sid) %>%
  arrange(desc(cryp.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, cryp.positive)

length(unique(cryp$sid)) 

# There are 529 individuals in the data set, as expected.

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.

TAC.nodata.cryp <- TAC.nodata.1year
TAC.nodata.cryp$cryp.positive <- 0 # All are designated as control.

cryp <- rbind(cryp, TAC.nodata.cryp)

length(unique(cryp$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.


table(cryp$cryp.positive) 

# For Cryptosporidium, there are 499 individuals designated as controls and 102 individuals designated as cases.

write.table(cryp, file = "cryp.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#### TYPICAL EPEC ##################################################################################################################
# Typical EPEC cases are defined as a positive TAC result from the "epecbf" column. First, I'll start with "epecbf" TAC column. 
# I'll make a new column which defines individuals with a TAC value of >35 as cases and <=35 as controls. Individuals with 
# "Missing Values" (-9) for their TAC data will also be considered controls.  

TAC.1year$epecbf.positive[TAC.1year$epecbf < 35 & TAC.1year$epecbf > 0] <- 1 
TAC.1year$epecbf.positive[TAC.1year$epecbf >= 35] <- 0 
TAC.1year$epecbf.positive[TAC.1year$epecbf == -9] <- 0 

table(TAC.1year$epecbf.positive) 

# All 1917 episodes are accounted for. There are 1417 negative and 500 positive episodes.  

# Next, I'll perform a count on this column, which will count the number of positive or negative diarrheal episodes per 
# individual. This produces a data set where individuals with both positive and negative episodes will appear as two 
# different rows. The "epecbf.positive" column lists the potential outcomes for this column (1 == positive; 0 == negative) 
# and the "freq" column counts how many times that outcome was observed per individual.  

# Because this data set marks individuals with both positive and negative episodes as both positive and negative, I need 
# to get rid of the negative episode, marking them as cases. To do this, I'll first group by SID, then arrange by the positive 
# column in descending order, then take only one slice per individual. By doing this, I'll take only the positive instance for 
# cases (because "1" is the higher number) but take the negative instance for all others (because "0" is the only number for 
#these individuals). There are 529 unique individuals, so I should get a data set with 529 rows.  

epecbf <- plyr::count(TAC.1year, c("sid", "epecbf.positive")) 

epecbf <- epecbf %>%
  group_by(sid) %>%
  arrange(desc(epecbf.positive)) %>% # Arranges in descending order
  slice(1L) %>%
  select(sid, epecbf.positive)

length(unique(epecbf$sid)) 

# There are 529 individuals in the data set, as expected.  

# Finally, I need to add individuals without TAC data as controls. These individuals were calculated above (see "TAC.nodata" 
# data set). There are 601 individuals who were followed for the entire first year of life and also have HLA info, so the 
# final data set should have 601 unique SIDs.  

TAC.nodata.epecbf <- TAC.nodata.1year
TAC.nodata.epecbf$epecbf.positive <- 0 # All are designated as control.

epecbf <- rbind(epecbf, TAC.nodata.epecbf)

length(unique(epecbf$sid))

# There are a total of 601 individuals who were followed for the first year of life and have HLA data. The total number of 
# individuals in our cases/controls data set is 601, meaning we haven't missed anyone.  


table(epecbf$epecbf.positive) 


# For Typical EPEC, there are 289 individuals designated as controls and 312 individuals designated as cases.  

write.table(epecbf, file = "epecbf.txt", sep = "\t", row.names = FALSE, quote = FALSE)
