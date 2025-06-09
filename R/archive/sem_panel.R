### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# This script prepares and compares structural equation models
# built to test whether diversity affects complexity or vice versa
# helpful weblink: <https://www.researchgate.net/post/How_to_test_the_bidirectional_relationship_in_a_SEM_model>

#### This script uses panel-level data, 
# and includes responses from different sampling dates (30,60,90 days) in the study
# It presents a single panel model with data from multiple times
# see Grace and Larson leafy spurge model
# see also Whalen et al. 2013 Supplementary Material

# packages
library(tidyverse)
library(readxl)
library(lavaan)
library(piecewiseSEM)

# # code chunk below prepare data from other scripts in the project ----- 
# d, richness is from point counts, Shannon diversity is for genus/species level for point counts
d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
names(d) <- tolower(names(d))
d$age <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
d$panel <- gsub( "90D","90d", d$panel)
d$panel <- gsub( "60D","60d", d$panel)

# transforming rugosity measurements
d <- d %>% mutate( rugosity_raw = rugosity, rug1 = 1/rugosity_raw, rug2 = 1-rugosity_raw) %>% 
  mutate( age_factor = factor(age))

# log transform rugosity
d$logrug <- log( d$rug2 )

# community data to grab open space and arborescent bryozoans
comm_raw <- read_xlsx("data/PCover_taxgroups.xlsx")
comm_select <- comm_raw %>% select(panel = Panel, site = Site, age = Age, ar_bryo, col_asc, sol_asc, sabellids, sponge, open_space)
comm_select$age <- as.numeric(gsub("([0-9]+).*$", "\\1", comm_select$age))

# merge
d <- left_join(d, comm_select)
# total_cover
d$total_cover <- 100 - d$open_space

# read metadata
meta <- read_csv("data/metadata.csv")

# ocean basin
d <- left_join( d, select(meta, site, ocean))
# Calculate averages at site level
dsite <- d %>% 
  group_by(site, age, lat, ocean) %>% 
  summarize( ar_bryo = mean(ar_bryo, na.rm = T), col_asc = mean(col_asc, na.rm = T), sol_asc = mean(col_asc, na.rm = T), sabellids = mean(sabellids, na.rm = T),
             richness = mean(richness, na.rm = T), #shannon = mean(shannon, na.rm = T), 
             rug2 = mean(rug2, na.rm = T), total_cover = mean(total_cover, na.rm = T),
             temp_mean = mean(temp, na.rm = T), sal_mean = mean(salinity, na.rm = T)) %>% 
  mutate( logrug = log(rug2))
d <- left_join( d, select(dsite, site, lat, age, richness_mean = richness))
#### -----------------------------

# pivot the data wider to separate times
# site-level averages
dsite <- d %>% 
  select( site, age, lat, richness, shannon, rug2, temp, salinity, total_cover, ar_bryo ) %>% 
  group_by(site, age, lat ) %>% 
  summarise( temp = mean(temp, na.rm = T ), salinity = mean(salinity, na.rm = T), 
             richness = mean(richness, na.rm = T ), shannon = mean(shannon, na.rm = T), 
             rug2 = mean(rug2, na.rm = T), total_cover = mean(total_cover),
             ar_bryo = mean(ar_bryo, na.rm = T)) %>% 
  mutate( logrug = log( rug2 ) )

# get average temps and keep the rest
dsitemean <- dsite %>% 
  group_by(site) %>% 
  summarize( temp_mean = mean(temp), sal_mean = mean(salinity) )

dsite <- left_join( dsite, dsitemean )


dwide <-  dsite %>%  
  select( site, age, temp_mean, sal_mean, richness, logrug, total_cover, ar_bryo) %>% 
  pivot_wider( names_from = age, values_from = c(richness, logrug, total_cover, ar_bryo))

# dwide <- read_csv("data/data_wide.csv")  
# # rename sites
# dwide$site <- unlist( lapply( strsplit(dwide$site,"-"), function(z) z[2] ) )
# # remove NA values for relevant variables
# dpick = dwide %>% filter( ! is.na(logrug_90), ! is.na(total_cover_90)  )
# 
# # read estimates of community growth (units = percent cover per day)
# slopes <- read_csv("data/cover_rate_slopes.csv")
# 
# # merge
# dpick <- left_join( dpick, slopes )
# write to disk
# write_csv(dpick, "data/data_sem.csv")

# -----



### Structural Equation Modeling
# start with a PANEL model with just days 30 and days 60
# based on correlations, species richness is highly autocorrelated over time in the study
# but rugosity is only correlated with adjacent sampling points (30-60 days; 60-90 days)


# link to path diagrams <https://app.diagrams.net/#Hmawhal%2Fpanels-complexity%2Fmain%2Fsem%2FPanels_SEM_Panel.drawio#%7B%22pageId%22%3A%22g-VMl2kAD7bW_ptmqm9F%22%7D>
# also found in "sem/Panels_SEM_Panel.drawio" in this project



# --------------------------
# Piecewise SEM
# SEM1 - diversity influences complexity
psem1 <- psem(
  lm( estimate ~ temp_mean + sal_mean, data = d),
  lm( richness_30 ~ temp_mean + sal_mean + estimate, data = d),
  lm( logrug_90 ~ temp_mean + estimate + richness_30 + ar_bryo_30, data = d),
  lm( ar_bryo_30 ~  temp_mean + sal_mean, data = d)
)
summary(psem1)
basisSet(psem1)
dSep(psem1)

# SEM2 - complexity influences diversity
psem2 <- psem(
  lm( estimate ~ temp_mean + sal_mean , data = d),
  lm( richness_90 ~ temp_mean + sal_mean + estimate + logrug_30, data = d),
  lm( logrug_30 ~  estimate + ar_bryo_30, data = d),
  lm( ar_bryo_30 ~  temp_mean + sal_mean, data = d)
  )
summary(psem2)
basisSet(psem2)
dSep(psem2)
plot(psem2)

anova(sem1, sem2)


