### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# This script prepares and compares structural equation models
# built to test whether diversity affects complexity or vice versa
# helpful weblink: <https://www.researchgate.net/post/How_to_test_the_bidirectional_relationship_in_a_SEM_model>

#### This script uses site-level data, and includes responses from different sampling dates (30,60,90 days) in the study

# packages
library(tidyverse)
library(lme4)
library(piecewiseSEM)

# # code chunk to prepare data from other scripts in the project -----
# ## this code chunk is copied from linear_models.R
# d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
# names(d) <- tolower(names(d))
# d$age <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
# d$panel <- gsub( "90D","90d", d$panel)
# d$panel <- gsub( "60D","60d", d$panel)
# 
# # transforming rugosity measurements
# d <- d %>% mutate( rugosity_raw = rugosity, rug1 = 1/rugosity_raw, rug2 = 1-rugosity_raw) %>%
#   mutate( age_factor = factor(age))
# 
# # log transform rugosity
# d$logrug <- log( d$rug2 )
# 
# # community data to grab open space and arborescent bryozoans
# comm_raw <- readxl::read_xlsx("data/PCover_taxgroups.xlsx")
# comm_select <- comm_raw %>% select(panel = Panel, site = Site, age = Age, ar_bryo, col_asc, sol_asc, sabellids, sponge, open_space)
# comm_select$age <- as.numeric(gsub("([0-9]+).*$", "\\1", comm_select$age))
# 
# # merge
# d <- left_join(d, comm_select)
# # total_cover
# d$total_cover <- 100 - d$open_space
# 
# # read metadata
# meta <- read_csv("data/metadata.csv")
# 
# #### Add functionality for ordering by richness or arranging by ocean basin
# # ocean basin
# d <- left_join( d, select(meta, site, ocean))
# # Calculate averages at site level
# dsite <- d %>%
#   group_by(site, lat, ocean) %>%
#   summarize( temp_mean = mean(temp, na.rm = T), sal_mean = mean(salinity, na.rm = T)) 
# d <- left_join( d, select(dsite, site, temp_mean, sal_mean) )
# d$site <- unlist( lapply( strsplit(d$site,"-"), function(z) z[2] ) )
# 
# ## get variables for community growth and site-level richness
# # read estimates of community growth (units = percent cover per day)
# slopes <- read_csv("data/cover_rate_slopes.csv")
# # merge
# d <- left_join( d, slopes )
# 
# # gamma diversity estimates
# # species list
# tlist <- read_csv("data/taxon_list.csv")
# totalrich <- tlist %>%
#   group_by( site ) %>%
#   summarize( total_richness = length(unique(taxon)) )
# totalrich$site <- unlist( lapply( strsplit(totalrich$site,"-"), function(z) z[2] ) )
# # merge
# d <- left_join(d, totalrich)
# # write to disk
# write_csv(d, "data/data_sem_long.csv")

# -----
# load merged data
d <- read_csv("data/data_sem_long.csv")

### Structural Equation Modeling
## compare two models, each with the same number of degrees of freedom, but different directionality
# focal variables are endogenous
# some models use data from particular dates, which may complicate model comparison


## Structural equation modeling

# link to path diagrams <https://app.diagrams.net/?src=about#Hmawhal%2Fpanels-complexity%2Fmain%2Fsem%2FPanels%20SEM.drawio#%7B%22pageId%22%3A%22D_jNqRS2Lb4KAGGym6pT%22%7D>
# also found in "sem/Panels SEM.drawio" in this project

# start with final time point (day 90)
d90 <- d %>% filter(age == 90) %>% filter( !is.na(logrug) )
d90site <- d90 %>% 
  group_by( site, temp_mean, sal_mean, estimate, total_richness ) %>% 
  summarize( richness = mean(richness), logrug = log(mean(rug2,na.rm=T)) )

# --------------------------
# Piecewise SEM
# SEM1 - diversity influences complexity
psem1 <- psem(
  lm( estimate ~ temp_mean + total_richness, data = d90site),
  lm( total_richness ~ temp_mean + sal_mean, data = d90site),
  lmer( richness ~ total_richness + (1|site), data = d90),
  # lmer( ar_bryo ~ temp_mean + estimate + (1|site), data = d90),
  lmer( logrug ~ total_richness + estimate + temp_mean + (1|site), data = d90)
  )
summary(psem1)
basisSet(psem1)
dSep(psem1)

# SEM2 - complexity influences diversity
psem2 <- psem(
  lm( estimate ~ temp_mean + total_richness, data = d90site),
  lm( total_richness ~ temp_mean + sal_mean, data = d90site),
  lmer( richness ~ total_richness + logrug + (1|site), data = d90),
  # lmer( ar_bryo ~ temp_mean + estimate + (1|site), data = d90),
  lmer( logrug ~ estimate + temp_mean + (1|site), data = d90)
)
summary(psem2)
basisSet(psem2)
dSep(psem2)
# plot(psem2)
anova(psem1, psem2)

with( d90, cor.test( richness, logrug) )
with( d90site, cor.test( richness, logrug) )
# 
d %>% select( age, richness, logrug ) %>% 
  filter(!is.na(logrug)) %>% 
  group_by( age ) %>% 
  summarize( cor( richness, logrug ) )
 
