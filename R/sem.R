### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# This script prepares and compares structural equation models
# built to test whether diversity affects complexity or vice versa
# helpful weblink: <https://www.researchgate.net/post/How_to_test_the_bidirectional_relationship_in_a_SEM_model>


# packages
library(tidyverse)
library(readxl)
library(lme4)


# prepare data

## the code chunk below is copied from linear_models.R
#### -----------------------------
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

# rename variables
d <- d %>% rename( shannon = sh_diversity )

# ## OUTLIER SITES in terms of environment and community growth
# # some outlier sites
# doutlier <- d[ (d$site %in% c("USA-WAS", "USA-ALD", "USA-MDA" )), ]
# # low salinity sites are Alabama and Maryland, while Washington had low cover throughout the study



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



#### Add functionality for ordering by richness or arranging by ocean basin
# ocean basin
d <- left_join( d, select(meta, site, ocean))
# Calculate averages at site level
dsite <- d %>% 
  group_by(site, age, temp, salinity, lat, ocean) %>% 
  summarize( ar_bryo = mean(ar_bryo, na.rm=T), col_asc = mean(col_asc, na.rm=T), sol_asc = mean(col_asc, na.rm=T), sabellids = mean(sabellids, na.rm=T),
             richness = mean(richness, na.rm=T), shannon = mean(shannon, na.rm=T), 
             rug2 = mean(rug2, na.rm=T), total_cover = mean(total_cover, na.rm=T)) %>% 
  mutate( logrug = log(rug2))
d <- left_join( d, select(dsite, site, age, richness_mean = richness))



# last time step
d90 <- d %>% filter(age == 90)

# read wide data
dwide <- read_csv("data/data_wide.csv")
dwide$site <- unlist( lapply( strsplit(dwide$site,"-"), function(z) z[2] ) )

dpick = dwide %>% filter( ! is.na(logrug_90), ! is.na(total_cover_90)  )

slopes <- read_csv("data/cover_rate_slopes.csv")

dpick <- left_join( dpick, slopes )


### Structural Equation Modeling
## compare two models, each with the same number of degrees of freedom, but different directionality
# focal variables are endogenous
# some models use data from particular dates, which may complicate model comparison


## Piecewise SEM with mixed effects models
library(piecewiseSEM)

#### site-level data

# create SEM from individual models
sem1 <- psem(
  lm( estimate ~ temp_mean + sal_mean, data = dpick),
  lm( richness_30 ~ temp_mean + sal_mean + estimate, data = dpick),
  lm( logrug_90 ~ temp_mean + estimate + richness_30, data = dpick)
)
summary(sem1)
basisSet(sem1)
dSep(sem1)

sem2 <- psem(
  lm( estimate ~ temp_mean + sal_mean , data = dpick),
  lm( richness_90 ~ temp_mean + sal_mean + estimate + logrug_30, data = dpick),
  lm( logrug_30 ~  estimate, data = dpick)
  )
summary(sem2)
basisSet(sem2)
dSep(sem2)

anova(sem1, sem2)



# ####
# # create SEM from individual models
# sem1 <- psem(
#   lm( total_cover_30 ~ temp_mean + sal_mean, data = dpick),
#   lm( richness_30 ~ temp_mean + sal_mean + total_cover_30, data = dpick),
#   lm( logrug_90 ~ temp_mean + total_cover_30 + richness_30, data = dpick)
# )
# summary(sem1)
# basisSet(sem1)
# dSep(sem1)
# 
# sem2 <- psem(
#   lm( total_cover_90 ~ temp_mean + sal_mean + logrug_30, data = dpick),
#   lm( richness_90 ~ temp_mean + sal_mean + total_cover_90 + logrug_30, data = dpick),
#   lm( logrug_30 ~ temp_mean, data = dpick)
#   )
# summary(sem2)
# basisSet(sem2)
# dSep(sem2)
# 
# anova(sem1, sem2)
# ####





### panel-level data using lmer
# create SEM from individual models
sem1 <- psem(
  lmer( total_cover ~ temp + (1|site/age), data = dpick),
  lmer( richness ~ temp + logrug + total_cover + (1|site/age), data = dpick),
  lmer( logrug ~ total_cover + (1|site/age), data = dpick),
  data = dpick
  )

summary(sem1)
basisSet(sem1)
dSep(sem1)



sem2 <- psem(
  lmer( total_cover ~ temp + richness + (1|site/age), data = dpick),
  lmer( richness ~ temp + (1|site/age), data = dpick),
  lmer( logrug ~ total_cover + richness + (1|site/age), data = dpick),
  data = dpick
)

summary(sem2)

anova(sem1,sem2)
