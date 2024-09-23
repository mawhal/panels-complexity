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
library(lavaan)
# library(piecewiseSEM)

# code chunk below prepare data from other scripts in the project -----
# read wide data from script "R/lagged_regression.R"
dwide <- read_csv("data/data_wide.csv")
# rename sites
dwide$site <- unlist( lapply( strsplit(dwide$site,"-"), function(z) z[2] ) )
# remove NA values for relevant variables
dpick = dwide %>% filter( ! is.na(logrug_90), ! is.na(total_cover_90)  )

# read estimates of community growth (units = percent cover per day)
slopes <- read_csv("data/cover_rate_slopes.csv")

# merge
dpick <- left_join( dpick, slopes )
# write to disk
write_csv(dpick, "data/data_sem.csv")

# -----

# load merged data
d <- read_csv("data/data_sem.csv")

# species list
tlist <- read_csv("data/taxon_list.csv")
totalrich <- tlist %>% 
  group_by( site ) %>% 
  summarize( total_richness = length(unique(taxon)) )
totalrich$site <- unlist( lapply( strsplit(totalrich$site,"-"), function(z) z[2] ) )

d <- left_join(d, totalrich)



### Structural Equation Modeling
## compare two models, each with the same number of degrees of freedom, but different directionality
# focal variables are endogenous
# some models use data from particular dates, which may complicate model comparison


## Structural equation modeling

# link to path diagrams <https://app.diagrams.net/?src=about#Hmawhal%2Fpanels-complexity%2Fmain%2Fsem%2FPanels%20SEM.drawio#%7B%22pageId%22%3A%22D_jNqRS2Lb4KAGGym6pT%22%7D>
# also found in "sem/Panels SEM.drawio" in this project

# log-transformed arborescenct bryozoan
d$log_ar_bryo_90 <- log10( d$ar_bryo_90+1 )

# create SEM using lavaan
##### SEM1 - diversity influences complexity as hypothesized
sem1 <- '
  # regressions
  lm_middle ~ temp_mean
  richness_30 ~ temp_mean + sal_mean
  logrug_90 ~ temp_mean + lm_middle + richness_30
  # ar_bryo_30 ~  temp_mean + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  lm_middle ~~ lm_middle
  richness_30 ~~ richness_30
  logrug_90 ~~ logrug_90
  # ar_bryo_30 ~~ ar_bryo_30
  # covariances of residuals
'
fit1 <- lavaan(sem1, data = d)
summary(fit1, fit.measures = T, standardized = T, rsquare = T)


##### include arborescent bryozoans
sem2a <- '
  # regressions
  lm_middle ~ temp_mean
  richness_30 ~ temp_mean + sal_mean
  logrug_90 ~ lm_middle + richness_30 + log_ar_bryo_90
  log_ar_bryo_90 ~   temp_mean  + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  lm_middle ~~ lm_middle
  richness_30 ~~ richness_30
  logrug_90 ~~ logrug_90
  log_ar_bryo_90 ~~ log_ar_bryo_90
  # covariances of residuals
'
fit2a <- lavaan(sem2a, data = d)
summary(fit2a, fit.measures = T, standardized = T, rsquare = T)
sem2 <- '
  # regressions
  lm_middle ~ temp_mean
  richness_30 ~ temp_mean + sal_mean
  logrug_90 ~ lm_middle + richness_30 + log_ar_bryo_90
  log_ar_bryo_90 ~   temp_mean + lm_middle + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  lm_middle ~~ lm_middle
  richness_30 ~~ richness_30
  logrug_90 ~~ logrug_90
  log_ar_bryo_90 ~~ log_ar_bryo_90
  # covariances of residuals
'
fit2a <- lavaan(sem2a, data = d)
summary(fit2a, fit.measures = T, standardized = T, rsquare = T)
anova(fit2, fit2a)

sem3 <- '
  # regressions
  lm_middle ~ temp_mean
  total_richness ~ temp_mean + sal_mean
  logrug_90 ~ lm_middle + total_richness + log_ar_bryo_90
  log_ar_bryo_90 ~   temp_mean + lm_middle + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  lm_middle ~~ lm_middle
  total_richness ~~ total_richness
  logrug_90 ~~ logrug_90
  log_ar_bryo_90 ~~ log_ar_bryo_90
  # covariances of residuals
'
fit3 <- lavaan(sem3, data = d)
summary(fit3, fit.measures = T, standardized = T, rsquare = T)


#
# model comparison
anova(fit2,fit2)
nonnest2::vuongtest( fit2, fit3, nested = FALSE )
#


# pairwise
ggplot( d, aes(x = temp_mean, y = glm) ) + geom_point()
ggplot( d, aes(x = temp_mean, y = lm_middle) ) + geom_smooth() + geom_point()
ggplot( d, aes(x = temp_mean, y = lm_initial) ) + geom_point()
ggplot( d, aes(x = temp_mean, y = lm_all) ) + geom_point()

ggplot( d, aes(x = glm, y = logrug_90) ) + geom_point()
ggplot( d, aes(x = lm_initial, y = logrug_90) ) + geom_point()
ggplot( d, aes(x = lm_middle, y = logrug_90) ) + geom_point()
ggplot( d, aes(x = lm_all, y = logrug_90) ) + geom_point()

ggplot( d, aes(x = total_cover_90, y = logrug_90) ) + geom_point()
ggplot( d, aes(x = total_cover_60, y = logrug_90) ) + geom_point()
ggplot( d, aes(x = total_cover_30, y = logrug_90) ) + geom_point()

d %>% select( temp_mean, lm_initial, lm_middle, glm, logrug_90 ) %>% 
  psych::pairs.panels(scale = T)

# strong correlations for initial growth rates with temperature
# 







# --------------------------
# Piecewise SEM
# SEM1 - diversity influences complexity
library(piecewiseSEM)
psem1 <- psem(
  lm( lm_middle ~ temp_mean, data = d),
  lm( richness_30 ~ temp_mean + sal_mean, data = d),
  lm( logrug_90 ~ lm_middle + richness_30 + log_ar_bryo_90, data = d),
  lm( log_ar_bryo_90 ~  temp_mean + lm_middle + sal_mean, data = d)
)
summary(psem1)
basisSet(psem1)
dSep(psem1)

# SEM2 - complexity influences diversity
psem2 <- psem(
  lm( lm_middle ~ temp_mean, data = d),
  lm( total_richness ~ temp_mean + sal_mean, data = d),
  lm( logrug_90 ~ lm_middle + total_richness + log_ar_bryo_90, data = d),
  lm( log_ar_bryo_90 ~  temp_mean + lm_middle + sal_mean, data = d)
)
summary(psem2)
basisSet(psem2)
dSep(psem2)
plot(psem2)

anova(sem1, sem2)


