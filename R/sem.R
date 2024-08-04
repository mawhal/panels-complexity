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
library(piecewiseSEM)

# prepare data
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


### Structural Equation Modeling
## compare two models, each with the same number of degrees of freedom, but different directionality
# focal variables are endogenous
# some models use data from particular dates, which may complicate model comparison


## Structural equation modeling

# link to path diagrams <https://app.diagrams.net/?src=about#Hmawhal%2Fpanels-complexity%2Fmain%2Fsem%2FPanels%20SEM.drawio#%7B%22pageId%22%3A%22D_jNqRS2Lb4KAGGym6pT%22%7D>
# also found in "sem/Panels SEM.drawio"




# create SEM from individual models
# SEM1 - diversity influences complexity
sem1 <- '
  # regressions
  estimate ~ temp_mean + sal_mean
  richness_30 ~ temp_mean + sal_mean + estimate
  logrug_90 ~ temp_mean + estimate + richness_30 + ar_bryo_30
  ar_bryo_30 ~  temp_mean + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  estimate ~~ estimate
  richness_30 ~~ richness_30
  logrug_90 ~~ logrug_90
  ar_bryo_30 ~~ ar_bryo_30
  # covariances of residuals
'
fit1 <- lavaan(sem1, data = dpick)
summary(fit1, fit.measures = T, standardized = T, rsquare = T)

# SEM2 - complexity influences diversity
sem2 <- '
  # regressions
  estimate ~ temp_mean + sal_mean
  richness_90 ~ temp_mean + sal_mean + estimate + logrug_30
  logrug_30 ~ estimate + ar_bryo_30
  ar_bryo_30 ~  temp_mean + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  estimate ~~ estimate
  richness_90 ~~ richness_90
  logrug_30 ~~ logrug_30
  ar_bryo_30 ~~ ar_bryo_30
  # covariances of residuals
'
fit2 <- lavaan(sem2, data = dpick)
summary(fit2, fit.measures = T, standardized = T, rsquare = T)
tidySEM::graph_sem(fit2)
#
# model comparison
nonnest2::vuongtest( fit1, fit2, nested = FALSE )
#








# --------------------------
# Piecewise SEM
# SEM1 - diversity influences complexity
psem1 <- psem(
  lm( estimate ~ temp_mean + sal_mean, data = dpick),
  lm( richness_30 ~ temp_mean + sal_mean + estimate, data = dpick),
  lm( logrug_90 ~ temp_mean + estimate + richness_30 + ar_bryo_30, data = dpick),
  lm( ar_bryo_30 ~  temp_mean + sal_mean, data = dpick)
)
summary(psem1)
basisSet(psem1)
dSep(psem1)

# SEM2 - complexity influences diversity
psem2 <- psem(
  lm( estimate ~ temp_mean + sal_mean , data = dpick),
  lm( richness_90 ~ temp_mean + sal_mean + estimate + logrug_30, data = dpick),
  lm( logrug_30 ~  estimate + ar_bryo_30, data = dpick),
  lm( ar_bryo_30 ~  temp_mean + sal_mean, data = dpick)
  )
summary(psem2)
basisSet(psem2)
dSep(psem2)
plot(psem2)

anova(sem1, sem2)


