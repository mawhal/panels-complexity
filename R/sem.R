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
library(readxl)
library(lavaan)



# load merged data
d <- read_csv("data/output/data_sem.csv")

# read metadata
meta <- read_csv("data/output/metadata.csv")
meta$site <- unlist( lapply( strsplit(meta$site,"-"), function(z) z[2] ) )



#### Add functionality for ordering by richness or arranging by ocean basin
# ocean basin
d <- left_join( d, select(meta, site, Lat, Long, ocean))


# log-transformed arborescenct bryozoan
d$log_ar_bryo_90 <- log10( d$ar_bryo_90+1 )




### Structural Equation Modeling
## compare two models, each with the same number of degrees of freedom, but different directionality
# focal variables are endogenous
# some models use data from particular dates, which may make compared models non-nested


## Structural equation modeling

# link to path diagrams <https://app.diagrams.net/?src=about#Hmawhal%2Fpanels-complexity%2Fmain%2Fsem%2FPanels%20SEM.drawio#%7B%22pageId%22%3A%22D_jNqRS2Lb4KAGGym6pT%22%7D>
# also found in "sem/Panels SEM.drawio" in this project


# create SEMs using lavaan
##### SEM1 - diversity influences complexity as hypothesized
# use panel-level total richness on day 30
sem1 <- '
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
fit1 <- lavaan(sem1, data = d)
summary(fit1, fit.measures = T, standardized = T, rsquare = T)

# include path from community growth rate to bryozoan cover
sem1a <- '
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
fit1a <- lavaan(sem1a, data = d)
summary(fit1a, fit.measures = T, standardized = T, rsquare = T)
anova(fit1, fit1a) # adding this path is supported

# arborescent bryozoans not included to test for inclusion of this variable
sem1b <- '
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
fit1b <- lavaan(sem1b, data = d)
summary(fit1b, fit.measures = T, standardized = T, rsquare = T)
anova(fit1a, fit1b) # two-degree of freedom chi-squared


##### SEM2 
# use site-level (total) richness
sem2 <- '
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
fit2 <- lavaan(sem2, data = d)
summary(fit2, fit.measures = T, standardized = T, rsquare = T)
summary(lm( total_richness ~ temp_mean+sal_mean, d))


#
# model comparison
nonnest2::vuongtest( fit1a, fit2, nested = FALSE )
#


##### SEM3 
# use site-level morphofunctional richness
sem3 <- '
  # regressions
  lm_middle ~ temp_mean
  mfrichness ~ temp_mean + sal_mean
  logrug_90 ~ lm_middle + mfrichness + log_ar_bryo_90
  log_ar_bryo_90 ~   temp_mean + lm_middle + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  lm_middle ~~ lm_middle
  mfrichness ~~ mfrichness
  logrug_90 ~~ logrug_90
  log_ar_bryo_90 ~~ log_ar_bryo_90
  # covariances of residuals
'
fit3 <- lavaan(sem3, data = d)
summary(fit3, fit.measures = T, standardized = T, rsquare = T)
summary(fit1a, fit.measures = T, standardized = T, rsquare = T)


##### SEM4 
# use site-level species richness
sem4 <- '
  # regressions
  lm_middle ~ temp_mean
  richness ~ temp_mean + sal_mean
  logrug_90 ~ lm_middle + richness + log_ar_bryo_90
  log_ar_bryo_90 ~   temp_mean + lm_middle + sal_mean
  # variances of exogenous variables
  sal_mean ~~ sal_mean
  temp_mean ~~ temp_mean
  # covariances of exogenous variables
  temp_mean ~~ sal_mean
  # residual variance for endogenous variables
  lm_middle ~~ lm_middle
  richness ~~ richness
  logrug_90 ~~ logrug_90
  log_ar_bryo_90 ~~ log_ar_bryo_90
  # covariances of residuals
'
fit4 <- lavaan(sem4, data = d)
summary(fit4, fit.measures = T, standardized = T, rsquare = T)




#
# model comparison
nonnest2::vuongtest( fit1a, fit1, nested = FALSE )
nonnest2::vuongtest( fit1a, fit2, nested = FALSE )
nonnest2::vuongtest( fit1a, fit3, nested = FALSE )
nonnest2::vuongtest( fit1a, fit4, nested = FALSE )
#


summary(fit1a, fit.measures = T, standardized = T, rsquare = T)
summary(fit3, fit.measures = T, standardized = T, rsquare = T)




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

ggplot( d, aes(x = richness_30, y = mfrichness, 
               fill = log_ar_bryo_90)) + 
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  geom_point(size = 3, pch = 21) +
  ylim( c(0,9)) + xlim(c(0,20) ) +
  ylab("Mean functional\ngroup richness") + 
  xlab("Specied richness (day 30)") +
  viridis::scale_fill_viridis( direction = -1)
ggsave("figs/mfrichness_richness30.svg", width = 5, height = 2)

d %>% select( temp_mean, lm_initial, lm_middle, glm, logrug_90 ) %>% 
  psych::pairs.panels(scale = T)

# strong correlations for initial growth rates with temperature
# 


### add latitude as a predictor variable

# include path from community growth rate to bryozoan cover
sem5 <- '
  # regressions
  temp_mean ~ Lat
  lm_middle ~ temp_mean + Lat
  richness_30 ~ temp_mean + sal_mean + Lat
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
fit5 <- lavaan(sem5, data = d)
summary(fit5, fit.measures = T, standardized = T, rsquare = T)
anova(fit1a, fit5)






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


