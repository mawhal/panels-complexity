### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# this script runs models of sessile invertebrate species richness at the site level and panel level 


# packages
library(tidyverse)
library(ggh4x)
library(lme4)
library(lmerTest)

# ------------
# prep data

# 
# read the processed data from script "R/lagged_regression.R"
# Richness is from species lists. Shannon diversity is for functional groups
#
# main data sheet containing metadata, rugosity data, richness data
# d, richness is from point counts, Shannon diversity is for genus/species level for point counts
d <- read_csv("data/output/data_long.csv")
# 




# data that provides mean conditions across each site and panel age category (30,60,90 days)
dsiteage <- d %>% 
  group_by(site,age,lat,salinity, temp, total_richness) %>% 
  summarise( richness_age = mean(richness), logrug_age = mean(logrug, na.rm=T) )
richness_30 <- dsiteage %>% filter( age == 30 ) %>% 
  ungroup() %>% 
  dplyr::select( site, richness_30 = richness_age, temp_30 = temp )
richness_60 <- dsiteage %>% filter( age == 60 ) %>% 
  ungroup() %>% 
  dplyr::select( site, richness_60 = richness_age, temp_60 = temp )
richness_90 <- dsiteage %>% filter( age == 90 ) %>% 
  ungroup() %>% 
  dplyr::select( site, richness_90 = richness_age, temp_90 = temp )
logrug_90 <- dsiteage %>% filter( age == 90 ) %>% 
  ungroup() %>% 
  dplyr::select( site, logrug_90 = logrug_age )

dsite <- d %>% 
  ungroup() %>% 
  group_by(site,lat, total_richness) %>% 
  summarise( richness = mean(richness), mfrichness = mean(mfrichness), 
             temp = mean(temp), sal = mean(salinity))


dsite <- left_join(dsite, richness_30 )
dsite <- left_join(dsite, richness_60 )
dsite <- left_join(dsite, richness_90 )
dsite <- left_join(dsite, logrug_90 )


# ------------
# data and distributions
hist(d$richness)
hist(dsiteage$richness_age)


# look at distributions of panel- and site-level data

## PANEL-LEVEL
# normal
# Plot histogram with normal distribution fit
ggplot(d, aes(x = richness)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, color = "white", fill = "skyblue") + # Use density for y-axis
  stat_function(fun = dnorm, color = "red", size = 1, # Overlay the normal density curve
                args = list(mean = mean(d$richness), sd = sd(d$richness))) +
  labs(title = "Normal Distribution Fit to Frequency Data",
       y = "Density",
       x = "Value")

# poisson
set.seed(42)

ggplot(d, aes(richness)) +
  geom_bar(aes(fill = ocean, y = after_stat(prop)),
           alpha = 0.5, width = 1, position = "identity") +
  stat_theodensity(aes(colour = ocean), distri = "pois") +
  stat_function(fun = dnorm, color = "red", size = 1, # Overlay the normal density curve
                args = list(mean = mean(d$richness), sd = sd(d$richness))) 
  


# linear mixed effects models
d$age_center = d$age/30 - 2
mnorm1 <- lmer( richness ~ lat + (1|site), data = d )
plot(mnorm1)
summary(mnorm1)
mnorm2 <- lmer( richness ~ lat + temp + salinity + (1|site), data = d )
plot(mnorm2)
summary(mnorm2)
mnorm3 <- lmer( richness ~ temp + salinity + (1|site), data = d )
plot(mnorm3)
summary(mnorm3)

mnorm4 <- lmer( richness ~ age_center +  (1|site), data = d )
summary(mnorm4)
mnorm5 <- lmer( richness ~ age_center + age_center:site + (1|site), data = d )
summary(mnorm5)
mnorm5.5 <- lmer( richness ~ age_center:site + (1|site), data = d )
summary(mnorm5.5)
fixef(mnorm5.5) * 2
mnorm6 <- lmer( richness ~ age_center + (age_center|site), data = d )
summary(mnorm6)



mpois1 <- glmer( richness ~ lat + (1|site), family = "poisson", data = d )
plot(mpois1)
summary(mpois1)
mpois2 <- glmer( richness ~ lat + temp + salinity + (1|site), family = "poisson", data = d )
plot(mpois2)
summary(mpois2)
mpois3 <- glmer( richness ~ temp + salinity + (1|site), family = "poisson", data = d )
plot(mpois3)
summary(mpois3)


# PANEL RICHNESS WITH SITE RANDOM EFFECT
mpois1 <- glmer( richness ~ age_center + (1|site), family = "poisson", data = d )
mpois2 <- glmer( richness ~ age_center + (age_center|site), family = "poisson", data = d )
anova(mpois1, mpois2)
plot(mpois2)
summary(mpois2)

newdat <- data.frame( age_center = c(-1,1) )


mpois3 <- glmer( richness ~ age_center:site + (1|site), family = "poisson", data = d )
summary(mpois3)
# predictions are roughly the same for each modeling approach (normal vs Poisson models)

ggplot( data = d, aes(x = age, y = richness, group = 1)) +
  geom_smooth( method = 'glm', se = T, method.args = list(family = "poisson")) +
  geom_point()
ggplot( data = filter(d, site == "MAD"), aes(x = age, y = richness, group = 1)) +
  geom_smooth( method = 'glm', se = T, method.args = list(family = "poisson")) +
  geom_smooth( method = 'lm', se = T) +
  geom_point()

## ------------
# show trend for each site (independently fitted trendlines)
ggplot( data = d, aes( y = richness, x = age)) +
  facet_wrap( ~ site ) +
  geom_smooth( method = 'lm', se = T ) +
  geom_point()
# sites treated as independent
mfixed1 <- lm( richness ~ age_center*site, data = d )
summary(mfixed1)
## ------------



# Mean PANEL richness (Site level analysis)
snorm1 <- lm( richness_age ~ age, data = dsiteage )
summary(snorm1)

ggplot( data = dsiteage, aes(x = age, y = richness_age, group = site)) +
  geom_smooth( method = 'lm', se = F) +
  # geom_smooth( method = 'glm', se = F, method.args = list(family = "poisson")) +
  geom_point()
# only Madeira
ggplot( data = filter(dsiteage, site == "MAD"), aes(x = age, y = richness_age, group = site)) +
  geom_smooth( method = 'glm', se = F, method.args = list(family = "poisson")) +
  geom_point()

ggplot( data = dsiteage, aes(x = age, y = richness_age, group = 1)) +
  geom_smooth( method = 'glm', se = T, method.args = list(family = "poisson")) +
  geom_point()
ggplot( data = filter(dsiteage, site == "MAD"), aes(x = age, y = richness_age, group = 1)) +
  geom_smooth( method = 'glm', se = T, method.args = list(family = "poisson")) +
  geom_point()

