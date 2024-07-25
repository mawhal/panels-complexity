### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# this script contains exploratory visualization and analysis
# of data from the project
# using processed data from panels and environmental measures
# including community indices, temperature, and salinity


# packages
library(tidyverse)
# mixed modeling packages
library(lme4)
library(lmerTest)
# colors
library(viridis)


# d, richness is from point counts, Shannon diversity is for genus/species level for point counts
d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
names(d) <- tolower(names(d))
d$age <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
d$panel <- gsub( "90D","90d", d$panel)
d$panel <- gsub( "60D","60d", d$panel)
# one panel appears to be missing: "2021_USA-WAS_30d_06"


# rugosity 
# boxplot(d$rugosity, ylim = c(0,1) )
# boxplot(1-d$rugosity, ylim = c(0,1) )
d <- d %>% mutate( rugosity_raw = rugosity, rugosity = 1-rugosity_raw) %>% 
  mutate( age_factor = factor(age))

# rename variables
d <- d %>% rename( shannon = sh_diversity )

# variables over time
ggplot( d, aes(x = age, y = temp, group = site)) + 
  geom_path() + geom_point()
ggplot( d, aes(x = age, y = salinity, group = site)) + 
  geom_path() + geom_point()
ggplot( d, aes(x = age, y = shannon, group = site)) + 
  facet_wrap(~site) + geom_point()
ggplot( d, aes(x = age, y = rugosity, group = site)) + 
  facet_wrap(~site) + geom_point()

# bivariate relationships
ggplot( d, aes( x = salinity, y = temp )) + geom_point()

ggplot(d, aes( x = richness, y = rugosity, col = age )) +
  facet_wrap(~site) +
  geom_point()
ggplot(d, aes( x = richness, y = rugosity, group = site )) +
  facet_wrap(~age) +
  geom_smooth( method = "lm", se = F, col = 'gray' ) +
  geom_smooth( aes(group = 1), method = "lm", se = F ) 
  # geom_point()
ggplot(d, aes( x = shannon, y = log(rugosity), col = site )) +
  facet_wrap(~age) +
  geom_point()

# which sites and times had zero shannon diversity
d %>% filter( shannon == 0 ) %>% 
  select(site, age, salinity, rugosity, richness)

#






# USA-MDA has low rugosity and low richness
# - this site was quickly dominated by encrusting bryozoans and is a low salinity site
# USA-ALD has low salinity like USA-MDA
#
# USA-WAS had low cover throughout the study and was dominated by a single 
# 

ggplot( d, aes(x = age, y = richness, group = site, color = rugosity)) + 
  facet_wrap(~site) + geom_point() + geom_smooth(method = "lm") 


# analysis 
ggplot(d, aes( x = richness, y = rugosity, fill = age_factor )) +
  facet_wrap(~site, scales = "free") +
  geom_smooth( aes(group = 1), method = "glm", method.args = list(family="binomial"), 
              se = F, col = "black", lwd = 0.75 ) +
  geom_smooth( aes(group = 1), method = "lm", se = F, col = "grey", lwd = 0.75 ) +
  geom_point( pch = 21 ) +
  scale_fill_manual( values = c("white","grey","black"))
my_breaks = c(1, 2, 4, 8, 16, 32)
ggplot(d, aes( x = age, y = rugosity, col = richness )) +
  facet_wrap(~site, scales = "free") +
  geom_smooth( method = "glm", method.args = list(family="binomial"), 
              se = F, col = "black", lwd = 0.75 ) +
  geom_smooth( method = "lm", se = F, col = "grey", lwd = 0.75 ) +
  geom_point() +
  scale_color_viridis( trans = "log",
                       breaks = my_breaks, labels = my_breaks)
  


# model of rugosity as a function of species richness and time
# simple linear model
m1 <- lmer( rugosity ~ richness*age + (1|site), data = d )
summary(m1)
# allow slopes to vary
m2 <- lmer( rugosity ~ richness*age + (richness|site), data = d )
summary(m2)
# allow slopes to vary
m2.5 <- lmer( rugosity ~ shannon*age + (shannon|site), data = d )
summary(m2.5) # effects depend on age, get stronger with age, but sites with higher diversity may have a weaker or negative relationship



#### just use the final time point
d90 <- d %>% filter(age == 90)
ggplot(d90, aes( x = (richness), y = rugosity, group = site )) +
  geom_smooth( aes(group = 1)) + 
  # geom_smooth( aes(group = 1), method = "lm") +
  geom_point()
  
ggplot(d90, aes( x = (shannon), y = rugosity, group = site )) +
  geom_smooth( aes(group = 1)) + 
  # geom_smooth( aes(group = 1), method = "lm") +
  geom_point()
  
ggplot(d90, aes( x = richness, y = rugosity, group = site )) +
  geom_smooth( method = "lm", se = F, col = 'gray' ) +
  geom_smooth( aes(group = 1), method = "lm", se = F ) +
  geom_point()
ggplot(d90, aes( x = shannon, y = rugosity, group = site )) +
  geom_smooth( method = "lm", se = F, col = 'gray' ) +
  geom_smooth( aes(group = 1), method = "lm", se = F ) +
  geom_point()


m3 <- lmer( rugosity ~ richness + (1|site), data = d90 )
summary(m3)
m4 <- lmer( rugosity ~ shannon + (1|site), data = d90 )
summary(m4)


# diversity of functional groups
ggplot(d90, aes( x = shannon, y = rugosity, group = site )) +
  geom_smooth( method = "lm", se = F, col = 'gray' ) +
  geom_smooth( aes(group = 1), method = "lm", se = F ) +
  geom_point()
ggplot(d90, aes( x = shannon, y = rugosity, col = richness )) +
  geom_smooth( aes(group = 1),  se = T ) +
  geom_point() +
  scale_color_viridis( trans = "log",
                       breaks = my_breaks, labels = my_breaks)
ggplot(d90, aes( x = richness, y = shannon )) +
  geom_smooth() +
  geom_point()

m4 <- lmer( rugosity ~ shannon + (1|site), data = d )
summary(m4)


##### Notes for MS #########
# - maybe not surprisingly, the relationship between diversity and complexity is complex
# - at several sites, a relationship between richness and rugosity emerges over time
#   presumably as communities go through succession and individuals/colonies grow
# - a few sites go against this trend, especially USA-IRL, 
#   where rugosity increases over time as species richness decreased
# - time is an important element. 
#   - As communities develop, both biodiversity and complexity/rugosity change in different ways



# Diversity-complexity as a two-way street
psych::pairs.panels( d %>% select(rugosity, sh_diversity, richness) )


q10 <- seq(0.1, 0.9, by = 0.1)
ggplot( data = d90, aes( x = rugosity, y = richness) ) +
  geom_point() + geom_quantile(quantiles = q10) +
  ylab("Taxon richness on day 90")
ggplot( data = d90, aes( x = rugosity, y = sh_diversity) ) +
  geom_vline(xintercept = 0.05) +
  geom_point() + geom_smooth( aes(group=site), method = "lm", se = F) +
  ylab("Taxon richness on day 90")

# summarize data at site level (not much going on within sites here)
d90site <- d90 %>% 
  group_by(site, lat) %>% 
  summarize( rugosity = mean(rugosity), richness = mean(richness), sh_diversity = mean(sh_diversity))
ggplot( data = d90site, aes( x = rugosity, y = sh_diversity) ) +
  geom_vline(xintercept = 0.05) +
  # geom_smooth( se = T) +
  geom_smooth( method = 'lm', se = T) +
  geom_smooth( data = d90site %>% filter(rugosity > 0.05), method = 'lm', se = T) +
  geom_point( alpha = 1) +
  ylab("Shannon diversity on day 90") + xlab("Rugosity on day 90")
  
