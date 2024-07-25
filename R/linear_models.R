### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# Linear models, multiple regression
# site level lms and panel level lmms


# packages
library(tidyverse)
library(ggtext) # for linear and quadratic axis tick labels
# mixed modeling packages
library(readxl)
library(lme4)
library(lmerTest)
# colors
library(viridis)
# AIC comparison
library(bbmle)
# simple model plots
library(modelsummary)




# d, richness is from point counts, Shannon diversity is for genus/species level for point counts
d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
names(d) <- tolower(names(d))
d$age <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
d$panel <- gsub( "90D","90d", d$panel)
d$panel <- gsub( "60D","60d", d$panel)

# transforming rugosity measurements
d <- d %>% mutate( rugosity_raw = rugosity, rug1 = 1/rugosity_raw, rug2 = 1-rugosity_raw) %>% 
  mutate( age_factor = factor(age))
# psych::pairs.panels(d[,c("rug1","rug2")])
# psych::pairs.panels( log(d[,c("rug1","rug2")]) )

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




##
dsite_richness <- dsite %>% 
  ungroup() %>% 
  group_by( site ) %>% 
  summarize( richness = mean(richness) )
d$site_order <- factor(d$site, levels = dsite_richness$site[rev(order(dsite_richness$richness))])

# as communities develop, cover and rugosity tend to increase. 
ggplot( d, aes( x = age, y = total_cover, fill = logrug )) + 
  facet_wrap( ~ site_order) +
  geom_smooth( aes(group = site), method = "lm", se = F, col = 'darkgrey', lwd = 0.75) + 
  geom_point(alpha = 1, size = 2, pch = 21) + 
  scale_fill_viridis(option = "D", direction = 1, name = "log(rugosity)") +
  scale_x_continuous(name = "Panel age (days)", breaks = c(30,60,90), limits = c(25,95)) +
  scale_y_continuous(name = "Total % cover", breaks = c(0,50,100), limits = c(-5,105)) +
  theme_classic() +
  theme(legend.title = element_text(angle = -90, hjust = 0.5, vjust = 0.5))
ggsave("figs/cover_age_site_rug.svg", width = 4, height = 4)





# bivariate relationships

ggplot(d, aes( x = richness, y = logrug, group = site )) +
  facet_wrap(~age) +
  geom_smooth( method = "lm", se = F, col = 'gray' ) +
  geom_smooth( aes(group = 1), method = "lm", se = F ) 
# site-level richness patterns - gradients in richness and rugosity
dsite$days <- factor( dsite$age, labels = c("30 days", "60 days", "90 days") )
ggplot(dsite, aes( x = richness, y = logrug )) +
  facet_wrap( ~ days ) +
  geom_smooth( method = "lm", se = T, col = "black", lwd = 0.75 ) +
  geom_point( col = "slateblue" ) +
  ylab("log(Rugosity)") + xlab("Species richness") +
  xlim(c(0,max(dsite$richness))) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank() )
ggsave( "figs/richness_time.svg", width = 5, height = 2.25 )


  # arborescent bryozoans
ggplot(d, aes( x = ar_bryo, y = logrug, group = site )) +
  facet_wrap(~age) +
  geom_smooth( method = "lm", se = F, col = 'gray' ) +
  geom_smooth( aes(group = 1), method = "lm", se = F ) +
  geom_point(aes(col = shannon)) +
  scale_color_viridis( trans = 'log' )
# ggsave("figs/rugosity_bryos.svg", width = 4.5, height = 2.25)



my_breaks = c(1, 2, 4, 8, 16, 32)
ggplot(d, aes( x = age, y = logrug, fill = ar_bryo )) +
  facet_wrap(~site) +
  # geom_smooth( method = "glm", method.args = list(family="binomial"), 
              # se = F, col = "black", lwd = 0.75 ) +
  geom_smooth( method = "lm", se = T, col = "grey", lwd = 0.75 ) +
  geom_point( pch = 21 ) +
  scale_fill_viridis( trans = "log",
                       breaks = my_breaks, labels = my_breaks)

 

d90 <- d %>% filter(age == 90)
ggplot(d90, aes( x = richness, y = logrug )) +
  facet_wrap(~site) + 
  geom_smooth( aes(group = 1), method = "lm", se = F, col = "grey", lwd = 1 ) +
  geom_point()



###### SITE-LEVEL ANALYSIS
# simple linear model
m1 <- lm( logrug ~ richness*age, data = dsite )
summary(m1)
# allow slopes to vary
m2 <- lm( logrug ~ shannon*age, data = dsite )
summary(m2) 
m3 <- lm( logrug ~ richness*total_cover+age, data = dsite )
summary(m3) 
m4 <- lm( logrug ~ total_cover+age, data = dsite )
summary(m4) 
# psych::pairs.panels( dsite %>% select(age, total_cover, richness, logrug))

ggplot( data = dsite, aes( y = logrug, x = shannon )) + 
  facet_wrap(~age) +
  geom_point() + geom_smooth(method = "lm")
ggplot( data = dsite, aes( y = logrug, x = total_cover )) + 
  facet_wrap(~age) + 
  geom_point() + geom_smooth(method = 'lm')
ggplot( data = dsite, aes( y = logrug, x = total_cover, col = age )) + 
  geom_point() + geom_smooth(aes(group = 1), method = 'lm')




# center and scale predictors
dsite$age_fac <- factor(dsite$age, ordered = T)
dsite$age_scale <- scale(dsite$age)[,1]
dsite$cover_scale = scale(dsite$total_cover)[,1]
dsite$rich_scale = scale(dsite$richness)[,1]
dsite$arbryo_scale = scale(dsite$ar_bryo)[,1]
dsite$temp_scale = scale(dsite$temp)[,1]
dsite$sal_scale = scale(dsite$salinity)[,1]
dsite$lat_scale = scale(dsite$lat)[,1]

mage_fac      <- lm( logrug ~ age_fac, data = dsite)
summary(mage_fac) # linear effect of age is strong, but no quadratic effect
mage          <- lm( logrug ~ age_scale, data = dsite)
mcover        <- lm( logrug ~ cover_scale, data = dsite)
mrichness     <- lm( logrug ~ rich_scale, data = dsite)
marbryo       <- lm( logrug ~ arbryo_scale, data = dsite)
mtemp         <- lm( logrug ~ temp_scale, data = dsite)
msal          <- lm( logrug ~ sal_scale, data = dsite)
mlat          <- lm( logrug ~ lat_scale, data = dsite)
mfull <- lm( logrug ~ age_scale+cover_scale+rich_scale+arbryo_scale+sal_scale+lat_scale, data = dsite)
mfull2 <- lm( logrug ~ age_scale+cover_scale+rich_scale+arbryo_scale+temp_scale+sal_scale, data = dsite)
car::vif(mfull)  
AICctab( mage,mcover,mrichness,
         marbryo, 
         mtemp,msal,mlat,
         mfull,mfull2,
         nobs = length(!is.na(dsite$logrug)))
b <- list(geom_vline(xintercept = 0, color = 'orange'))
modelplot(mfull, background = b, coef_omit = "(Intercept)")
summary(mfull)$r.squared
#

# Only show the final time points
dsite90 <- dsite %>% filter(age == 90)
ggplot( data = dsite90, aes( y = logrug, x = richness )) + 
  geom_point() + geom_smooth(method = "lm")

# center and scale predictors
dsite90$cover_scale = scale(dsite90$total_cover)[,1]
dsite90$rich_scale = scale(dsite90$richness)[,1]
dsite90$arbryo_scale = scale(dsite90$ar_bryo)[,1]
dsite90$temp_scale = scale(dsite90$temp)[,1]
dsite90$sal_scale = scale(dsite90$salinity)[,1]
dsite90$lat_scale = scale(dsite90$lat)[,1]

mcover        <- lm( logrug ~- cover_scale, data = dsite90)
mrichness     <- lm( logrug ~ rich_scale, data = dsite90)
marbryo       <- lm( logrug ~ arbryo_scale, data = dsite90)
mtemp         <- lm( logrug ~ temp_scale, data = dsite90)
msal          <- lm( logrug ~ sal_scale, data = dsite90)
mlat          <- lm( logrug ~ lat_scale, data = dsite90)
mfull <- lm( logrug ~ cover_scale+rich_scale+arbryo_scale+sal_scale+lat_scale, data = dsite90)
mfull2 <- lm( logrug ~ cover_scale+rich_scale+arbryo_scale+temp_scale+sal_scale, data = dsite90)
car::vif(mfull)
AICctab( mcover,mrichness,
         marbryo,
         mtemp,msal,mlat,
         mfull,mfull2,
         nobs = length(!is.na(dsite90$logrug)))
b <- list(geom_vline(xintercept = 0, color = 'orange'))
modelplot(mfull2, background = b, coef_omit = "(Intercept)")
summary(mfull)$r.squared
summary(mrichness)
#





####  Linear Mixed Effects models of panel-level data

# model comparison
# factors of interest include age, diversity, cover of particular groups (e.g., arboresenct bryozoans, colonial tunicates)
d$age_fac <- factor(d$age, ordered = T)
d$age_scale <- scale(d$age)[,1]
d$cover_scale = scale(d$total_cover)[,1]
d$rich_scale = scale(d$richness)[,1]
d$arbryo_scale = scale(d$ar_bryo)[,1]
d$temp_scale = scale(d$temp)[,1]
d$sal_scale = scale(d$salinity)[,1]
d$lat_scale = scale(d$lat)[,1]


# consider removing Washington and low salinity sites
# They are outliers in terms of cover and environmental stress (both potentially limiting diversity)
drem <- d[ !(d$site %in% c("USA-WAS", "USA-ALD", "USA-MDA" )), ]

dmodel <- d

# simple models reflecting hypothetical relationships
mage          <- lmer( logrug ~ age_fac  + (1|site), data = dmodel, REML = F)
summary(mage)
mcover        <- lmer( logrug ~ cover_scale  + (1|site), data = dmodel, REML = F)
mrichness     <- lmer( logrug ~ rich_scale  + (1|site), data = dmodel, REML = F)
marbryo       <- lmer( logrug ~ arbryo_scale  + (1|site), data = dmodel, REML = F)
mtemp       <- lmer( logrug ~ temp_scale  + (1|site), data = dmodel, REML = F)
msal       <- lmer( logrug ~ sal_scale  + (1|site), data = dmodel, REML = F)
mlat       <- lmer( logrug ~ lat_scale  + (1|site), data = dmodel, REML = F)
mrichnessage  <- lmer( logrug ~ rich_scale*age_fac  + (1|site), data = dmodel, REML = F)
mcoverage     <- lmer( logrug ~ cover_scale*age_fac  + (1|site), data = dmodel, REML = F)
mfull  <- lmer( logrug ~ age_scale+cover_scale+rich_scale+arbryo_scale+sal_scale+lat_scale + (1|site), data = dmodel, REML = F)
mfull2  <- lmer( logrug ~ age_scale+cover_scale+rich_scale+arbryo_scale+sal_scale+temp_scale + (1|site), data = dmodel, REML = F)

AICctab( mage,mcover,mrichness,
         marbryo, mtemp, msal, mlat,
         mrichnessage,mcoverage,
         mfull, mfull2,
         nobs = length(!is.na(dmodel$logrug)))

b <- list(geom_vline(xintercept = 0, color = 'orange'))
modelplot(mfull2, background = b,  coef_omit = "Intercept|Observations")
ggsave("figs/model_results_lmer_full.svg", width = 4, height = 4)

# refit the model with REML and scale variables
mchoose <- update(mfull, REML=T)
MuMIn::r.squaredGLMM(mchoose)
# mchoose <- update(mfullf, REML=T)
plot(mchoose)
car::vif(mchoose)
ests <- data.frame( coefs = fixef(mchoose), se = sqrt(diag(vcov(mchoose))) )
ests$lower = with(ests, coefs - 2*se )
ests$upper = with(ests, coefs + 2*se )
ests <- ests[-1,]
ests$parameter <- rownames(ests)
# # paramter names
# ests$parameter <- factor(ests$parameter, 
#                          levels = (c("Total_cover","Ar_bryo","Sol_asc",
#                                     "Richness",
#                                     "Richness:Total_cover", "Day.L","Day.Q")))
# mylabs <- rev(c("total % cover","arb. bryozoans","sol. ascidians","species richness","richness x cover","day (linear)","day (quadratic)"))
# ests$parameter <- factor(ests$parameter, 
#                          levels = c("Ar_bryo","Sol_asc","Total_cover",
#                                     "Total_cover:Day60","Total_cover:Day90","Day60","Day90"),
#                          labels = c("arb. bryozoans","sol. ascidians","total cover","total cover x day 60","total cover x day 90","day 60","day 90"))
# colors
colors <- c("")
# graph
ggplot( ests, (aes( x = parameter, y = coefs))) +
  geom_errorbar( aes(ymin = lower, ymax = upper), width = 0.5 ) +
  geom_point() +
  geom_hline(yintercept = 0, lwd = 0.25) +
  ylab("Coefficient") + xlab("Parameter")+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
  ) 

ggplot( ests, (aes( x = parameter, y = coefs))) +
  geom_errorbar( aes(ymin = lower, ymax = upper), width = 0.5 ) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, lwd = 0.5, col = "orange") +
  ylab("Standardized coefficient") + xlab("")+
  theme_classic() + 
  # scale_x_discrete(labels = mylabs, limits = rev) +
  theme(axis.text.x = element_markdown()) +
  coord_flip()

ggsave("figs/final_model_lmer_full.svg", width = 3, height = 3)



## Repeat with day 90 data
d90 <- d %>% filter(age==90)
dmodel <- d90
dmodel <- left_join()

#
mcover        <- lmer( logrug ~ cover_scale  + (1|site), data = dmodel, REML = F)
mrichness     <- lmer( logrug ~ rich_scale  + (1|site), data = dmodel, REML = F)
marbryo       <- lmer( logrug ~ arbryo_scale  + (1|site), data = dmodel, REML = F)
mtemp       <- lmer( logrug ~ temp_scale  + (1|site), data = dmodel, REML = F)
msal       <- lmer( logrug ~ sal_scale  + (1|site), data = dmodel, REML = F)
mlat       <- lmer( logrug ~ lat_scale  + (1|site), data = dmodel, REML = F)

mfull  <- lmer( logrug ~ cover_scale+rich_scale+arbryo_scale+sal_scale+lat_scale + (1|site), data = dmodel, REML = F)
mfull2  <- lmer( logrug ~ cover_scale+rich_scale+arbryo_scale+sal_scale+temp_scale + (1|site), data = dmodel, REML = F)

AICctab( mcover,mrichness,
         marbryo, mtemp, msal, mlat,
         mfull, mfull2,
         nobs = length(!is.na(dmodel$logrug)))

b <- list(geom_vline(xintercept = 0, color = 'orange'))
modelplot(mfull2, background = b, coef_omit = "Observations" )
modelplot(mfull2, background = b, coef_omit = "Intercept|Observations" )
ggsave("figs/model_results_lmer_90.svg", width = 3.2, height = 3)


# 
# # plot the model results
# ## Show parameter estimates for model
# # from <https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions>
# newdat <- expand.grid(
#   age = c(30,60,90)
#   , shannon = seq(0,2.5, by = 0.1)
#   , rugosity = 0
# )
# newdat$logrug <- predict(m2,newdat,re.form=NA)
# mm <- model.matrix(terms(m2),newdat)
# ## or newdat$distance <- mm %*% fixef(fm1)
# pvar1 <- diag(mm %*% tcrossprod(vcov(m2),mm))
# tvar1 <- pvar1+VarCorr(m2)$site[1]  ## must be adapted for more complex models
# cmult <- 2 ## could use 1.96
# newdat <- data.frame(
#   newdat
#   , plo = newdat$logrug-cmult*sqrt(pvar1)
#   , phi = newdat$logrug+cmult*sqrt(pvar1)
#   , tlo = newdat$logrug-cmult*sqrt(tvar1)
#   , thi = newdat$logrug+cmult*sqrt(tvar1)
# )
# newdat$agefac <- as.factor(newdat$age)
# dshan$agefac <- as.factor(dshan$age)
# # plot confidence - fixed effects only
# g0 <- ggplot(newdat, aes(x=shannon, y=logrug, colour = age))+geom_point()
# g0 + geom_pointrange(aes(ymin = plo, ymax = phi))+
#   labs(title="CI based on fixed-effects uncertainty ONLY")
# ggplot(newdat, aes(x = shannon, y = logrug, color = agefac)) + 
#   geom_point(data = dshan, aes( x = shannon, y = logrug, color = agefac), alpha = 0.5) +
#   geom_ribbon(aes(ymin = plo, ymax = phi, group = agefac),  
#               alpha = 0.1, col = "gray") +
#   geom_path( aes(group = agefac), lwd = 1) +
#   xlab("Shannon diversity") + ylab("log(Rugosity)") +
#   theme_classic() +
#   scale_color_manual( values = c("dodgerblue","blue","black"), name = "Panel\nage")
# ggsave("figs/model_rugosity_shannon_time.svg", width = 4, height = 3)

