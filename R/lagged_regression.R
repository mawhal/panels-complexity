### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

## The potentially recursive relationship between 
# diversity and complexity might be better understood
# if conditions at the start of an experiment 

# because communities on panels were destructively sampled
# we use data points as independent conditions with
# only site-level pairing among data points

# packages
library(tidyverse)
library(readxl)

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
#### -----------------------------



# relationships to investigate
# rugosity early -> diversity later - this might be better to look at with species associated with primary substrate
# diversity early - > rugosity later
# local diversity was fairly steady, seemed to saturate early with respect to primary settlers
ggplot( d, aes( x = age, y = shannon, col = site )) + 
  geom_point() + geom_smooth(method = 'lm', se = F)
ggplot( d, aes( x = age, y = richness, col = site )) + 
  geom_point() + geom_smooth(method = 'lm', se = F)


# pivot the data wider to separate times
# site-level averages
dsite <- d %>% 
  select( site, age, lat, richness, shannon, rug2, temp, salinity, total_cover, ar_bryo ) %>% 
  group_by(site, age, lat ) %>% 
  summarise( temp = mean(temp, na.rm = T ), salinity = mean(salinity, na.rm = T), 
             richness = mean(richness, na.rm=T ), shannon = mean(shannon, na.rm=T), 
             rug2 = mean(rug2, na.rm=T), total_cover = mean(total_cover),
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

# write to disk for other analyses
write_csv( dwide, "data/data_wide.csv" )

# richness -> complexity
ggplot( dwide, aes( x = richness_30, y = logrug_30 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_30, y = logrug_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_30, y = logrug_90 )) +
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_60, y = logrug_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_60, y = logrug_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_90, y = logrug_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()


# compexity -> richness

ggplot( dwide, aes( x = logrug_30, y = richness_30 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_30, y = richness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_30, y = richness_90 )) +
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_60, y = richness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_90, y = richness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_90, y = richness_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()



# linear models
lm1 <- lm( logrug_30 ~ scale(richness_30), data = dwide )
lm2 <- lm( logrug_60 ~ scale(richness_30), data = dwide )
lm3 <- lm( logrug_90 ~ scale(richness_30), data = dwide )
lm4 <- lm( logrug_60 ~ scale(richness_60), data = dwide )
lm5 <- lm( logrug_90 ~ scale(richness_60), data = dwide )
lm6 <- lm( logrug_90 ~ scale(richness_90), data = dwide )

lm7 <- lm( richness_30 ~ scale(logrug_30), data = dwide )
lm8 <- lm( richness_60 ~ scale(logrug_30), data = dwide )
lm9 <- lm( richness_90 ~ scale(logrug_30), data = dwide )
lm10 <- lm( richness_60 ~ scale(logrug_60), data = dwide )
lm11 <- lm( richness_90 ~ scale(logrug_60), data = dwide )
lm12 <- lm( richness_90 ~ scale(logrug_90), data = dwide )


l1 <- list( lm1, lm2, lm3, lm4, lm5, lm6 )
l2 <- list( lm7, lm8, lm9, lm10, lm11, lm12 )
allmods <- c(l1,l2)
do.call( rbind, lapply(l1, function(z) summary(z)$coeff[c(2,4)] ))
ests <- data.frame( name = paste0("lm",1:12), 
                    do.call( rbind, lapply(allmods, function(z) summary(z)$coeff[c(2,4)] )),
                    do.call( rbind, lapply(allmods, function(z) confint(z)[2,] )) 
                    )
names(ests) <- c("model","est","se","lcl","ucl")

# compare AIC
library(bbmle)
aic <- AICctab( allmods, nobs = nrow(dwide) )
# model 3, 2
# model 6, 5, 4, 1
# models 7-12
aictable <- data.frame( model = rownames(as.data.frame(aic)), aic )
aictable$direction <- as.character( gl( 2, 6, labels = c("orange","blue") ) )
plot( aic$dAICc )
points( x = 1:12, y = aic$dAICc, col = aictable$direction )

aictable$model[1:2]
lm2
lm3

# pairs
# windows(5,5)
psych::pairs.panels( dwide[4:9], scale = T, ellipses = T, smooth = F, stars = F,
                     method = "pearson", 
                     hist.col = "lightcoral",
                     cex.cor = 1.75, cex = 1.5,
                     )

## Prepare the figure
ests$direction <- gl( 2, 6, labels = c("richness->complexity","complexity->richness"))
ests$comparison <- rep(c("30-30","30-60","30-90","60-60","60-90","90-90"), 2)
ests$focus_lag <- rep(c("30_lag0","60_lag1","90_lag2","60_lag0","90_lag1","90_lag0"), 2)
ests$lag <- rep(c(0,1,2,0,1,0), 2)
ests$focus <- rep(c("30 days","60 days","90 days","60 days","90 days","90 days"), 2)

ggplot(data = ests, aes(x = focus, y = est)) +
  facet_grid(lag~direction, scales = "free") +
  geom_hline( yintercept = 0, col = "orange" ) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.25 ) +
  geom_point() +
  theme_classic() +
  coord_flip()
