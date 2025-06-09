### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# this script contains visualization and analysis
# focused on gradients in measured variables across latitude and temperature


# packages
library(tidyverse)
# models
library(lmerTest)
library(bbmle)
library(car)
# colors
library(viridis)
# tables
library(memisc)
library(readxl)
# aesthetics
library(ggrepel)


# 
# read the processed data from script "R/lagged_regression.R"
# Richness is from species lists. Shannon diversity is for functional groups
#
# main data sheet containing metadata, rugosity data, richness data
# d, richness is from point counts, Shannon diversity is for genus/species level for point counts
d <- read_csv("data/output/data_long.csv")
# 
# species/morphofunctional list for each panel
tlist <- read_csv("data/taxon_list.csv")
totalrich <- tlist %>%
  group_by( site ) %>%
  summarize( total_richness = length(unique(taxon)) )

d <- left_join(d, totalrich)

# update site names
d$site <- unlist( lapply( strsplit(d$site,"-"), function(z) z[2] ) )



# data that provides mean conditions across each site and panel age category (30,60,90 days)
dsiteage <- d %>% 
  group_by(site,age,lat,salinity, temp, total_richness) %>% 
  summarise( richness_age = mean(richness), logrug_age = mean(logrug, na.rm=T) )
richness_30 <- dsiteage %>% filter( age == 30 ) %>% 
  ungroup() %>% 
  dplyr::select( site, richness_30 = richness_age )
logrug_90 <- dsiteage %>% filter( age == 90 ) %>% 
  ungroup() %>% 
  dplyr::select( site, logrug_90 = logrug_age )

dsite <- d %>% 
  ungroup() %>% 
  group_by(site,lat, total_richness) %>% 
  summarise( richness = mean(richness), mfrichness = mean(mfrichness), 
               temp = mean(temp), sal = mean(salinity))


dsite <- left_join(dsite, richness_30 )
dsite <- left_join(dsite, logrug_90 )

# pivot longer
drichscale <- dsite %>% 
  dplyr::select( site, lat, temp, sal, richness_30, richness, total_richness, mfrichness )
drichscale$richness_30 = c(scale(drichscale$richness_30))
drichscale$richness = c(scale(drichscale$richness))
drichscale$total_richness = c(scale(drichscale$total_richness))
drichscale$mfrichness = c(scale(drichscale$mfrichness)) 
drichscale <- drichscale %>% 
  pivot_longer( !c(site, lat, temp, sal), names_to ="measurement" )


# bivariate relationships
# entire dataset - all sampling dates combined
ggplot( drichscale, aes( x = temp, y = value, col = measurement )) + 
  geom_smooth(se = F)+
  geom_point()

# pairs on site-level data
d_pairs <- dsite %>% ungroup() %>% dplyr::select(lat, temp, sal, total_richness, richness, mfrichness)
psych::pairs.panels(d_pairs)


# define average temperature and salinity for each site
d <- d %>% group_by( site) %>% 
  mutate( sal_mean = mean(salinity), temp_mean = mean(temp) )

dmax <- d %>% 
  group_by( site, lat, temp_mean, sal_mean, total_richness ) %>% 
  summarize( richness = max(richness) )

ggplot( data = dmax, aes(x = lat, y = richness )) +
  geom_smooth( aes(group = 1), method = "lm", se = T) +
  geom_point() +
  geom_text_repel( aes(label = site) ) +
  theme_classic()
ggplot( data = dmax, aes( x = lat, y = total_richness, col = sal_mean )) +
  geom_smooth( aes(group = 1), method = "lm", se = F, lwd = 0.75, col = "black") +
  geom_point( size = 3) +
  # geom_text_repel( aes(label = site), col = "slateblue" ) +
  ylab("Total species richness") + xlab("Latitude") +
  scale_color_viridis() +
  theme_classic() 
# ggsave("figs/richness_latitude.svg", width = 2.5, height = 2.5)
# temperature and salinity
b <- ggplot( data = dsite, aes( x = temp, y = richness_30, fill = sal )) +
  geom_point( pch = 21, size = 3) +
  # geom_text_repel( aes(label = site), col = "slateblue" ) +
  ylab("Species richness\n(day 30)") + xlab(expression(paste("Temperature (", degree, "C)"))) +
  scale_fill_gradient(guide = F) +
  theme_classic() 
ggplot( data = d, aes( x = richness, y = logrug, col = age )) +
  geom_smooth(  aes(group = age))+
  # geom_smooth( method = "lm", aes(group = age))+
  geom_point( size = 3) +
  ylab("log(Rugosity)") + xlab("Richness") +
  scale_fill_viridis() +
  theme_classic() 
dsem <- read_csv("data/output/data_sem.csv")
c <- ggplot( data = dsem, aes( x = richness_30, y = logrug_90, fill = sal_mean )) +
  # geom_smooth(method = "lm")+
  geom_point( pch = 21, size = 3) +
  ylab("log(Rugosity)\n ") + xlab("Species richness (day 30)") +
  scale_fill_gradient() +
  theme_classic() 
cowplot::plot_grid(b,c, rel_widths = c(1,1.5))
c1 <- ggplot( data = dsem, aes( x = richness_30, y = logrug_90, fill = sal_mean )) +
  # geom_smooth(method = "lm")+
  geom_text_repel( aes(label = site) ) +
  geom_point( pch = 21, size = 3) +
  ylab("log(Rugosity)\n ") + xlab("Species richness (day 30)") +
  scale_fill_gradient() +
  theme_classic()
cowplot::plot_grid(b,c, align = "hv",nrow = 2 )




# ranges of morphofunctional (largely taxonomic) richness numbers
# pooled total richness at a site
range(dmax$total_richness)
# richness on a give panel 
range(d$richness)
# average richness of sessile invertebrates on panels from each site
range(dsite$richness)



# look at residual effect of temperature after accounting for salinity and latitude
mt <-  lm( total_richness ~ lat+sal, dsite)
plot( resid(mt) ~ temp, data = dsite )
ms <-  lm( richness ~ lat+sal, dsite)
plot( resid(ms) ~ temp, data = dsite )
m30 <- lm( richness_30 ~ lat+sal, dsite)
plot( resid(m30) ~ temp, data = dsite )
mf <-  lm( mfrichness ~ lat+sal, dsite)
plot( resid(mf) ~ temp, data = dsite )
# just salinity
mt <-  lm( total_richness ~ sal, dsite)
plot( resid(mt) ~ temp, data = dsite )
ms <-  lm( richness ~ sal, dsite)
plot( resid(ms) ~ temp, data = dsite )
m30 <- lm( richness_30 ~ sal, dsite)
plot( resid(m30) ~ temp, data = dsite )
mf <-  lm( mfrichness ~ sal, dsite)
plot( resid(mf) ~ temp, data = dsite )
#
mt <-  lm( total_richness ~ salinity, d)
plot( resid(mt) ~ temp, data = d )
ms <-  lm( richness ~ salinity, d)
plot( resid(ms) ~ temp, data = d )
mf <-  lm( mfrichness ~ salinity, d)
plot( resid(mf) ~ temp, data = d )


# test for quadratic term
lm1 <- lm( total_richness ~ temp + sal, dsite)
lm2 <- lm( total_richness ~ temp + I(temp^2) + sal, dsite)
AICctab(lm1,lm2, nobs = nrow(dsite))
summary(lm2)

ggplot( data = d, aes( x = richness, y = mfrichness, col = as.factor(age)) ) +
  geom_point( alpha=0.5 ) + geom_smooth(se = F)
ggplot( data = dsite, aes( x = richness, y = mfrichness) ) +
  geom_point( ) + geom_smooth(se = F, method = 'lm')



# total richness
m <- lm( total_richness ~ lat+temp_mean+sal_mean, dmax)
#  latitutde and salinity
ma <- lm( total_richness ~ lat+sal_mean, dmax)
# temperature and salinity
mb <- lm( total_richness ~ temp_mean+sal_mean, dmax)
# latitude alone
mc <- lm( total_richness ~ lat, dmax)
# temperature alone
md <- lm( total_richness ~ temp_mean, dmax)
# salinity alone
me <- lm( total_richness ~ sal_mean, dmax)
# compare models with AIC
maic <- as.data.frame( AICctab( m, ma, mb, mc, md, me, nobs = nrow(dmax),
         weights = T, base = T, logLik = T )  )
maic$model <- rownames(maic)
write_csv(maic, "tables/AIC_site_richness.csv")

summary(m)
vif(lm( total_richness ~ lat+temp_mean+sal_mean, data = dmax))
summary(ma)
vif(m)


# total richness
m <- lm( richness_30 ~ lat+temp+sal, dsite)
#  latitutde and salinity
ma <- lm( richness_30 ~ lat+sal, dsite)
# temperature and salinity
mb <- lm( richness_30 ~ temp+sal, dsite)
# latitude alone
mc <- lm( richness_30 ~ lat, dsite)
# temperature alone
md <- lm( richness_30 ~ temp, dsite)
# salinity alone
me <- lm( richness_30 ~ sal, dsite)
# compare models with AIC
maic2 <- as.data.frame( AICctab( m, ma, mb, mc, md, me, nobs = nrow(dmax),
                                weights = T, base = T, logLik = T )  )
maic2



#### mapping 

# consider point color is temperature 
meta <- read_csv("data/output/metadata.csv")
# update site names
meta$site <- unlist( lapply( strsplit(meta$site,"-"), function(z) z[2] ) )
meta <- meta[!is.na(meta$Lat),]
meta <- meta %>% 
  group_by(site) %>% 
  summarize( Lat = mean(Lat), Long = mean(Long) )
meta$rowid = 1
meta$region = 1

# add temperature
meta <- left_join( meta, dplyr::select(ungroup(dsite), site, temp, sal, total_richness) )

#
library(ggthemes)

world_map = map_data("world") %>% 
  filter(! long > 10, ! lat < 0, ! lat > 70 ) %>% 
  # distinct( region ) %>% 
  rowid_to_column()

world_map %>% 
  ggplot(aes(fill = rowid, map_id = region)) +
  geom_map(map = world_map,  color="black", fill="white", size=0.25) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("albers", lat0 = 5, lat1 = 60) +
  geom_point( data = meta, mapping = aes(x = Long, y = Lat), col = "black", size = 3 ) +
  geom_point( data = meta, mapping = aes(x = Long, y = Lat, col = temp), size = 2.5 ) +
  geom_text_repel(  data = meta, aes(x = Long, y = Lat, label = site), col = "slateblue", box.padding = 0.33  ) +
  scale_color_viridis(name = expression(paste(degree,"C")), option = "C") +
  theme_map() +  theme(legend.position = "top") +
  guides( fill = "none", labels = "temperature" )
ggsave("figs/map_temp.svg", width = 6, height = 4)

world_map %>% 
  ggplot(aes(fill = rowid, map_id = region)) +
  geom_map(map = world_map,  color="black", fill="white", size=0.25) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("albers", lat0 = 5, lat1 = 60) +
  # geom_point( data = meta, mapping = aes(x = Long, y = Lat), col = "black", size = 3 ) +
  geom_point( data = meta, mapping = aes(x = Long, y = Lat, fill = temp, size = total_richness),
              pch = 21 ) +
  geom_text_repel(  data = meta, aes(x = Long, y = Lat, label = site), col = "slateblue", box.padding = 0.33  ) +
  scale_fill_viridis(name = expression(paste(degree,"C")), option = "D", limits = range(dsite$temp)) +
  theme_map() +  theme(legend.position = "top") +
  theme( panel.grid.major = element_line(colour = "grey") )
ggsave("figs/map_rich.svg", width = 6, height = 4)
