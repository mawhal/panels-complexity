### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# this script contains visualization and analysis
# focused on gradients in measured variables across latitude and temperature


# packages
library(tidyverse)
# models
library(car)
# colors
library(viridis)
# # tables
# library(memisc)
# library(readxl)
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




# data that provides mean conditions across each site and panel age category (30,60,90 days)
dsiteage <- d %>% 
  group_by(site,age,lat,salinity, temp, total_richness) %>% 
  summarise( richness_age = mean(richness), logrug_age = mean(logrug, na.rm=T) )
richness_30 <- dsiteage %>% filter( age == 30 ) %>% 
  ungroup() %>% 
  dplyr::select( site, richness_30 = richness_age, temp_30 = temp )
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
# psych::pairs.panels(d_pairs)


# define average temperature and salinity for each site
d <- d %>% group_by( site) %>% 
  mutate( sal_mean = mean(salinity), temp_mean = mean(temp) )

dmax <- d %>% 
  group_by( site, lat, temp_mean, sal_mean, total_richness ) %>% 
  summarize( richness = max(richness) )

dsem <- read_csv("data/output/data_sem.csv")
all(dsem$site == dsite$site)
all(dsem$total_richness == dsite$total_richness)


a <- ggplot( data = dmax, aes(x = lat, y = total_richness, fill = sal_mean )) +
  geom_smooth( aes( group = 1 ), method = "lm", col = "black", lwd = 0.75 ) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  # geom_text_repel( aes(label = site), col = "slateblue" ) +
  xlab("Latitude (degrees N)") + ylab("Pooled\ntaxonomic richness") +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 
a

a1 <- ggplot( data = dmax, aes(x = lat, y = total_richness, fill = sal_mean )) +
  # geom_smooth( aes( group = 1 ), method = "lm", col = "black", lwd = 0.75 ) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   point.padding = 0.01,
                   label.padding = 0.7) +
  xlab("Latitude (degrees N)") + ylab("Pooled\ntaxonomic richness") +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 
a1

## richness and latitude at panel and site levels
# pivot longer with richness estimates as the response
dlong <- dsem %>% dplyr::select( site_name = site, lat, sal_mean, `day 30` = richness_30, `day 60` = richness_60, `day 90` = richness_90, site = total_richness )
dlong <- dlong %>% 
  pivot_longer( `day 30`:site, names_to = "richness")


a <- ggplot( data = dlong, aes(x = lat, y = value, color = richness )) +
  geom_smooth(method = 'lm', se = F, lwd = 0.5) +
  geom_smooth(data = filter(dlong, richness == "site"), 
              method = 'lm', se = T) +
  geom_point( data = filter(dlong, richness == "site"), aes( fill = sal_mean), 
              pch = 21, size = 3, show.legend = F ) +
  # geom_text_repel( data = filter(dlong, measure == "total_richness"), aes(label = site), size = 3, 
                   # point.padding = 0.01, show.legend = T)  +
  xlab("Latitude (degrees North)") + ylab("Morphospecies richness") +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() + theme( legend.position = "top") +
  guides(color = guide_legend(nrow=2, byrow = T))
a

ggplot( data = dmax, aes( x = lat, y = total_richness, col = sal_mean )) +
  geom_smooth( aes(group = 1), method = "lm", se = F, lwd = 0.75, col = "black") +
  geom_point( size = 3) +
  # geom_text_repel( aes(label = site), col = "slateblue" ) +
  ylab("Total species richness") + xlab("Latitude") +
  scale_color_viridis() +
  theme_classic() 
# ggsave("figs/richness_latitude.svg", width = 2.5, height = 2.5)
# temperature and salinity
b <- ggplot( data = dsite, aes( x = temp_30, y = richness_30, fill = sal )) +
  geom_smooth( aes(group = 1), method = "lm", se = F, lty = 2, col = "black", lwd = 0.75) +
  geom_point( aes(group = sal), pch = 21, size = 3) +
    xlim(c(10,30))+ ylim(c(0,20))+
  ylab("Panel morphospecies\nrichness (day 30)") + xlab(expression(paste("Temperature (", degree, "C) days 1-30"))) +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 
b1 <- ggplot( data = dsite, aes( x = temp_30, y = richness_30, fill = sal )) +
  geom_point( pch = 21, size = 3) +
  geom_text_repel( aes(label = site), size = 3, 
                   point.padding = 0.01,
                   label.padding = 0.5 ) +
  xlim(c(10,30))+ ylim(c(0,20))+
  ylab("Taxonomic richness\nday 30") + xlab(expression(paste("Temperature (", degree, "C) days 1-30"))) +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 
b1
 
c <- ggplot( data = dsem, aes( x = richness_30, y = logrug_90, fill = sal_mean )) +
  geom_smooth( aes(group = 1), method = "lm", col = "black", lwd = 0.75)+
  geom_point( pch = 21, size = 3) +
  xlim(c(0,20))+
  ylab("log(Rugosity)\n (day 90)") + xlab("Panel richness (day 30)") +
  scale_fill_gradientn(limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() + theme( legend.position = "top")
cowplot::plot_grid(b,c, rel_widths = c(1,1.4))
cowplot::plot_grid(a,b,c, ncol = 3, align = 'hv' )
# cowplot::plot_grid(b,c, align = "hv",nrow = 2 )
ggsave("figs/richness_panels.svg", width = 9, height = 3.5)
c1 <- ggplot( data = dsem, aes( x = richness_30, y = logrug_90, fill = sal_mean )) +
  geom_text_repel( aes(label = site), size = 3, 
                   point.padding = 0.01,
                   label.padding = 0.5 ) +
  geom_point( pch = 21, size = 3) +
  ylab("log(Rugosity)\nday 90") + xlab("Taxonomic richness\nday 30") +
  scale_fill_gradientn(limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic()
c1

d1 <- ggplot( data = dmax, aes(x = lat, y = sal_mean )) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   point.padding = 0.01,
                   label.padding = 0.7) +
  xlab("Latitude (degrees N)") + ylab("Pooled\ntaxonomic richness") +
  # scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 

e1 <- ggplot( data = dsem, aes(x = temp_mean, y = lm_middle, fill = sal_mean )) +
  # geom_smooth(se = T) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   box.padding = 0.3,
                   label.padding = 0.7) +
  xlab(expression(paste("Mean temperature (", degree, "C)"))) + ylab("Community growth rate\n(percent cover per day)") +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 
e1


# log-transformed arborescenct bryozoan
dsem$log_ar_bryo_90 <- log10( dsem$ar_bryo_90+1 )

f1 <- ggplot( data = dsem, aes(x = lm_middle, y = log_ar_bryo_90, fill = sal_mean )) +
  # geom_smooth(se = T) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   box.padding = 0.3,
                   label.padding = 0.7) +
  xlab("Community growth rate\n(percent cover per day)") + ylab("log(arborescentbryozoan\n% cover) day 90") + 
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic()
f1

h1 <- ggplot( data = dsem, aes(x = log_ar_bryo_90, y = logrug_90, fill = sal_mean )) +
  # geom_smooth(se = T) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   box.padding = 0.3,
                   label.padding = 0.7) +
  xlab("log(arborescent bryozoan\n% cover) day 90") + ylab("log(Rugosity)\nday 90") +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic()
h1

g1 <- ggplot( data = dsem, aes(x = lat, y = temp_mean )) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   point.padding = 0.01,
                   label.padding = 0.7) +
  xlab("Latitude (degrees N)") + ylab(expression(paste("Mean temperature (", degree, "C)"))) + 
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic() 

i1 <- ggplot( data = dsem, aes(x = lm_middle, y = logrug_90, fill = sal_mean )) +
  # geom_smooth(se = T) +
  geom_point(  aes( fill = sal_mean), pch = 21, size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   box.padding = 0.3,
                   label.padding = 0.7) +
  xlab("Community growth rate\n(percent cover per day)") + ylab("log(Rugosity)\nday 90") +
  scale_fill_gradientn(guide = F, limits = c(0,36), colors = c("White","magenta","darkmagenta"), name = "salinity") +
  theme_classic()
i1

j1 <- ggplot( data = dsem, aes(x = sal_mean, y = richness_30 )) +
  geom_point( size = 3 ) +
  geom_text_repel( aes(label = site), size = 3, 
                   box.padding = 0.4) +
  xlab("Mean salinity") + ylab("Taxonomic richness\nday 30") +
  theme_classic()

# windows(14,10)
panels <- cowplot::plot_grid(g1,b1,j1, a1,e1,f1,
                             h1,i1,c1,
                   ncol = 3, align = 'hv')
ggsave("figs/richness_panels_labels.svg", panels, width = 10, height = 7)


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

# ggplot( data = d, aes( x = richness, y = mfrichness, col = as.factor(age)) ) +
#   geom_point( alpha=0.5 ) + geom_smooth(se = F)
# ggplot( data = dsite, aes( x = richness, y = mfrichness) ) +
#   geom_point( ) + geom_smooth(se = F, method = 'lm')



# total richness
m <- lm( total_richness ~ lat+sal_mean+temp_mean, dsem)
#  latitutde and salinity
ma <- lm( total_richness ~ lat+sal_mean, dsem)
# temperature and salinity
mb <- lm( total_richness ~ temp_mean+sal_mean, dsem)
# # latitude alone
# mc <- lm( total_richness ~ lat, dmax)
# # temperature alone
# md <- lm( total_richness ~ temp_mean, dmax)
# # salinity alone
# me <- lm( total_richness ~ sal_mean, dmax)
# # compare models with AIC
# maic <- as.data.frame( AICctab( m, ma, mb, mc, md, me, nobs = nrow(dmax),
#          weights = T, base = T, logLik = T )  )
# maic$model <- rownames(maic)
# write_csv(maic, "tables/AIC_site_richness.csv")


vif(m)
summary(m)
vif(ma)
summary(ma)
vif(mb)
summary(mb)


# day 30
#  richness
m <- lm( richness_30 ~ lat+temp_mean+sal_mean, dsem)
#  latitutde and salinity
ma <- lm( richness_30 ~ lat+sal_mean, dsem)
# temperature and salinity
mb <- lm( richness_30 ~ temp_mean+sal_mean, dsem)
# # latitude alone
# mc <- lm( richness_30 ~ lat, dsite)
# # temperature alone
# md <- lm( richness_30 ~ temp, dsite)
# # salinity alone
# me <- lm( richness_30 ~ sal, dsite)
# # compare models with AIC
# maic2 <- as.data.frame( AICctab( m, ma, mb, mc, md, me, nobs = nrow(dmax),
#                                 weights = T, base = T, logLik = T )  )
# maic2


vif(m)
summary(m)
vif(ma)
summary(ma)
vif(mb)
summary(mb)





# ------------------------------------------


#### mapping 
mapping <- d %>% 
  group_by(site) %>% 
  summarize( Lat = mean(lat), Long = mean(long) )
mapping$rowid = 1
mapping$region = 1

# add temperature
mapping <- left_join( mapping, dplyr::select(ungroup(dsite), site, temp, sal, total_richness) )

#
library(ggthemes)

world_map = map_data("world") %>% 
  filter(! long > 20, ! lat < 0, ! lat > 70 ) %>% 
  # distinct( region ) %>% 
  rowid_to_column()

world_map %>% 
  ggplot(aes(fill = rowid, map_id = region)) +
  geom_map(map = world_map,  color="black", fill="white", size=0.25) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  coord_map("albers", lat0 = 5, lat1 = 60) +
  geom_point( data = mapping, mapping = aes(x = Long, y = Lat), col = "black", size = 3 ) +
  geom_point( data = mapping, mapping = aes(x = Long, y = Lat, col = temp), size = 2.5 ) +
  geom_text_repel(  data = mapping, aes(x = Long, y = Lat, label = site), col = "slateblue", box.padding = 0.33  ) +
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
  geom_point( data = mapping, mapping = aes(x = Long, y = Lat, fill = temp), # size = total_richness,
              pch = 21, size = 4 ) +
  geom_text_repel(  data = mapping, aes(x = Long, y = Lat, label = site), col = "slateblue", 
                    label.padding = 1 ) +
  scale_fill_viridis(name = expression(paste(degree,"C")), option = "D", limits = range(dsite$temp)) +
  theme_map() +  theme(legend.position = "top") +
  theme( panel.grid.major = element_line(colour = "grey") )
ggsave("figs/map_temp2.svg", width = 6, height = 4)
