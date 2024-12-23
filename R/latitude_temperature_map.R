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
# colors
library(viridis)
# tables
library(memisc)
library(readxl)
# aesthetics
library(ggrepel)



# read the processed data. Richness is from species lists. Shannon divesity is for functional groups
# d, richness is from point counts, Shannon diversity is for genus/species level for point counts
d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
names(d) <- tolower(names(d))
d$age <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
d$panel <- gsub( "90D","90d", d$panel)
d$panel <- gsub( "60D","60d", d$panel)
d <- d %>% dplyr::mutate( rugosity_raw = rugosity, rugosity = 1-rugosity_raw) %>% 
  mutate( age_factor = factor(age))
# rename variables
d <- d %>% dplyr::rename( shannon = sh_diversity )

# meta <- read_csv("data/metadata.csv")

# community data to grab open space and arborescent bryozoans
comm_raw <- read_xlsx("data/PCover_taxgroups.xlsx")
comm_select <- comm_raw %>% dplyr::select(panel = Panel, site = Site, age = Age, ar_bryo, open_space)
comm_select$age <- as.numeric(gsub("([0-9]+).*$", "\\1", comm_select$age))


# add metadata and cover data
# meta_ocean <- meta %>% select(site, ocean)
# d <- left_join(d, meta_ocean)
d <- left_join(d, comm_select)

# species list
tlist <- read_csv("data/taxon_list.csv")
totalrich <- tlist %>% 
  group_by( site ) %>% 
  summarize( total_richness = length(unique(taxon)) )

d <- left_join(d, totalrich)




# update site names
d$site <- unlist( lapply( strsplit(d$site,"-"), function(z) z[2] ) )

# site level data
dsite <- d %>% 
  group_by(site,age,lat,salinity, temp, total_richness) %>% 
  summarise( richness = mean(richness), shannon = mean(shannon))




# bivariate relationships
# entire dataset - all sampling dates combined
d_pairs <- dsite %>% ungroup() %>% dplyr::select(lat, temp, salinity, total_richness, richness)
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
ggsave("figs/richness_latitude.svg", width = 2.5, height = 2.5)




# ranges of species richness
range(dmax$total_richness)
range(d$richness)
range(dsite$richness)



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




#### mapping 

# consider point color is temperature 
meta <- read_csv("data/metadata.csv")
# update site names
meta$site <- unlist( lapply( strsplit(meta$site,"-"), function(z) z[2] ) )
meta <- meta[!is.na(meta$Lat),]
meta <- meta %>% 
  group_by(site) %>% 
  summarize( Lat = mean(Lat), Long = mean(Long) )
meta$rowid = 1
meta$region = 1

# add temperature
dsite_means <- dsite %>%
  group_by(site) %>% 
  summarize( temp = mean(temp), salinity = mean(salinity), richness = mean(richness) )
summary(lm())
meta <- left_join( meta, dsite_means )

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
ggsave("figs/map.svg", width = 6, height = 4)
