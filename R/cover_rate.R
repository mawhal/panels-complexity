
# packages
library(tidyverse)
library(ggtext) # for linear and quadratic axis tick labels
library(readxl)
library(viridis)



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

d$site <- unlist( lapply( strsplit(d$site,"-"), function(z) z[2] ) )


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

# add zeros to the dataset
dzeros <- bind_rows( data.frame( site = rep(unique(dsite$site), each = 6), age = 0, total_cover = 0, logrug = -8, richness = 0  ),
           d )
dzeros$site_order <- factor(dzeros$site, levels = dsite_richness$site[rev(order(dsite_richness$richness))])

# regressions to get growth rate (% area per day)
library(broom)

regressions <- dzeros %>%
  nest(data = -site) %>%
  mutate(
    fit = map(data, ~ lm(total_cover ~ age + 0, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

slopes <- regressions %>%
  unnest(tidied) %>% 
  filter(term == "age") %>% 
  select( site, estimate )

write_csv(slopes,"data/cover_rate_slopes.csv")

dzeros <- left_join(dzeros, slopes)
dzeros$site_order <- factor(dzeros$site, levels = slopes$site[rev(order(slopes$estimate))])


ggplot( dzeros, aes( x = age, y = total_cover )) + 
  facet_wrap( ~ site_order) +
  geom_smooth( aes(group = site), method = "lm", formula = y ~ x + 0, se = F, col = 'darkgrey', lwd = 0.75) + 
  geom_point(alpha = 1, size = 2, pch = 21) + 
  scale_fill_viridis(option = "D", direction = 1, name = "log(rugosity)") +
  scale_x_continuous(name = "Panel age (days)", breaks = c(0,30,60,90), limits = c(0,95)) +
  scale_y_continuous(name = "Total % cover", breaks = c(0,50,100), limits = c(-5,105)) +
  theme_classic() 


# compare to other data
dsite90 <- dsite %>% filter( age == 90 )
dsite90 <- left_join( dsite90, slopes )
ggplot(dsite90, aes(x = estimate, y = logrug)) +
  geom_point()
ggplot(dsite90, aes(x = temp, y = estimate)) +
  geom_smooth() +
  geom_point() 
