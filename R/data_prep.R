### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###


# scripts used to prepare data for various uses


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
comm_raw$Panel[comm_raw$Panel == "2021_USA-LIS_90D_17"] <- "2021_USA-LIS_90d_17"
# pull out the relevant columns
comm_select <- comm_raw %>% select(panel = Panel, site = Site, age = Age, ar_bryo, col_asc, sol_asc, sabellids, sponge, open_space)
# reformat panel sampling dates (panel ages of 30, 60, and 90 days)
comm_select$age <- as.numeric(gsub("([0-9]+).*$", "\\1", comm_select$age))
# merge
d <- left_join(d, comm_select)
# community data to calcuclate morphospecies richness
# extract the community data set 
comm <- comm_raw %>% dplyr::select( algae:sponge) %>% dplyr::select(-open_space)
comm_meta <- comm_raw[1:3]
names(comm_meta) <- tolower(names(comm_meta))
### compare morphofunctional richness to that of species richness
# convert cover data to presence/absence
comm_pa <- ifelse(comm == 0, 0, 1)
comm_meta$mfrichness <- rowSums(comm_pa)
# richness data
comm_meta$age <- as.numeric(gsub("([0-9]+).*$", "\\1", comm_meta$age))
d <- left_join(d, comm_meta)

# total_cover
d$total_cover <- 100 - d$open_space

# read metadata - see script "R/cover_data_metadata.R" for source code
meta <- read_csv("data/output/metadata.csv")


## changes made to csv generated from originally received "metadata.xlsx" file
# # add Ocean Basin
# meta$ocean <- ifelse( meta$Long < -98, "Pacific", "Atlantic") # this works even for Panama and Texas sites
# 
# # rename Dauphin Island site based on the Pier site, which seems to be the one listed in the data summaries
# meta$site[ meta$site == "USA-ALD-PIER"] <- "USA-ALD"
# # salinity and temperature data
# salinity <- meta %>%
#   pivot_longer( `Salinity - 30 days`:`Salinity - 90 days`,
#                 names_to = "name", values_to = "salinity") %>%
#   mutate( age = gsub( "\\D", "", name ) ) %>%
#   select( site = `Site Code`, age, salinity)
# 
# temperature <- meta %>%
#   pivot_longer( `Ave Temp - 30 days`:`Ave Temp - 90 days`,
#                 names_to = "name", values_to = "temperature_c") %>%
#   mutate( age = gsub( "\\D", "", name ) ) %>%
#   select( site = `Site Code`, age, temperature_c)
# 
# abiotic <- full_join( salinity, temperature)
# 
# # write temperature and salinity data to disk
# write_csv(abiotic, "csv/abiotic.csv")
# 
# meta <- meta %>%
#   select( site = `Site Code`, site_name = Site, Lat, Long,
#           `Date Range - 30 days`, `Date Range - 60 days`, `Date Range - 90 days`)
# # rename Santa Barabara site code from SBC to SBH
# meta <- meta %>%
#   mutate( site = gsub( pattern = "SBC", replacement = "SBH", site) )
# write metadata to disk
# write_csv( meta, "data/output/metadata.csv" )



#### Add functionality for ordering by richness or arranging by ocean basin
# ocean basin
d <- left_join( d, select(meta, site, ocean, long = Long))



# species/morphofunctional list for each panel
tlist <- read_csv("data/taxon_list.csv")
totalrich <- tlist %>%
  group_by( site ) %>%
  summarize( total_richness = length(unique(taxon)) )

d <- left_join(d, totalrich)

# recode sites after merging
d$site <- unlist( lapply( strsplit(d$site,"-"), function(z) z[2] ) )
d$site[d$site == "VAS"] <- "VCR"

# write to disk
write_csv(d, "data/output/data_long.csv")


# ------------------------------------------------------------------------------


# Calculate averages at site level
dsite <- d %>% 
  group_by(site, lat, ocean, age, temp, salinity, total_richness) %>% 
  summarize( ar_bryo = mean(ar_bryo, na.rm=T), col_asc = mean(col_asc, na.rm=T), sol_asc = mean(col_asc, na.rm=T), sabellids = mean(sabellids, na.rm=T),
             richness = mean(richness, na.rm=T), shannon = mean(shannon, na.rm=T), 
             mfrichness = mean(mfrichness),
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
  select( site, age, lat, richness, total_richness , shannon, mfrichness, rug2, temp, salinity, total_cover, ar_bryo ) %>% 
  group_by(site, age, lat, total_richness ) %>% 
  summarise( temp = mean(temp, na.rm = T ), salinity = mean(salinity, na.rm = T), 
             richness = mean(richness, na.rm=T ), shannon = mean(shannon, na.rm=T), 
             mfrichness = mean(mfrichness, na.rm=T),
             rug2 = mean(rug2, na.rm=T), total_cover = mean(total_cover),
             ar_bryo = mean(ar_bryo, na.rm = T)) %>% 
  mutate( logrug = log( rug2 ) )

write_csv(dsite, "data/output/data_site.csv")


# ------------------------------------------------------------------------------
## Make the data wide

# get average temps and keep the rest
dsitemean <- dsite %>% 
  group_by(site) %>% 
  summarize( temp_mean = mean(temp), sal_mean = mean(salinity) )

dsite <- left_join( dsite, dsitemean )


dwide <-  dsite %>%  
  select( site, total_richness, age, temp_mean, sal_mean, richness, mfrichness, logrug, total_cover, ar_bryo) %>% 
  pivot_wider( names_from = age, values_from = c(richness, mfrichness, logrug, total_cover, ar_bryo))

# write to disk for other analyses
write_csv( dwide, "data/output/data_wide.csv" )




# # ------------------------------------------------------------------------------
# 
# 
# # calculate total cover
# total_cover <- rowSums(comm)
# 
# 
# # combine total cover with open space
# cover <- data.frame( comm_select[1:3], total_cover, open_space = comm_select$open_space)
# cover$site <- unlist( lapply( strsplit(cover$site,"-"), function(z) z[2] ) )
# 
# # write to disk
# write_csv( cover, "data/output/pcover_total.csv" )
# read_csv("data/output/pcover_total.csv")

# ------------------------------------------------------------------------------


# 
# # merge
# dsel <- d %>% select(panel, site, age, rugosity = rug2, richness, sh_diversity)
# cover <- left_join( cover, dsel  )
# 


# ------------------------------------------------------------------------------
## estimate rate of change of visual percent cover as slopes from regressions
#  units = percent cover per day


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

# repeat, omitting day 90 data to get an initial growth rate
dzeros90 <- dzeros %>% filter( age < 90 )
regressions2 <- dzeros90 %>%
  nest(data = -site) %>%
  mutate(
    fit = map(data, ~ lm(total_cover ~ age + 0, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

slopes2 <- regressions2 %>%
  unnest(tidied) %>% 
  filter(term == "age") %>% 
  select( site, estimate )

# repeat, omitting days 60 and 90 data to get an initial growth rate
dzeros30 <- dzeros %>% filter( age < 60 )
regressions3 <- dzeros30 %>%
  nest(data = -site) %>%
  mutate(
    fit = map(data, ~ lm(total_cover ~ age + 0, data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

slopes3 <- regressions3 %>%
  unnest(tidied) %>% 
  filter(term == "age") %>% 
  select( site, estimate )

# use quasibinomial models
dzeros$cover_prop <- dzeros$total_cover/100
dzeros90$cover_prop <- dzeros90$total_cover/100
glms <- dzeros %>%
  nest(data = -site) %>%
  mutate(
    fit = map(data, ~ glm(cover_prop ~ age, family = quasibinomial(), data = .x)),
    tidied = map(fit, tidy),
    glanced = map(fit, glance),
    augmented = map(fit, augment)
  )

slopes_glm <- glms %>%
  unnest(tidied) %>% 
  filter(term == "age") %>% 
  select( site, estimate )
# hist((slopes_glm$estimate))


# merge the two slope estimates
slopes$method <- "lm_all"; slopes2$method <- "lm_middle"; slopes3$method = "lm_initial"; slopes_glm$method <- "glm" 
slopes <- bind_rows(slopes,slopes2,slopes3,slopes_glm)
slopes <- slopes %>% pivot_wider( names_from = method, values_from = estimate )
plot(lm_all ~ lm_initial, data = slopes); abline(a = 0, b = 1)
plot(glm ~ lm_initial, data = slopes); abline(a = 0, b = 1)

# write to disk
write_csv(slopes,"data/output/cover_rate_slopes.csv")
# psych::pairs.panels(slopes[-1], breaks = 10)





# ------------------------------------------------------------------------------
# prepare data used for SEM




# remove NA values for relevant variables
dpick = dwide %>% filter( ! is.na(logrug_90), ! is.na(total_cover_90)  )


# merge
dslope <- left_join( dpick, slopes )

# # species list to calculate pooled site-level species richness (all species observed during the study at each site)
# tlist <- read_csv("data/taxon_list.csv")
# totalrich <- tlist %>%
#   group_by( site ) %>%
#   summarize( total_richness = length(unique(taxon)) )
# totalrich$site <- unlist( lapply( strsplit(totalrich$site,"-"), function(z) z[2] ) )

# # add mean morphofunctional richness and mean species richness for each site
# d <- d %>%
#   group_by(site) %>%
#   summarize( mfrichness = mean(mfrichness),
#              richness = mean(richness))
# # merge
# drich <- left_join(drich,d)

# write to disk
write_csv(dslope, "data/output/data_sem.csv")
