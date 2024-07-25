### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# this script uses percent cover data from the panels
# to calculate total cover and relative cover 

# packages
library(tidyverse)
library(readxl)
library(viridis)

# read the data
cover_raw <- read_excel("data/PCover_taxgroups.xlsx", sheet = 1 )
# misspelling
cover_raw$Age[ cover_raw$Age %in% "90D" ] <- "90d"
# numeric age in days
cover_raw$age_num <- as.numeric(gsub("([0-9]+).*$", "\\1", cover_raw$Age))


# read metadata
meta <- read_csv("data/metadata.csv")
names(meta)

# add Ocean Basin
meta$ocean <- ifelse( meta$Long < -98, "Pacific", "Atlantic") # this works even for Panama and Texas sites

# rename Dauphin Island site based on the Pier site, which seems to be the one listed in the data summaries
meta$site[ meta$site == "USA-ALD-PIER"] <- "USA-ALD"



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
# # write metadata to disk
write_csv( meta, "data/metadata.csv" )


# # write cover data to file, so it can be used for other community analyses
# cover_save <- cover %>% select(Panel_ID:Forams)
# write_csv(cover_save, "csv/community_cover.csv")

# extract the community data set 
comm <- cover_raw %>% select( algae:sponge) %>% select(-open_space)
# calculate total cover
total_cover <- rowSums(comm)


# combine total cover with open space
cover <- data.frame( cover_raw[1:3],age_num = cover_raw$age_num, total_cover, open_space = cover_raw$open_space)
cover <- cover %>% 
  select( panel = Panel, site = Site, age = Age, age_num,
          total_cover, open_space)

# write to disk
write_csv( cover, "data/pcover_total.csv" )

# add latitude and longitude
cover <- left_join( cover, meta )

# get rugosity
d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
names(d) <- tolower(names(d))
d$age_num <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
d$panel <- gsub( "90D","90d", d$panel)
d$panel <- gsub( "60D","60d", d$panel)
d <- d %>% mutate( rugosity_raw = rugosity, rug1 = 1/rugosity_raw, rug2 = 1-rugosity_raw) %>% 
  mutate( age_factor = factor(age))

# merge
dsel <- d %>% select(panel, site, age, rugosity = rug2, richness, sh_diversity)
cover <- left_join( cover, dsel  )


## plotting data

# cover over time
ggplot( cover, aes( x = age, y = total_cover)) + 
  geom_point(alpha = 0.25) + geom_smooth( method = "glm", method.args = list(family="binomial"), se = T)
# outliers are USA-WAS

ggplot( cover, aes( x = age_num, y = total_cover, fill = log(rugosity))) + 
  facet_wrap( ~ site) +
  geom_smooth( aes(group = site), method = "lm", se = F, col = 'darkgrey', lwd = 0.75) + 
  geom_point(alpha = 1, size = 2, pch = 21) + 
  scale_fill_viridis(option = "D", direction = 1) +
  scale_x_continuous(name = "Panel age (days)", breaks = c(30,60,90), limits = c(25,95)) +
  scale_y_continuous(name = "Total % cover", breaks = c(0,50,100), limits = c(-5,105))


ggplot( cover, aes( x = total_cover, y = log(rugosity), fill = age)) + 
  facet_wrap( ~ site) +
  geom_smooth( aes(group = site), method = "lm", se = F, col = 'darkgrey', lwd = 0.75) + 
  geom_point(alpha = 1, size = 2, pch = 21) 
  # scale_fill_viridis(option = "D", direction = -1) +
  # scale_x_continuous(name = "Panel age (days)", breaks = c(30,60,90), limits = c(25,95)) +
  # scale_y_continuous(name = "Total % cover", breaks = c(0,50,100), limits = c(-5,105))

ggplot( cover, aes( x = sh_diversity, y = log(rugosity), fill = total_cover)) + 
  facet_wrap( ~ site) +
  geom_smooth( aes(group = site), method = "lm", se = F, col = 'darkgrey', lwd = 0.75) + 
  geom_point(alpha = 1, size = 2, pch = 21)

psych::pairs.panels( cover %>% select(age,total_cover,rugosity))

# order based on latitude
# order based on final total cover
cover90 <- cover %>% filter(age == 90) %>% group_by(site) %>% 
  summarize( mean_cover_90 = mean(total_cover) )
cover30 <- cover %>% filter(age == 30) %>% group_by(site) %>% 
  summarize( mean_cover_30 = mean(total_cover) )
cover <- left_join(left_join(cover, cover90), cover30)
cover <- cover %>% mutate( site_order = fct_reorder2(site, mean_cover_90, -mean_cover_30) )
ggplot( cover, aes( x = age, y = cover01, group = site)) + 
  facet_wrap(~site_order) + 
  geom_point() + geom_smooth( method = "glm", method.args = list(family="binomial"), se = T)
# ggsave("figs/cover_age_site.png", width = 6, height = 4)

# add zero time point data
cover_meta <- cover %>% select(site, site_name, site_order, Lat, Long) %>% distinct()
cover0 <- cover_meta %>% mutate( age = 0, total_cover = 0, open_space = 100, cover01 = 0 )
coverall <- bind_rows(cover, cover0)
ggplot( coverall, aes( x = age, y = cover01, group = site)) + 
  facet_wrap(~site_order) + 
  geom_point() + geom_smooth(se = F) + 
  geom_smooth( method = "glm", method.args = list(family="binomial"), 
               se = F, col = "black", lwd = 0.75 )
# ggsave("figs/cover_age_site.png", width = 6, height = 4)

