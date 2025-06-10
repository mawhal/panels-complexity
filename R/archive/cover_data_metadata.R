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

# read main tidied data
d <- read_csv("data/output/data_long.csv")
# read metadata
meta <- read_csv("data/output/metadata.csv")

# read cover data
cover <- read_csv("data/output/pcover_total.csv")

# add latitude and longitude
cover <- left_join( cover, meta )

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

