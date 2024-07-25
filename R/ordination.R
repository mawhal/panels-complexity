### Smithsonian MarineGEO PANELS project
# conducted in 2021
# lead for the project: Dean Janiak
###

# this script uses multivariate analyses to visualize and test hypotheses 
# about drivers of community structure and how structure might relate to complexity

library(tidyverse)
library(vegan)
library(readxl)


# read data
cover_raw <- read_excel("data/PCover_taxgroups.xlsx", sheet = 1 )

# extract the community data set 
comm <- cover_raw %>% select( algae:sponge) %>% select(-open_space)
comm_meta <- cover_raw[1:3]
names(comm_meta) <- tolower(names(comm_meta))

# meta data and complexity data
# meta <- read_csv("data/metadata.csv")

d <- readxl::read_xlsx("data/data_community.xlsx", sheet = "Sheet1")
names(d) <- tolower(names(d))
d$age <- as.numeric(gsub("([0-9]+).*$", "\\1", d$age))
d$panel <- gsub( "90D","90d", d$panel)
d$panel <- gsub( "60D","60d", d$panel)
# convert rugosity  to that large values are high rugosity
d <- d %>% mutate( rugosity = 1 - rugosity)
dselect <- d %>% select(panel, lat, temp, salinity, rugosity, richness )
comm_meta <- left_join(comm_meta, dselect)

# add ocean
meta_og <- read_csv("data/metadata.csv")
comm_meta <- left_join( comm_meta, meta_og %>% select(site,ocean,lon=Long) %>% distinct()  )

# remove miss
missing = which(is.na(comm_meta$rugosity))
comm_meta <- comm_meta[-missing,]
comm <- comm[-missing,]

# nmds


# capscale/dbrda
comm_meta$temperature <- comm_meta$temp
ord1 <- capscale( comm ~ rugosity + temperature + salinity, comm_meta, dist = "bray")
ord2 <- capscale( comm ~ rugosity, comm_meta, dist = "bray")
ord <- ord1
plot(ord)
anova(ord)
summary(ord)
taxa <- scores(ord, display = 'species')
envfit(ord, env = comm_meta[,c("rugosity")] )
scores2 <- as.data.frame(taxa[   order( taxa[,1], taxa[,2] ), ])
scores2 <- scores2 %>% mutate(taxon = rownames(scores2))  %>% arrange(-CAP1)
# 

# nicer plot

## get species vectors (can be thought of correlations with CAP axes?)
ord.v <- data.frame( family=row.names(ord$CCA$v), ord$CCA$v ) %>% arrange(CAP1)

# extract axes
nax <- 1:2
scaling = 0
sitescore <- scores(ord,choices = nax, scaling = scaling)$sites
taxascore <- scores(ord,choices = nax, scaling = scaling)$species







# customize an RDA figure
# points make sure rate points to the right
cap <- sitescore
# cap[,1] <- cap[,1] * vec.dir
capspec <- taxascore
# rates
cap <- data.frame(cap,comm_meta)
coef(ord) # Function coef will give the regression coefficients from centered environmental variables (constraints and conditions) to linear combination scores. The coefficients are for unstandardized environmental variables. The coefficients will be NA for aliased effects.
# proportion explained
RsquareAdj(ord)  # explains 20-25% of variation in consumption rate?
R2 <- eigenvals(ord)/sum(eigenvals(ord))
R2
summary(ord)
cummr2 <- R2
for(i in 2:length(eigenvals(ord))){
  cummr2[i] <- cummr2[i]+cummr2[i-1]  
}

# join with meta data
sr <- data.frame( comm_meta, sitescore )


library(grid)
library(ggrepel)

vec.sp <- envfit(ord, perm = 1000, env = comm_meta[,c("temperature","salinity","rugosity")] )
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)

# filter data for the big three taxa
taxhi <- as.data.frame(taxascore) %>% filter(rownames(taxascore) %in% c("ar_bryo","col_asc","en_bryo") )
taxhi$taxnames <- c("Arborescent\nbryozoans","Colonial\nascidians","Encrusting\nbryozoans")
ggplot( data=sr, aes(x = CAP1,y = CAP2) ) + 
  geom_point(size=2, alpha = 0.5) +
  geom_point(data = taxascore, aes(x = CAP1, y = CAP2), fill = "orange", pch = 21, size = 2.5) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=CAP1,y=0,yend=CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour="dodgerblue") + 
  geom_text_repel(data=vec.sp.df,aes(x = CAP1,y = CAP2,label = species),size=5, col = "dodgerblue") +
  geom_text_repel(data = taxhi, aes(x = CAP1, y = CAP2, label = taxnames), col = "darkorange") +
  theme_test() +
  coord_fixed()
ggsave("figs/capscale.svg", width = 5, height = 3)  

