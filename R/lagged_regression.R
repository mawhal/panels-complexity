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

library(tidyverse)


## read data
dwide <- read_csv("data/output/data_wide.csv")


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
ggplot( dwide, aes( x = logrug_90, y = richness_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_90, y = richness_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()



# mfrichness -> complexity
ggplot( dwide, aes( x = mfrichness_30, y = logrug_30 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = mfrichness_30, y = logrug_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = mfrichness_30, y = logrug_90 )) +
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = mfrichness_60, y = logrug_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = mfrichness_60, y = logrug_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = mfrichness_90, y = logrug_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()


# compexity -> mfrichness

ggplot( dwide, aes( x = logrug_30, y = mfrichness_30 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_30, y = mfrichness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_30, y = mfrichness_90 )) +
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_60, y = mfrichness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_90, y = mfrichness_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = logrug_90, y = mfrichness_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()

ggplot( dwide, aes( x = richness_30, y = mfrichness_30 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_30, y = mfrichness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_30, y = mfrichness_90 )) +
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_60, y = mfrichness_60 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_60, y = mfrichness_90 )) + 
  geom_smooth( method = 'lm' ) +
  geom_point()
ggplot( dwide, aes( x = richness_90, y = mfrichness_90 )) + 
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

# pairs plot for lagged responses
dpairs <- dwide[c(4:6,10:12)]
names(dpairs) <- c("Richness\n            day 30", "Richness\n            day 60","Richness\n            day 90",
                      "ln(Rugosity)\n            day 30","ln(Rugosity)\n            day 60", "ln(Rugosity)\n            day 90")

# windows(6,6)
psych::pairs.panels( dpairs, scale = T, ellipses = T, smooth = F, stars = F,
                     method = "pearson", 
                     hist.col = "lightcoral",
                     cex.cor = 1.75, cex = 1.5
                     # diag.panel = panel.hist 
                     )



# ## ### consider models with morphfunctional richness
# # linear models
# lm1 <- lm( logrug_30 ~ scale(mfrichness_30), data = dwide )
# lm2 <- lm( logrug_60 ~ scale(mfrichness_30), data = dwide )
# lm3 <- lm( logrug_90 ~ scale(mfrichness_30), data = dwide )
# lm4 <- lm( logrug_60 ~ scale(mfrichness_60), data = dwide )
# lm5 <- lm( logrug_90 ~ scale(mfrichness_60), data = dwide )
# lm6 <- lm( logrug_90 ~ scale(mfrichness_90), data = dwide )
# 
# lm7 <- lm( mfrichness_30 ~ scale(logrug_30), data = dwide )
# lm8 <- lm( mfrichness_60 ~ scale(logrug_30), data = dwide )
# lm9 <- lm( mfrichness_90 ~ scale(logrug_30), data = dwide )
# lm10 <- lm( mfrichness_60 ~ scale(logrug_60), data = dwide )
# lm11 <- lm( mfrichness_90 ~ scale(logrug_60), data = dwide )
# lm12 <- lm( mfrichness_90 ~ scale(logrug_90), data = dwide )
# 
# 
# l1 <- list( lm1, lm2, lm3, lm4, lm5, lm6 )
# l2 <- list( lm7, lm8, lm9, lm10, lm11, lm12 )
# allmods <- c(l1,l2)
# do.call( rbind, lapply(l1, function(z) summary(z)$coeff[c(2,4)] ))
# ests <- data.frame( name = paste0("lm",1:12), 
#                     do.call( rbind, lapply(allmods, function(z) summary(z)$coeff[c(2,4)] )),
#                     do.call( rbind, lapply(allmods, function(z) confint(z)[2,] )) 
# )
# names(ests) <- c("model","est","se","lcl","ucl")
# 
# 
# ## Prepare the figure
# ests$direction <- gl( 2, 6, labels = c("mfrichness->complexity","complexity->mfrichness"))
# ests$comparison <- rep(c("30-30","30-60","30-90","60-60","60-90","90-90"), 2)
# ests$focus_lag <- rep(c("30_lag0","60_lag1","90_lag2","60_lag0","90_lag1","90_lag0"), 2)
# ests$lag <- rep(c(0,1,2,0,1,0), 2)
# ests$focus <- rep(c("30 days","60 days","90 days","60 days","90 days","90 days"), 2)
# 
# ggplot(data = ests, aes(x = focus, y = est)) +
#   facet_grid(lag~direction, scales = "free") +
#   geom_hline( yintercept = 0, col = "orange" ) +
#   geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.25 ) +
#   geom_point() +
#   theme_classic() +
#   coord_flip()
