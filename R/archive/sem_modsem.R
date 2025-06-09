# SEM using quadratic terms for latent variables

# packages
library(tidyverse)
library(readxl)
library(lavaan)
library(modsem)

# load merged data
d <- read_csv("data/data_sem.csv")
# log-transformed arborescenct bryozoan
d$log_ar_bryo_90 <- log10( d$ar_bryo_90+1 )

# model
model <-'
# Measurement Mode
GROW  =~ lm_middle
ID    =~ log_ar_bryo_90
DIV   =~ richness_30 #+ mfrichness + total_richness
# FUNCT =~ mfrichness
THERM =~ temp_mean
OSMO  =~ sal_mean
COMPL =~ logrug_90
# Structural Model
DIV   ~ THERM + OSMO #+ THERM:THERM
GROW  ~ THERM + THERM:THERM
ID    ~ THERM + OSMO + GROW
COMPL  ~ DIV + GROW + ID
'
est_dca <- modsem(model, data = d)
summary(est_dca, fit.measures = T, standardized = T, rsquare = T)
# User Model versus Baseline Model:
#   
#   Comparative Fit Index (CFI)                    0.966
#   Tucker-Lewis Index (TLI)                       0.910