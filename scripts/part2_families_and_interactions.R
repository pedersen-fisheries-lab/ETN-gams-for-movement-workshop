library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gratia)
panther_steps <- read.csv("data/panther_stepsel.csv")

n_angles <- 8
n_steplengths <- 4
test_model <- gam(cbind(times, step_stratum)~s(log10(step_length),k=5)+s(dist_to_center, m=1), 
                  data=panther_steps,
                  family =cox.ph, 
                  gamma= n_angles*n_steplengths + 1, 
                  weights= panther_steps$is_obs)
