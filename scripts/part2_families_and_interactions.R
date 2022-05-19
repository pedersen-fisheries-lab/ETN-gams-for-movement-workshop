## 0. Getting set up ####
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gratia)

##1. New families with 1-d smoothers ####

##1.1 Modelling speed with a Gamma distribution ----


#data is from Freitas, et al. 2021 "Sea temperature effects on depth use and
#habitat selection in a marine fish community", and is shared based on the
#conditions of the Creative Commons CC0 1.0 Universal license. ETN co-author:
#David Villegas-Ríos. Thanks to the authors for making the code avaiable!

cod_move <- read.csv("data/cod_move.csv") %>%
  mutate(cod = case_when(indiv=="A69-9002-11803"~"A",
                         indiv=="A69-9002-11807"~"B",
                         indiv=="A69-9004-15547"~"C",
                         indiv=="A69-9004-15551"~"D"),
         cod = factor(cod),
         habitat = factor(habitat)
  )


View(cod_move)
#Key variables:
#- days_from_tag: fractional days since fish was tagged
#- depth_m: depth of the fish at detection time
#- bdepth_m: depth of the bottom at detection
#- T1m-t33m: Temperature at depth at the location: e.g. T1m = temp in degrees C at 1 m
#- dist_from_bottom_m: how far off the bottom the cod is


## Visualize speed distribution

cod_speed_hist <- ggplot(cod_move, aes(x= speed_mpmin))+
  geom_histogram(aes(y= after_stat(density)))+
  labs(x = "Speed (meters per minute)",
       y = "Density (1/m/min)")

print(cod_speed_hist)

print(cod_speed_hist + scale_x_log10())

## Let's regress speed on time since the start of the study for one cod

speed_mod1 <- gam(speed_mpmin~ s(days_from_tag),
                  data= filter(cod_move, cod=="A"),
                  family = Gamma(link="log"),
                  method = "REML")

plot(speed_mod1 )

## We can allow for more temporal variability by increasing the number of basis
## functions

speed_mod2 <- gam(speed_mpmin~ s(days_from_tag,k=50),
                  data= filter(cod_move, cod=="A"),
                  family = Gamma(link="log"),
                  method = "REML")

plot(speed_mod2)


##1.2 Modelling distance from bottom as a function of temperature with a Tweedie distribution ----



## Visualize distance from bottom distribution

cod_dist_hist <- ggplot(cod_move, aes(x= dist_from_bottom_m))+
  geom_histogram(aes(y= after_stat(density)))+
  labs(x = "Distance from bottom (m)",
       y = "Density (1/m)")

print(cod_dist_hist)

#Fraction of zero observations
print(mean(cod_move$dist_from_bottom_m==0))
#We need to deal with Zero-inflation!


## Modelling distance from bottom

bottom_dist_mod <- gam(dist_from_bottom_m ~ s(week, bs="cc") + s(cod, bs= "re"),
                       family =tw, data= cod_move)

plot(bottom_dist_mod)

##1.3 1D Step selection analysis with Cox PH family ---- 


panther_steps <- read.csv("data/panther_stepsel.csv")%>%
  mutate(habitat = factor(habitat))



#This model just fits a step-selection model to the step length and turn angle
#key elements here: c(bind(times, steps_stratum)): times will be a constant.
#step_stratum identifies which observations correspond to the observed step and
#alternative steps for a given step

#weights: corresponds to whether an "event" occurred or not: weight of 1 here
#corresponds to whether a given observation is a step or not
null_model <- gam(cbind(times, step_stratum)~s(log10(step_length),k=5)+s(rel_angle), 
                  data=panther_steps,
                  family =cox.ph, 
                  weights= panther_steps$is_obs)

summary(null_model)
plot(null_model,page=1, scale=0)


# We can add in a smoother for habitat type as a random effect, to estimate
# habitat preferences

habitatpref_model <- gam(cbind(times, step_stratum)~s(log10(step_length),k=5)+s(rel_angle)+
                      s(habitat,bs="re"), 
                    data=panther_steps,
                    family =cox.ph, 
                    weights= panther_steps$is_obs)



summary(habitatpref_model)
plot(habitatpref_model,page=1, scale=0)

pred_data <- tibble(habitat= unique(panther_steps$habitat),
                    step_length=1000,
                    rel_angle = 0
                    )
habitat_pred <- pred_data %>%
bind_cols(predict(habitat_pref,newdata = pred_data,type = "terms"))%>%
  select(habitat, `s(habitat)`)

print(habitat_pred) 

#model the tendency to return toward the center of the range 
#dist_to_center: distance in meters to the centroid of the observations 
#*excluding the given step*
homerange_model <- gam(cbind(times, step_stratum)~s(log10(step_length),k=5)+s(rel_angle)+
                         s(dist_to_center)+s(wetforest) + s(dryforest), 
                       data=panther_steps,
                       family =cox.ph, 
                       weights= panther_steps$is_obs)

summary(homerange_model)
plot(homerange_model,page=1, scale=0)



##2. Handling issues with big and autocorrelated data ####

##2.1 using `bam` for big data ----

# bam can handle big data sets with less memory overhead: 

speed_mod2_bam <- bam(speed_mpmin~ s(days_from_tag,k=50),
                  data= filter(cod_move, cod=="A"),
                  family = Gamma(link="log"),
                  method = "fREML")

plot(speed_mod2_bam)

#This does not give the exact same model as `gam` though! 

# You can also set `discrete = TRUE` to speed up calculations (not just memory
# usage). be careful, because this is still somewhat in development

speed_mod3_bam <- bam(speed_mpmin~ s(days_from_tag,k=50),
                      data= filter(cod_move, cod=="A"),
                      family = Gamma(link="log"),
                      discrete = TRUE,
                      method = "fREML")

plot(speed_mod3_bam)

## 2.2 Using a spatial smoother to account for dependency ----

#more about the 2D smoother in a bit!
bottom_dist_mod2 <- gam(dist_from_bottom_m ~ s(week, bs="cc") + 
                         s(cod, bs= "re")+ 
                         s(x,y, bs="tp"),
                       family =tw, data= cod_move)


# Compare the weekly smoother when accounting for space vs. not:

par(mfrow=c(1,2))

plot(bottom_dist_mod,select = 1)

plot(bottom_dist_mod2,select = 1)


## 2.3 Using `gamma` to account for too large a sample size ----

#With the step-selection data: we have many alternative turn angles and step
#lengths that we added; this isn't adding true new samples, so we need to
#account for this!

#The number of angles and step lengths
n_angles <- 8
n_steplengths <- 4

n_alternates <- n_angles*n_steplengths

#Larger gamma = each data point is downweighted
step_sel_gamma = n_alternates + 1


homerange_model2 <- gam(cbind(times, step_stratum)~s(log10(step_length),k=5)+s(rel_angle)+
                         s(dist_to_center)+s(wetforest) + s(dryforest), 
                       data=panther_steps,
                       family =cox.ph, 
                       gamma = step_sel_gamma,
                       weights= panther_steps$is_obs)

summary(homerange_model)
plot(homerange_model2,page=1, scale=0)

#compare with the original:
plot(homerange_model,page=1, scale=0)


