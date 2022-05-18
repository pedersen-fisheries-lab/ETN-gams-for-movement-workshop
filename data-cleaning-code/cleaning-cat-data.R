
## This code is adapted from Fletcher and Fortin 2019, 
## Spatial Ecology and Conservation Modeling: Applications with R
## Chapter 8: Space use and Resource Selection
## Code from: https://ufdc.ufl.edu/IR00010770/00001

#load packages

library(stars)           #for raster covariate data
library(starsExtra)
library(sf)
library(dplyr)         #for re-formatting data
library(tidyr)
library(ggplot2)
library(mgcv)
library(resample)


#Function to calculate the mean of a value for a group, excluding each oberservation
calc_mean_holdout <- function(x) {
  xbar <- jackknife(x,statistic = mean)$replicates[,1]
  return(xbar)
}


###################################################
#1. Prepping panther track data ####
###################################################

#load panther data
panther <- st_read("data/raw-panther-data/panthers.shp")%>%
  select(-X, -Y)%>%
  #Get UTM coordinates for each location from the shapefile
  bind_cols(st_coordinates(.))%>%
  mutate(CatID = factor(CatID))%>%
  group_by(CatID) %>% 
  mutate(
    is_obs = 1, 
    X_prior = lag(X),
    Y_prior = lag(Y),
    step_length = sqrt((X-X_prior)^2+ (Y-Y_prior)^2),
    Xbar_jack = calc_mean_holdout(X),
    Ybar_jack = calc_mean_holdout(Y))%>%
  #removing a single zero-step length observation
  filter(!step_length==0)%>%
  mutate(
    step = 0:(n()-1),
    abs_angle = atan2(y =  Y-Y_prior,x = X-X_prior),
    #deals with edge case where the first point does not have an angle
    abs_angle = ifelse(step==0, 0,abs_angle),
    abs_angle_prior = lag(abs_angle)
  )%>%
  #remove the first observation for each panther, as it has no info for step
  #selection
  filter(step !=0)

#used for plotting trajectory
panther_traj <- panther %>%
  summarize(do_union=FALSE) %>%
  st_cast("LINESTRING") 

#inspect
summary(panther)
levels(panther$CatID) #the unique cat IDs
head(panther)

#create as a shape object for plotting
panther_paths <- panther%>%
  group_by(CatID)%>%
  summarize()%>%
  st_cast("LINESTRING")

panther_crs <- st_crs(panther)


###########################################
#Find alternative steps
###########################################

#plot trajectories
ggplot(panther_traj)+
  geom_sf(aes(color=CatID))

set.seed(2022)

# Generate new steps: for each step, generate 8 angles around the compass rose,
# and four log-scaled distances
n_angles <- 8
n_steplengths <- 4
panther_altsteps <- panther %>%
  select(CatID, step,Juldate,AgeClass,X_prior,Y_prior,Xbar_jack,Ybar_jack, abs_angle_prior)%>%
  #create a combination of new angles and step lengths for each cat
  expand_grid(is_obs =0,
              abs_angle = 2*pi*(0:(n_angles-1))/n_angles,
              step_length = 10^seq(1,4, length=n_steplengths))%>%
  group_by(CatID,step_length)%>%
  mutate(
    #randomly rotate the sampled angles for each distance class to increase
    #coverage of potential turn angles 
    abs_angle = abs_angle + runif(1, 0,2*pi),
    abs_angle = abs_angle%/%(2*pi),
    #randomly perturb the step lengths to increase coverage of the step-length
    #distribution sampled
    step_length = step_length*rlnorm(n(), 0,0.1))%>%
  ungroup()%>%
  #find the physical location of each step
  mutate(X = X_prior + step_length*cos(abs_angle),
         Y = Y_prior + step_length*sin(abs_angle))

panther_steps <-panther %>%
  #remove the spatial geometry from the prior data
  st_set_geometry(NULL)%>%
  #combine the two data sets
  bind_rows(panther_altsteps)%>%
  #calculate relative angles of each turn
  mutate(rel_angle = abs_angle - abs_angle_prior,
         rel_angle = case_when(rel_angle >  pi ~ rel_angle - 2*pi,
                               rel_angle < -pi ~ rel_angle + 2*pi,
                               step == 1 ~ 0,
                               TRUE~rel_angle),
         #this finds how far each step or alternative step is to the average
         #center of the observations for the cat *excluding the current step
         #from the mean*. This is necessary because otherwise there would be a
         #bias towards points closer to the mean, as a point that is far away
         #from the mean will pull the mean toward it
         dist_to_center = sqrt((X-Xbar_jack)^2 +(Y-Ybar_jack)^2),
         #need step stratum to be distinct number for each combination of cat
         #and step
         step_stratum = factor(paste(step, CatID, sep = "-")),
         step_stratum = as.numeric(step_stratum),
         #need a constant time variable for the Cox PH model
         times = 1)%>%
  #turn back into a spatial data set
  st_as_sf(coords = c("X","Y"),
           remove=FALSE,
           crs = panther_crs)

test_model <- gam(cbind(times, step_stratum)~s(log10(step_length),k=5)+s(rel_angle,bs="cc"), 
                  data=panther_steps,
                  family =cox.ph, 
                  weights= panther_steps$is_obs)

## Get environmental covariates ####


#landcover source: fwc/fnai
land <- stars::read_stars("data/raw-panther-data/panther_landcover.grd")%>%
  rename(cover_type = panther_landcover.grd)%>%
  st_transform_proj(crs = panther_crs)%>%
  #reduce the size of the map
  st_downsample(n = 2)


#label projection for later use
land_crs <- st_crs(land)

land_crs == panther_crs

#plot
ggplot()+
  geom_stars(data=land)+
  geom_sf(data=panther, aes(col=CatID))

#load reclassification table for reclassifying map
classification <- read.table("data/raw-panther-data/landcover reclass.txt", header=TRUE,stringsAsFactors = TRUE)%>%
  arrange(Landcover)

#inspect
head(classification)
levels(classification$Description)    #original classification
levels(classification$Description2)   #re-class


#format for reclassify function;
class <- as.matrix(classification[,c(1,3)])

land_sub <- cut(land,c(0,class[,1]),labels = class[,2])



#5 km moving window to get neighborhood proportion since the grid size is 1000 m,
#this gives a circle with a radius of 5 km around the focal location This is
#used for calculating a spatial average of habitat type
in_circle <- function(x,y,radius = 5) ifelse(sqrt(x^2+y^2) <= radius, 1,0 )

#uses the in_circle function to create the filter
fw  <- outer((-5:5), (-5:5),FUN = in_circle)

#scaling the filter to sum to 1
fw <- fw/sum(fw)

#plot the filter to check if it works:
image(fw)
dim(fw) == c(21,21)



dry_forest <- land %>%
  mutate(cover_type = case_when(cover_type==10|cover_type==12~1,
                                TRUE ~ 0))%>%
  focal2(x = ., w= fw,fun = "sum",na.rm = TRUE)%>%
  rename(dryforest = cover_type)

wet_forest <- land %>%
  mutate(cover_type = case_when(cover_type==9|cover_type==11~1,
                                TRUE ~ 0))%>%
  focal2(x = ., w= fw,fun = "sum",na.rm = TRUE)%>%
  rename(wetforest = cover_type)


#merge into a single raster stack
layers <- c(land_sub, wet_forest, dry_forest)


landscape_class_map <- ggplot() + 
  geom_stars(data= layers, aes(fill = cover_type))+
  scale_fill_discrete(breaks= classification$ChangeTo,
                      labels= classification$Description2) +
  geom_sf(data = panther_traj, aes(group=CatID))+
  labs(title= "All landscape classes")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom")


wetforest_map <- ggplot() + 
  geom_stars(data= layers, aes(fill = wetforest))+
  scale_fill_viridis_c() +
  geom_sf(data = panther_traj, aes(group=CatID,color=CatID))+
  labs(title= "Average abundance of wet forest\nwithin 5 km")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom")



dryforest_map <- ggplot() + 
  geom_stars(data= layers, aes(fill = dryforest))+
  scale_fill_viridis_c(option = "magma")+
  geom_sf(data = panther_traj, aes(group=CatID,color=CatID))+
  labs(title= "Average abundance of dry forest\nwithin 5 km")+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom")




classification_labels <- classification%>%
  transmute(cover_type = factor(ChangeTo),
         cover_label = Description2)%>%
  group_by(cover_type)%>%
  summarize( cover_label = cover_label[1])

#extracts covariates
panther_steps_withcovs <- st_join(panther_steps,
                                  st_as_sf(layers,as_points = FALSE))%>%
  left_join(classification_labels)%>%
  #filtering out observations occurring in open water
  filter(!cover_type==10)%>%
  #making data nonspatial to save more easily
  st_set_geometry(NULL)%>%
  #filter out any alternative steps where the true step was filtered
  group_by(  step_stratum )%>%
  filter(any(is_obs==1))

write.csv( panther_steps_withcovs, "data/panther_stepsel.csv",row.names = FALSE)

ggsave(filename = "figures/habitat_map.png",plot = landscape_class_map,width = 8,height=9,units = "in",dpi = 400)
ggsave(filename = "figures/dryforest_map.png",plot = dryforest_map,width = 8,height=9,units = "in",dpi = 400)
ggsave(filename = "figures/wetforest_map.png",plot = wetforest_map,width = 8,height=9,units = "in",dpi = 400)


