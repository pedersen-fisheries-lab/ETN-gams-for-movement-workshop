library(dplyr)
library(lubridate)

#data is from Freitas, et al. 2021 "Sea temperature effects on depth use and
#habitat selection in a marine fish community", and is shared based on the
#conditions of the Creative Commons CC0 1.0 Universal license

#data downloaded from: https://zenodo.org/record/4689686#.YoVlGahBziD 
#paper: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2656.13497

cod_move <-  read.csv(unz("data/raw-cod-data/VPS_locations_cod__pollack__ballan_10_best.zip","VPS_locations_cod__pollack__ballan_10_best.csv"))%>%
  filter(species =="Cod")

cod_counts <- cod_move%>%
  group_by(DETECTEDID)%>%
  summarize(n=n())%>%
  arrange(desc(n))

seconds_per_day <- 24*60*60

cod_filtered <- cod_move%>%
  filter(DETECTEDID %in% cod_counts$DETECTEDID[1:4])%>%
  mutate(datetime = lubridate::as_datetime(DATETIME))%>%
  group_by(DETECTEDID)%>%
  mutate(days_from_tag = time_length(datetime - min(datetime),unit="days"),
         days_from_tag = as.numeric(days_from_tag),
         week = as.numeric(week(datetime)))%>%
  dplyr::select(datetime,
                days_from_tag, 
                jul,
                week,
                DETECTEDID,
                x,
                y,
                LON,
                LAT,
                DEPTH, 
                Bottomdepth,
                Habitat, 
                T1m:T33m)%>%
  dplyr::rename(x = x,
                y= y, 
                lon= LON, 
                lat = LAT,
         indiv = DETECTEDID,
         depth_m = DEPTH,
         bdepth_m = Bottomdepth,
         habitat = Habitat,
         juldate = jul)%>%
  mutate(depth_m = ifelse(depth_m<0,-depth_m, depth_m),
         indiv   = factor(indiv))%>%
  #remove a couple very separated from the rest of observations from one fish
  filter(!(indiv=="A69-9002-011784"&days_from_tag>200))%>%
  group_by(indiv)%>%
  arrange(days_from_tag)%>%
  mutate(
    #this is useful for segregating observations into time series for
    #autoregressive modelling in bam
    start_index= c(TRUE, rep(FALSE,times= n()-1)),
    dist_from_bottom_m = bdepth_m - depth_m,
    #negative values will get re-coded as zeros
    dist_from_bottom_m = ifelse(dist_from_bottom_m<0,0,dist_from_bottom_m),
    speed_mps = sqrt((x-dplyr::lag(x))^2+(y-dplyr::lag(y))^2)/(days_from_tag-lag(days_from_tag))*(1/seconds_per_day))
    

write.csv(cod_filtered,file = "data/cod_move.csv")
  
