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

cod_filtered <- cod_move%>%
  filter(DETECTEDID %in% cod_counts$DETECTEDID[1:4])%>%
  mutate(datetime = lubridate::as_datetime(DATETIME))%>%
  group_by(DETECTEDID)%>%
  mutate(days_from_tag = time_length(datetime - min(datetime),unit="days"),
         days_from_tag = as.numeric(days_from_start),
         week = as.numeric(week(datetime)))%>%
  dplyr::select(datetime,
                days_from_tag, 
                jul,
                week,
                DETECTEDID,
                LON,
                LAT,
                DEPTH, 
                Bottomdepth,
                Habitat, 
                T1m:T33m)%>%
  dplyr::rename(lon= LON,
         lat = LAT,
         indiv = DETECTEDID,
         depth_m = DEPTH,
         bdepth_m = Bottomdepth,
         habitat = Habitat,
         juldate = jul)%>%
  mutate(depth_m = ifelse(depth_m<0,-depth_m, depth_m),
         indiv   = factor(indiv))%>%
  #remove a couple very separated observations from one fish
  filter(!(indiv=="A69-9002-011784"&days_from_start>200))%>%
  group_by(indiv)%>%
  mutate(start_index= c(TRUE, rep(FALSE,times= n()-1)))
  
