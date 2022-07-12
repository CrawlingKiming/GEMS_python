
# clear the environment
rm(list=ls()) 

#--------------------------------------#
#---- 1. Load the required library ----#
#--------------------------------------#

library(udunits2)
if(!require(ggplot2)) install.packages('ggplot2'); library(ggplot2)
if(!require(sf)) install.packages('sf'); library(sf)
if(!require(rworldmap)) install.packages('rworldmap'); library(rworldmap)
if(!require(rgeos)) install.packages('rgeos'); library(rgeos)
library(udunits2)
library(tidyr)
library(stringr)
library(rstudioapi)


#saveloc <- aa

cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
saveloc = cur_dir

cur_dir
#--------------------------#
#---- 3. Function part ----#
#--------------------------#

plot_save <- function(NO2,lat,lon,plot_name,plot_title="",scale=16,threshold=5){
  
  # range of Korean Peninsula : lon: 123~131 / lat: 30~43
  #lon_min <- 123;lon_max <- 131
  lon_min <- 80;lon_max <- 146
  lat_min <- -6;lat_max <- 46
  
  dir.create(saveloc,showWarnings = F)
  newmap <- getMap(resolution = "low")
  par(mar = c(0,0,0,0))
  
  koreamap <- st_as_sf(newmap) # import map and change file format.
  #class(koreamap)
  
  rainbow.color = c("#910000FF","#FF0000FF",'#FF4300FF', "#FF9900FF",'#FFFF00FF', "#CCFF00FF",
                    "#00FF66FF","#00FFFFFF", "#0066FFFF",'#0000FFFF', "#3300FFFF",'#00008BFF')
  
  # map before marking dots 
  p <- ggplot(data = koreamap) +
    xlab("") + ylab("") +
    geom_sf(color = "black", fill = "white") +
    coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
    theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5),
          panel.background = element_rect(fill = "aliceblue"))
  
  pointsize <- 1.5 # Size of the estimated NO2 point
  shapenum <-  15 # shape of the point
  
  # make dataset 
  df1 = data.frame(longitude = as.vector(as.matrix(lon)), 
                   latitude = as.vector(as.matrix(lat)),
                   NO2 = as.vector(as.matrix(NO2)))
  
  # remove NA
  df1 = df1[complete.cases(df1),] 
  #df1$NO2[df1$NO2 <= 0] <- NA 
  df1 = df1[complete.cases(df1),] 
  # change scale wrt scale (ex. 16 => 10^16)
  df1$NO2 = df1$NO2 / (10^scale)
  # apply threshold as maximum
  df1$NO2[df1$NO2 >= threshold] <- threshold 
  df2 = rbind(df1, c(-1,-1,0))
  df2 = rbind(c(1000,1000,threshold),df2)
  
  # all area of GOCI
  g.t=paste0('GEMS L3 NO2 - ',plot_title %>% substr(1,100),' UTC')
  
  p.plot <- p +
    geom_point(data = df2, aes(x = longitude, y = latitude, color = NO2), size = pointsize, shape = shapenum) +
    scale_colour_gradientn(colours = rainbow.color, trans="reverse", breaks=seq(0,threshold,by=.5)) +
    borders("world", xlim = c(lon_min,lon_max), ylim = c(lat_min,lat_max)) +
    theme(legend.position="bottom") + # location change is possible
    guides(color = guide_colourbar(title = bquote("["*10^.(scale)~molecules~cm^{-2}*"]"),
                                   title.hjust = 0.5,title.position = "top",
                                   barwidth = 25, barheight = 1, reverse=T))+
    ggtitle(g.t)+theme(plot.title=element_text(hjust=0.5))
  
  eval(parse(text = paste0("ggsave('",saveloc,"/",plot_name,".png',width=30,height=20,unit='cm')")))
}

#-----------------------#
#---- 3. plot part ----#
#----------------------#
get_data <- function(ncin,dname){
  data <- ncvar_get(ncin,dname)
  fillvalue <- ncatt_get(ncin,dname,"_FillValue") # NA values
  data[data==fillvalue$value] <- NA
  return(data)
}

plot_ncf <- function(ncfname,plot_name,saveloc){
  require(ncdf4)
  #ncfname="GK2_GEMS_L2_20210322_0145_NO2_FC_DPRO_ORI"
  #size_chr <- str_pad(str_remove(as.character(size),"\\."),3,side="right","0")
  
  ncin <- nc_open(paste(ncfname,".nc",sep=""))
  lat <- get_data(ncin, "Geolocation Fields/Latitude") # latitude
  lon <- get_data(ncin, "Geolocation Fields/Longitude") # longitude
  NO2 <- get_data(ncin, "Data Fields/ColumnAmountNO2") # NO2
  #NO2_NA <- get_data(ncin, "Data Fields/ColumnAmountNO2_NA") # NO2_NA
  
  plot_title <- str_replace(plot_name,"_"," : ")
  plot_save(NO2,lat,lon,str_c(plot_name),plot_title,scale=16,threshold=5)
  #plot_save1(NO2_NA,lat,lon,str_c(plot_name,"_",size_chr),plot_title,threshold=1)
}

plot_ncf(ncfname="GK2_GEMS_L3_20210321_2345_NO2_WGS84",plot_name="L3_202103_NO2_tessellation_212345",saveloc=aa)
# 2021/3
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_mean_all",plot_name="L3_202103_NO2_tessellation_mean_all",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_mean_weekday",plot_name="L3_202103_NO2_tessellation_mean_weekday",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_mean_weekend",plot_name="L3_202103_NO2_tessellation_mean_weekend",saveloc=aa)

plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_idw_mean_all",plot_name="L3_202103_NO2_idw_mean_all",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_idw_mean_weekday",plot_name="L3_202103_NO2_idw_mean_weekday",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_idw_mean_weekend",plot_name="L3_202103_NO2_idw_mean_weekend",saveloc=aa)

plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_max_all",plot_name="L3_202103_NO2_tessellation_max_all",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_max_weekday",plot_name="L3_202103_NO2_tessellation_max_weekday",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_max_weekend",plot_name="L3_202103_NO2_tessellation_max_weekend",saveloc=aa)

plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_min_all",plot_name="L3_202103_NO2_tessellation_min_all",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_min_weekday",plot_name="L3_202103_NO2_tessellation_min_weekday",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_min_weekend",plot_name="L3_202103_NO2_tessellation_min_weekend",saveloc=aa)

plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_median_all",plot_name="L3_202103_NO2_tessellation_median_all",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_median_weekday",plot_name="L3_202103_NO2_tessellation_median_weekday",saveloc=aa)
plot_ncf(ncfname="GK2_GEMS_L3_202103_NO2_tessellation_median_weekend",plot_name="L3_202103_NO2_tessellation_median_weekend",saveloc=aa)












