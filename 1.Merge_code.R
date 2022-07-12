rm(list=ls())
library(ncdf4)
library(matrixStats)
library(qlcMatrix)
library(rstudioapi)

#---> 1. Get the data

# get function
get_data <- function(ncin,dname){
  data <- ncvar_get(ncin,dname)
  fillvalue <- ncatt_get(ncin,dname,"_FillValue") # NA values
  data[data==fillvalue$value] <- NA
  return(data)
}

## select algorithm
method <- "tesellation"
stat_n <- "mean"

# nc file location
#dir <- paste('Z:/♬ 나만의 공간_업뎃20190819/#8.  JY_연구/1.JY_연구/2.환경과제_강상욱교수님/2021_new/알고리즘분석중/20211201_IDW_Tesell차이점/data/',method,'/',sep="")
cur_dir = dirname(getSourceEditorContext()$path)
dir = cur_dir 
len <- list.files(path=dir)
dir
setwd(cur_dir)
saveloc = cur_dir
list.files(path=dir)
# put nc files on the list 
NO2_nc <- list()
lat_nc = list()
#for(i in 1:length(len)){
for(i in 11:17){
  dir_f <- paste(dir,len[i],sep="/")
  NO2_nc[[i-10]] <- get_data(nc_open(dir_f), "Data Fields/ColumnAmountNO2") # NO2
  cat(i,"/",length(len),'\n')
}
for(i in 11:17){
  dir_f <- paste(dir,len[i],sep="/")
  lat_nc[[i-10]] <- get_data(nc_open(dir_f), "Geolocation Fields/Latitude") # NO2
  cat(i,"/",length(len),'\n')
}
for(i in 11:17){
  dir_f <- paste(dir,len[i],sep="/")
  lon_nc[[i-10]] <- get_data(nc_open(dir_f), "Geolocation Fields/Longitutud") # NO2
  cat(i,"/",length(len),'\n')
}
dir_f
p =get_data(nc_open(dir_f), "Geolocation Fields/Latitude")

#---> 2. Calculate the average and count missing
NO2_dat <- NULL
len <- length(NO2_nc)

for(k in 1:len){
  NO2_dat <- cbind(NO2_dat,as.vector((as.matrix(NO2_nc[[k]])))) # merge NO2
}


## NO2 cal
NO2_dat_VV <- apply(NO2_dat,1,mean ,na.rm=T)
NO2_dat_VV[which(NO2_dat_VV==-Inf)] <- NA
NO2_dat_VV[which(NO2_dat_VV==Inf)] <- NA
NO2_dat_VV <- matrix(NO2_dat_VV,ncol=1040)
NO2_dat_NA <- rowSums(is.nan(NO2_dat))/dim(NO2_dat)[2]
#mean(NO2_dat_NA) # 0.4988943 for IDW
mean(NO2_dat_NA) # 0.5158957 for Tessel

#---> 3. convert to nc file

# lat, lon 
degree <- 0.05 # grid degree
lon_min <- 80;lon_max <- 146
lat_min <- -6;lat_max <- 46
lon <- seq(lon_min+(degree/2),lon_max-(degree/2),by=degree)
lat <- seq(lat_min+(degree/2),lat_max-(degree/2),by=degree)
lon_mat <- matrix(rep(sort(lon,decreasing = T),length(lat)),nrow=length(lat),byrow=T)
lat_mat <- matrix(rep(sort(lat,decreasing = T),length(lon)),ncol=length(lon),byrow=F)
lat <- get_data(nc_open(dir_f), "Geolocation Fields/Latitude") # latitude
lon <- get_data(nc_open(dir_f), "Geolocation Fields/Longitude")
# Save .nc
#path2 <- 'Z:/♬ 나만의 공간_업뎃20190819/#8.  JY_연구/1.JY_연구/2.환경과제_강상욱교수님/2021_new/알고리즘분석중/20211201_IDW_Tesell차이점/Result/'
#setwd(path2)
image <- ncdim_def("image","image",seq(1,dim(lat_mat)[1]))
spatial <- ncdim_def("spatial","spatial",seq(1,dim(lat_mat)[2]))

fill_value <- NA
lat_def <- ncvar_def(name="Geolocation Fields/Latitude",
                     units="degrees_north",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")
lon_def <- ncvar_def(name="Geolocation Fields/Longitude",
                     units="degrees_north",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")
no2_def <- ncvar_def(name="Data Fields/ColumnAmountNO2",
                     units="Molec per cm^2",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")
no2_NA <- ncvar_def(name="Data Fields/ColumnAmountNO2_NA",
                     units="percent",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")


ncfname <- paste("GK2_GEMS_L3_2021","03","_NO2_",method,"_",stat_n,"_day3.nc",sep="")
ncout <- nc_create(ncfname,list(lon_def,lat_def,no2_def,no2_NA),force_v4=TRUE)
ncvar_put(ncout,lon_def,lon)#as.matrix(lon_mat))
ncvar_put(ncout,lat_def,lat)#as.matrix(lat_mat))
ncvar_put(ncout,no2_def,as.matrix(NO2_dat_VV))
ncvar_put(ncout,no2_NA,as.matrix(NO2_dat_NA))
nc_close(ncout)


###############################################################################
#######################  week day and weekend ##################################
###############################################################################
rm(list=ls())
library(ncdf4)

#---> 1. Get the data

# get function
get_data <- function(ncin,dname){
  data <- ncvar_get(ncin,dname)
  fillvalue <- ncatt_get(ncin,dname,"_FillValue") # NA values
  data[data==fillvalue$value] <- NA
  return(data)
}

## select algorithm
method <- "idw"
stat_n <- "mean"

# nc file location
dir <- paste('Z:/♬ 나만의 공간_업뎃20190819/#8.  JY_연구/1.JY_연구/2.환경과제_강상욱교수님/2021_new/알고리즘분석중/20211201_IDW_Tesell차이점/data/',method,'/',sep="")
len <- list.files(path=dir)
len_vec <- c(which(substr(len,17,20)=="0306"),
which(substr(len,17,20)=="0307"),
which(substr(len,17,20)=="0313"),
which(substr(len,17,20)=="0314"),
which(substr(len,17,20)=="0320"),
which(substr(len,17,20)=="0321"),
which(substr(len,17,20)=="0327"),
which(substr(len,17,20)=="0328"))
len_weekend <- len[len_vec]
len_weekday <- len[-len_vec]

# dd <- "06" 
# path_d <- paste(dir,dd,"/",sep="")
# len <- list.files(path=path_d)
# len_vec <- c(which(substr(len,17,20)=="0605"),
#             which(substr(len,17,20)=="0606"),
#             which(substr(len,17,20)=="0612"),
#             which(substr(len,17,20)=="0613"),
#             which(substr(len,17,20)=="0619"),
#             which(substr(len,17,20)=="0620"),
#             which(substr(len,17,20)=="0626"),
#             which(substr(len,17,20)=="0627"))
# len_weekend <- len[len_vec]
# len_weekday <- len[-len_vec]

# dd <- "08" 
# path_d <- paste(dir,dd,"/",sep="")
# len <- list.files(path=path_d)
# len_vec <- c(which(substr(len,17,20)=="0807"),
#             which(substr(len,17,20)=="0808"),
#             which(substr(len,17,20)=="0814"),
#             which(substr(len,17,20)=="0815"),
#             which(substr(len,17,20)=="0821"),
#             which(substr(len,17,20)=="0822"),
#             which(substr(len,17,20)=="0828"),
#             which(substr(len,17,20)=="0829"))
# len_weekend <- len[len_vec]
# len_weekday <- len[-len_vec]


# put nc files on the list 
NO2_nc_weekend <- list()
NO2_nc_weekday <- list()

for(i in 1:length(len_weekend)){
  dir_f_weekend <- paste(dir,len_weekend[i],sep="")
  NO2_nc_weekend[[i]] <- get_data(nc_open(dir_f_weekend), "Data Fields/ColumnAmountNO2") # NO2
  cat(i,"/",length(len_weekend),'\n')
}

for(i in 1:length(len_weekday)){
  dir_f_weekday <- paste(dir,len_weekday[i],sep="")
  NO2_nc_weekday[[i]] <- get_data(nc_open(dir_f_weekday), "Data Fields/ColumnAmountNO2") # NO2
  cat(i,"/",length(len_weekday),'\n')
}

#---> 2. Calculate the average and count missing

NO2_dat_weekday <- NULL
NO2_dat_weekend <- NULL

for(k in 1:length(len_weekend)){
  NO2_dat_weekend <- cbind(NO2_dat_weekend,as.vector((as.matrix(NO2_nc_weekend[[k]])))) # merge NO2
}
for(k in 1:length(len_weekday)){
  NO2_dat_weekday <- cbind(NO2_dat_weekday,as.vector((as.matrix(NO2_nc_weekday[[k]])))) # merge NO2
}

## NO2 cal
NO2_dat_VV_weekend <- apply(NO2_dat_weekend,1,stat_n,na.rm=T)
NO2_dat_VV_weekend[which(NO2_dat_VV_weekend==-Inf)] <- NA
NO2_dat_VV_weekend[which(NO2_dat_VV_weekend==Inf)] <- NA

NO2_dat_VV_weekend <- matrix(NO2_dat_VV_weekend,ncol=320)
NO2_dat_NA_weekend <- rowSums(is.na(NO2_dat_weekend))/dim(NO2_dat_weekend)[2]

## NO2 cal
NO2_dat_VV_weekday <- apply(NO2_dat_weekday,1,stat_n,na.rm=T)
NO2_dat_VV_weekday[which(NO2_dat_VV_weekday==-Inf)] <- NA
NO2_dat_VV_weekday[which(NO2_dat_VV_weekday==Inf)] <- NA
NO2_dat_VV_weekday <- matrix(NO2_dat_VV_weekday,ncol=320)
NO2_dat_NA_weekday <- rowSums(is.na(NO2_dat_weekday))/dim(NO2_dat_weekday)[2]


#---> 3. convert to nc file

# lat, lon 
degree <- 0.05 # grid degree
lon_min <- 115;lon_max <- 131
lat_min <- 30;lat_max <- 43
lon <- seq(lon_min+(degree/2),lon_max-(degree/2),by=degree)
lat <- seq(lat_min+(degree/2),lat_max-(degree/2),by=degree)
lon_mat <- matrix(rep(sort(lon,decreasing = F),length(lat)),nrow=length(lat),byrow=T)
lat_mat <- matrix(rep(sort(lat,decreasing = T),length(lon)),ncol=length(lon),byrow=F)

# Save .nc
path2 <- 'Z:/♬ 나만의 공간_업뎃20190819/#8.  JY_연구/1.JY_연구/2.환경과제_강상욱교수님/2021_new/알고리즘분석중/20211201_IDW_Tesell차이점/Result/'
setwd(path2)
image <- ncdim_def("image","image",seq(1,dim(lat_mat)[1]))
spatial <- ncdim_def("spatial","spatial",seq(1,dim(lat_mat)[2]))

fill_value <- NA
lat_def <- ncvar_def(name="Geolocation Fields/Latitude",
                     units="degrees_north",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")
lon_def <- ncvar_def(name="Geolocation Fields/Longitude",
                     units="degrees_north",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")
no2_def <- ncvar_def(name="Data Fields/ColumnAmountNO2",
                     units="Molec per cm^2",
                     dim=list(image,spatial),
                     missval=fill_value,
                     prec="single")
no2_NA <- ncvar_def(name="Data Fields/ColumnAmountNO2_NA",
                    units="percent",
                    dim=list(image,spatial),
                    missval=fill_value,
                    prec="single")


ncfname <- paste("GK2_GEMS_L3_2021","03","_NO2_",method,"_",stat_n,"_weekend.nc",sep="")
ncout <- nc_create(ncfname,list(lon_def,lat_def,no2_def,no2_NA),force_v4=TRUE)
ncvar_put(ncout,lon_def,as.matrix(lon_mat))
ncvar_put(ncout,lat_def,as.matrix(lat_mat))
ncvar_put(ncout,no2_def,as.matrix(NO2_dat_VV_weekend))
ncvar_put(ncout,no2_NA,as.matrix(NO2_dat_NA_weekend))
nc_close(ncout)


ncfname <- paste("GK2_GEMS_L3_2021","03","_NO2_",method,"_",stat_n,"_weekday.nc",sep="")
ncout <- nc_create(ncfname,list(lon_def,lat_def,no2_def,no2_NA),force_v4=TRUE)
ncvar_put(ncout,lon_def,as.matrix(lon_mat))
ncvar_put(ncout,lat_def,as.matrix(lat_mat))
ncvar_put(ncout,no2_def,as.matrix(NO2_dat_VV_weekday))
ncvar_put(ncout,no2_NA,as.matrix(NO2_dat_NA_weekday))
nc_close(ncout)


