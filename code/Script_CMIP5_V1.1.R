#########################################################
##CMIP5 output assessment
##created by MJ on Jul-26-2016
##Script to process downloaded nc data
#########################################################
##libraries
library(ncdf4)
library(spatstat)
library(fields)
library(raster)
require(RCurl)
library(lattice)
library(fields)
#########################################################
##Run wget command and download CMIP5 data

##Currently not working
setwd("~/Downloads/")
sourceDir <- getwd()
DatFiles <- list.files(path = sourceDir, pattern = "\\.sh")

for (i in 1:length(DatFiles))
{
  command <- paste("bash ", DatFiles[i], " -H")
  system(command)
}


#########################################################
##Process downloaded CMIP5 GPP data
setwd("~/Documents/PostDoc/CMIP5_assessment/data/")

#########################################################
##Obtain a polygon file for Australia region
##Also made available: lon and lat min and max

ncpath <- "~/Documents/PostDoc/GDAY/precip/raw/monthly/bom-rain_month-19000101-19001231.nc"
corNC <- nc_open(ncpath)

r <- raster(ncpath, varname = "rain_month")
#plot(r)

corDF <- rasterToPoints(r)
corDF <- as.data.frame(corDF, stringsAsFactors = F)
colnames(corDF) <- c("x","y","value")
corDF$value <- 1
coordinates(corDF)=~x+y
gridded(corDF) = T
r_processed <- raster(corDF)

poly <- rasterToPolygons(r_processed, dissolve=T)

lat.min <- min(corDF$y)
lat.max <- max(corDF$y)
lon.min <- min(corDF$x)
lon.max <- max(corDF$x)


#######################################  Historic #######################################  
##Check for # of files and prepare the directory

sourceDir <- paste(getwd(), "CMIP5_data/historic/", sep="/")
destDir <- paste(getwd(), "extracted/historic", sep="/")

if (file.exists(destDir)){
  destDir <- destDir
} else {
  dir.create(file.path(destDir))
  destDir <- destDir
}

DatFiles <- list.files(path = sourceDir, pattern = "\\.nc")
#########################################################
##Document spatial resolution into a summary table

output <- as.data.frame(matrix(nrow=length(DatFiles), ncol=3))
colnames(output) <- c("model","nlon","nlat")

output <- spatres()

output$lon_res <- round(360/output$nlon,2)
output$lat_res <- round(180/output$nlat,2)

write.table(output, paste(getwd(), "/extracted/historic_spatial_resolution.csv", sep=""),
            col.names=T, row.names=F, sep=",")

#########################################################
##gpp_Lmon_CMCC-CESM_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 1
CMIPraster()

#########################################################
##gpp_Lmon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc

##Read in data
i <- 2
CMIPslice3()


#########################################################
##gpp_Lmon_GFDL-ESM2M_historical_r1i1p1_200101-200512.nc

##Read in data
i <- 3
CMIPslice3()


#########################################################
##gpp_Lmon_GISS-E2-H_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 4
CMIPslice3()


#########################################################
##gpp_Lmon_GISS-E2-H-CC_historical_r1i1p1_195101-201012.nc

##Read in data
i <- 5
CMIPslice4()


#########################################################
##gpp_Lmon_GISS-E2-R_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 6
CMIPslice3()


#########################################################
##gpp_Lmon_GISS-E2-R-CC_historical_r1i1p1_200101-201012.nc

##Read in data
i <- 7
CMIPslice4()


#########################################################
##gpp_Lmon_HadCM3_historical_r1i1p1_198412-200512.nc

##Read in data
i <- 8
CMIPslice3()


#########################################################
##gpp_Lmon_HadGEM2-CC_historical_r1i1p1_198412-200511.nc

##Read in data
i <- 9
CMIPslice3()


#########################################################
##gpp_Lmon_HadGEM2-ES_historical_r1i1p1_198412-200511.nc

##Read in data
i <- 10
CMIPslice3()


#########################################################
##gpp_Lmon_inmcm4_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 11
CMIPslice3()


#########################################################
##gpp_Lmon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 12
CMIPslice3()


#########################################################
##gpp_Lmon_IPSL-CM5A-MR_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 13
CMIPslice3()


#########################################################
##gpp_Lmon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 14
CMIPslice3()


#########################################################
##gpp_Lmon_MIROC-ESM_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 15
CMIPslice3()

#########################################################
##gpp_Lmon_MIROC-ESM-CHEM_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 16
CMIPslice3()

#########################################################
##gpp_Lmon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc

##Read in data
i <- 17
CMIPslice3()

#########################################################
##save modis onto individual CMIP resolution



#########################################################

## Downscale CMIP5 to modis resolution
setwd("~/Documents/PostDoc/CMIP5_assessment/data/")

cmipDir <- paste(getwd(), "/extracted/historic", sep="")
modisDir <- "/Volumes/JIANG/data/MODIS_GPP/raw/"

destDir <- paste(getwd(), "/modis/historic",sep="")

cmipFiles <- list.files(path = cmipDir, pattern = "\\.csv")
modisFiles <- list.files(path = modisDir, pattern = "\\.nc")

modisToCmip()

#########################################################
##Compute % diff at gridded level for 5-yr mean data

setwd("~/Documents/PostDoc/CMIP5_assessment/data/")

modis <- read.table("~/Documents/PostDoc/GDAY/MODIS_GPP/extracted/MODIS_GPP_200101-200512.csv",
                    header=T,sep=",")

cmipsourceDir <- "~/Documents/PostDoc/CMIP5_assessment/data/evaluation/historic/consist_reso/"
DatFiles <- list.files(path = cmipsourceDir, pattern = "\\.csv")

destDir1 <- paste(getwd(), "/evaluation/historic/monthlypercent/", sep="")
destDir2 <- paste(getwd(), "/evaluation/historic/annualstats/", sep="")

destDir3 <- paste(getwd(), "/evaluation/historic/monthlyanomaly/", sep="")
plotdest <- paste(getwd(), "/evaluation/historic/plot/", sep="")

compareF()

#########################################################
##Global Koppen climate classification

inName <- "~/Documents/PostDoc/CMIP5_assessment/Koppen_climate/wcs"
ncDF <- nc_open(inName)
r <- raster(inName, varname="crs")
##


