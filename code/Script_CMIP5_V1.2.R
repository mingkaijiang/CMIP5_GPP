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
library(utils)
#########################################################
##Obtain downloading path from raijin platform, then run script to donwload the data

setwd("~/Documents/PostDoc/CMIP5_assessment/data/CMIP5_data/historical")
myDF <- read.table(paste(getwd(), "/out_to_download.csv",sep=""),
                   header=T,sep=",")

myDF$file.url <- as.character(myDF$file.url)
myDF$ensemble <- as.character(myDF$ensemble)

downDF <- subset(myDF, ensemble == "r1i1p1")

for (i in 1:length(downDF$file.url))
{
  file.path <- downDF$file.url[i]
  tperiod <- substr(file.path, (nchar(file.path) - 15), nchar(file.path))
  
  t.start <- substr(file.path, (nchar(file.path) - 15), (nchar(file.path)-10))
  t.end <- substr(file.path, (nchar(file.path) - 8), (nchar(file.path)-3))
  t.s.yr <- 
  
  outName <- paste(downDF$var[i], downDF$mip_table[i],downDF$model[i],downDF$experiment[i], downDF$ensemble[i],
                   tperiod, sep="_")
  
  download.file(file.path, outName)
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

setwd("~/Documents/PostDoc/CMIP5_assessment/data/")

cmipDir <- paste(getwd(), "/extracted/historic", sep="")
modisDir <- "/Volumes/JIANG/data/MODIS_GPP/raw/"

destDir <- paste(getwd(), "/modis/historic",sep="")

cmipFiles <- list.files(path = cmipDir, pattern = "\\.csv")
modisFiles <- list.files(path = modisDir, pattern = "\\.nc")

modisToCmip()

#########################################################
##Compute % diff at gridded level for 5-yr mean data

##Model list
#model <- c("CMCC-CESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","GISS-E2-H-CC",
#           'GISS-E2-R','HadCM3',"HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR",
#           "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM","MPI-ESM-LR")

model <- c("CMCC-CESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","GISS-E2-H-CC",
           'GISS-E2-R','HadCM3',"HadGEM2-ES","inmcm4","IPSL-CM5A-LR",
           "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM","MPI-ESM-LR")

setwd("~/Documents/PostDoc/CMIP5_assessment/data/")

modisDir <- paste(getwd(), "/modis/historic", sep="")
cmipsourceDir <- paste(getwd(), "/extracted/historic",sep="")

cmipFiles <- list.files(path = cmipsourceDir, pattern = "\\.csv")
modisFiles <- list.files(path = modisDir, pattern = "\\.csv")


destDir1 <- paste(getwd(), "/evaluation/historic/monthlypercent/", sep="")
destDir2 <- paste(getwd(), "/evaluation/historic/annualstats/", sep="")
destDir3 <- paste(getwd(), "/evaluation/historic/monthlystats/", sep="")
plotdest <- paste(getwd(), "/evaluation/historic/plot/", sep="")

compareF()

#########################################################
##Prepare modis and cmip data for taylor diagram
setwd("~/Documents/PostDoc/CMIP5_assessment/data/")
modisDir <- paste(getwd(), "/modis/historic", sep="")
cmipsourceDir <- paste(getwd(), "/extracted/historic",sep="")

cmipFiles <- list.files(path = cmipsourceDir, pattern = "\\.csv")
modisFiles <- list.files(path = modisDir, pattern = "\\.csv")
destDir <- paste(getwd(), "/evaluation/historic/monthlystats/", sep="")

DataPrepForTaylor()

#########################################################
##Plot Taylor diagram using all data
##Note: CMIP individual datasets are at different spatial resolution
##Therefore, modis are at respective resolution
##Also, the taylor diagram compares monthly data, not only about spatial variability

setwd("~/Documents/PostDoc/CMIP5_assessment/data/")
plotdest <- paste(getwd(), "/evaluation/historic/plot/", sep="")

i <- 3

inName <- paste(sourceDir, "/evaluation/historic/monthlystats/CMIP5_models_time_series_data.csv",sep="")
plotDF <- read.table(inName, sep=",", header=T)

colorlist <- c(1:ncol(plotDF)-2)


pdf(paste(plotdest, "/Taylor_GPP_interannual.pdf", sep=""))
par(mar = c(5.1,5.1,4.1,4.1))

taylor.diagram(plotDF$modis, plotDF[,i], main = "Normalized Taylor Diagram of GPP Interannual variability",
               pch = 16, xlab = "Standard deviation", normalize=T, col = colorlist[i-2])

for (i in 4:ncol(plotDF))
{
  taylor.diagram(plotDF$modis, plotDF[i], normalize=T, main = NA, xlab = NA, ylab = NA, col = colorlist[i-2])
}

legend("topright", c("MODIS", model[i]), col = c("black", "red"), pch = c(1,16))



#########################################################
##Global Koppen climate classification

inName <- "~/Documents/PostDoc/CMIP5_assessment/Koppen_climate/wcs"
ncDF <- nc_open(inName)
r <- raster(inName, varname="crs")
##


