#########################################################
##CMIP5 output assessment
##created by MJ on Jul-26-2016
##Function to process downloaded nc data
#########################################################
##libraries
library(ncdf4)
library(spatstat)
library(fields)
library(raster)
require(RCurl)
library(lattice)
library(maps)
require(plotrix)
#########################################################
##Check for spatial resolution
spatres <- function()
{
  for (i in 1:length(DatFiles))
  {
    inName <- file.path(sourceDir, DatFiles[i], fsep = .Platform$file.sep)
    check <- nc_open(inName)
    tmp.array <- ncvar_get(check, "gpp")
    nlon <- dim(tmp.array)[1]
    nlat <- dim(tmp.array)[2]
    nc_close(check)
    
    outName <- DatFiles[i]
    outName <- sub("*\\.nc", "", outName)
    outNamef <- substr(outName, 1,(nchar(outName)-14))
    
    output[i,"model"] <- outNamef
    output[i, "nlon"] <- nlon
    output[i, "nlat"] <- nlat
  }
  return(output)
}

#########################################################

##Function to extract Australia (200101-200512 based on raster method)
#
CMIPraster <- function()
{
  inName <- file.path(sourceDir, DatFiles[i], fsep = .Platform$file.sep)
  check <- nc_open(inName)
  
  ##Obtain time coverage information
  t <- ncvar_get(check, "time")
  nt <- dim(t)
  tunits <- ncatt_get(check, "time", "units")
  t.start <- substr(tunits$value,12,21)
  d2001 <- as.numeric(as.Date("2001-01-01") - as.Date(t.start))
  d2005 <- as.numeric(as.Date("2005-12-31") - as.Date(t.start))
  t1 <- t - d2001
  nt1 <- length(t1[t1<0]) + 1
  t2 <- t - d2005
  nt2 <- length(t2[t2<0]) 
  
  ##Get lon & lat information
  lon <- ncvar_get(check, "lon")
  nlon <- dim(lon)
  lonlist <- lon[,1]
  lat <- ncvar_get(check, "lat")
  nlat <- dim(lat)
  latlist <- lat[1,]
  nc_close(check)
  
  ##Project onto right resolution and coordinates
  d <- brick(inName, varname = "gpp")
  n <- nlayers(d)
  outDF <- rasterToPoints(d)
  outDF <- as.data.frame(outDF, stringsAsFactors=F)
  outDF$x <- outDF$x * 3.75 - 3.75
  outDF$y <- -(outDF$y * 3.75 - 3.75) + 90
  outDF1 <- subset(outDF, x <= 180)
  outDF2 <- subset(outDF, x > 180)
  outDF2$x <- outDF2$x - 360
  outDF3 <- rbind(outDF1, outDF2)
  
  ##Extract for Australia, and for the period of 200101-200512
  outDF <- subset(outDF, x <= lon.max & x >= lon.min & y <= lat.max & y >= lat.min)
  outcor <- outDF[,1:2]
  outdat <- outDF[,(nt1+2):(nt2+2)]
  out <- as.data.frame(cbind(outcor, outdat))
  
  ##Save data as excel file
  outName <- DatFiles[i]
  outName <- sub("*\\.nc", "", outName)
  outNamef <- substr(outName, 1,(nchar(outName)-14))
  outName <- paste(outNamef, "_200101-200512.csv", sep="")
  
  write.table(out, paste(destDir, "/", outName, sep=""),
              row.names=F, col.names=T, sep=",") 
}


#########################################################
##Function to extract Australia (200101-200512 based on irregular point method)
CMIPslice <- function()
{
  inName <- file.path(sourceDir, DatFiles[i], fsep = .Platform$file.sep)
  check <- nc_open(inName)
  
  ##Get lon lat information
  lon <- ncvar_get(check, "lon")
  nlon <- dim(lon)
  lonlist <- lon
  lat <- ncvar_get(check, "lat")
  nlat <- dim(lat)
  latlist <- lat
  
  ##Obtain time coverage information
  t <- ncvar_get(check, "time")
  nt <- dim(t)
  tunits <- ncatt_get(check, "time", "units")
  t.start <- substr(tunits$value,12,21)
  d2001 <- as.numeric(as.Date("2001-01-01") - as.Date(t.start))
  d2005 <- as.numeric(as.Date("2005-12-31") - as.Date(t.start))
  t1 <- t - d2001
  nt1 <- length(t1[t1<0])
  t2 <- t - d2005
  nt2 <- length(t2[t2<0]) 
  
  ##create 3d matrix
  tmp.array <- ncvar_get(check, "gpp")
  nc_close(check)
  
  ##extract lon lat according to time j
  j <- 1
  
  #dim(tmp.array)
  sliceDF <- tmp.array[,,j]
  sliceDF <- t(sliceDF)
  
  DF <- im(sliceDF, xcol = lon, 
           yrow = lat) 
  
  DF <- as.data.frame(DF)
  outDF1 <- subset(DF, x <= 180)
  outDF2 <- subset(DF, x > 180)
  outDF2$x <- outDF2$x - 360
  DF <- rbind(outDF1, outDF2)
  
  for (j in 2:nt)
  {
    sliceDF <- tmp.array[,,j]
    sliceDF <- t(sliceDF)
    
    iDF <- im(sliceDF, xcol = lon, 
              yrow = lat) 
    
    iDF <- as.data.frame(iDF)
    outDF1 <- subset(iDF, x <= 180)
    outDF2 <- subset(iDF, x > 180)
    outDF2$x <- outDF2$x - 360
    iDF <- rbind(outDF1, outDF2)
    
    DF[,j+2] <- iDF[,3]
  }
  
  colnames(DF) <- c("x","y", t)
  
  ##Extract for Australia, and for the period of 200101-200512
  outDF <- subset(DF, x <= lon.max & x >= lon.min & y <= lat.max & y >= lat.min)
  outcor <- outDF[,1:2]
  outdat <- outDF[,(nt1+2):(nt2+2)]
  out <- as.data.frame(cbind(outcor, outdat))
  
  ##Save data as excel file
  outName <- DatFiles[i]
  outName <- sub("*\\.nc", "", outName)
  outNamef <- substr(outName, 1,(nchar(outName)-14))
  outName <- paste(outNamef, "_200101-200512.csv", sep="")
  
  write.table(out, paste(destDir, "/", outName, sep=""),
              row.names=F, col.names=T, sep=",")
}

#########################################################
##Function to extract Australia (200101-200512 based on irregular point method)
##inconsistent starting time
CMIPslice2 <- function()
{
  inName <- file.path(sourceDir, DatFiles[i], fsep = .Platform$file.sep)
  check <- nc_open(inName)
  
  ##Get lon lat information
  lon <- ncvar_get(check, "lon")
  nlon <- dim(lon)
  lonlist <- lon
  lat <- ncvar_get(check, "lat")
  nlat <- dim(lat)
  latlist <- lat
  
  ##Obtain time coverage information
  t <- ncvar_get(check, "time")
  nt <- dim(t)
  tunits <- ncatt_get(check, "time", "units")
  t.start <- substr(tunits$value,12,21)
  d2001 <- as.numeric(as.Date("2001-01-01") - as.Date(t.start))
  d2005 <- as.numeric(as.Date("2005-12-31") - as.Date(t.start))
  t1 <- t - d2001
  nt1 <- length(t1[t1<0])
  t2 <- t - d2005
  nt2 <- length(t2[t2<0]) 
  
  ##create 3d matrix
  tmp.array <- ncvar_get(check, "gpp")
  nc_close(check)
  
  ##extract lon lat according to time j
  j <- 1
  
  #dim(tmp.array)
  sliceDF <- tmp.array[,,j]
  sliceDF <- t(sliceDF)
  
  DF <- im(sliceDF, xcol = lon, 
           yrow = lat) 
  
  DF <- as.data.frame(DF)
  outDF1 <- subset(DF, x <= 180)
  outDF2 <- subset(DF, x > 180)
  outDF2$x <- outDF2$x - 360
  DF <- rbind(outDF1, outDF2)
  
  for (j in 2:nt)
  {
    sliceDF <- tmp.array[,,j]
    sliceDF <- t(sliceDF)
    
    iDF <- im(sliceDF, xcol = lon, 
              yrow = lat) 
    
    iDF <- as.data.frame(iDF)
    outDF1 <- subset(iDF, x <= 180)
    outDF2 <- subset(iDF, x > 180)
    outDF2$x <- outDF2$x - 360
    iDF <- rbind(outDF1, outDF2)
    
    DF[,j+2] <- iDF[,3]
  }
  
  colnames(DF) <- c("x","y", t)
  
  ##Extract for Australia, and for the period of 200101-200512
  outDF <- subset(DF, x <= lon.max & x >= lon.min & y <= lat.max & y >= lat.min)
  outcor <- outDF[,1:2]
  outdat <- outDF[,(nt1+3):(nt2+2)]
  out <- as.data.frame(cbind(outcor, outdat))
  
  ##Save data as excel file
  outName <- DatFiles[i]
  outName <- sub("*\\.nc", "", outName)
  outNamef <- substr(outName, 1,(nchar(outName)-14))
  outName <- paste(outNamef, "_200101-200512.csv", sep="")
  
  write.table(out, paste(destDir, "/", outName, sep=""),
              row.names=F, col.names=T, sep=",")
}

#########################################################
##Function to extract Australia (200101-200512 based on irregular point method)
##count # of months backwards for data ended by 200512

CMIPslice3 <- function()
{
  inName <- file.path(sourceDir, DatFiles[i], fsep = .Platform$file.sep)
  check <- nc_open(inName)
  
  ##Get lon lat information
  lon <- ncvar_get(check, "lon")
  nlon <- dim(lon)
  lonlist <- lon
  lat <- ncvar_get(check, "lat")
  nlat <- dim(lat)
  latlist <- lat
  
  ##Obtain time coverage information
  t <- ncvar_get(check, "time")
  nt <- dim(t)
  t.start <- nt - 59
  
  ##create 3d matrix
  tmp.array <- ncvar_get(check, "gpp")
  nc_close(check)
  
  ##extract lon lat according to time j
  j <- t.start
  
  #dim(tmp.array)
  sliceDF <- tmp.array[,,j]
  sliceDF <- t(sliceDF)
  
  DF <- im(sliceDF, xcol = lon, 
           yrow = lat) 
  
  DF <- as.data.frame(DF)
  outDF1 <- subset(DF, x <= 180)
  outDF2 <- subset(DF, x > 180)
  outDF2$x <- outDF2$x - 360
  DF <- rbind(outDF1, outDF2)
  
  for (j in (t.start+1):nt)
  {
    sliceDF <- tmp.array[,,j]
    sliceDF <- t(sliceDF)
    
    iDF <- im(sliceDF, xcol = lon, 
              yrow = lat) 
    
    iDF <- as.data.frame(iDF)
    outDF1 <- subset(iDF, x <= 180)
    outDF2 <- subset(iDF, x > 180)
    outDF2$x <- outDF2$x - 360
    iDF <- rbind(outDF1, outDF2)
    
    DF[,(j-t.start+3)] <- iDF[,3]
  }
  
  colnames(DF) <- c("x","y", c(t.start:nt))
  
  ##Extract for Australia, and for the period of 200101-200512
  outDF <- subset(DF, x <= lon.max & x >= lon.min & y <= lat.max & y >= lat.min)

  
  ##Save data as excel file
  outName <- DatFiles[i]
  outName <- sub("*\\.nc", "", outName)
  outNamef <- substr(outName, 1,(nchar(outName)-14))
  outName <- paste(outNamef, "_200101-200512.csv", sep="")
  
  write.table(outDF, paste(destDir, "/", outName, sep=""),
              row.names=F, col.names=T, sep=",")
}

#########################################################
##Function to extract Australia (200101-200512 based on irregular point method)
##count # of months backwards for data ended by 201012

CMIPslice4 <- function()
{
  inName <- file.path(sourceDir, DatFiles[i], fsep = .Platform$file.sep)
  check <- nc_open(inName)
  
  ##Get lon lat information
  lon <- ncvar_get(check, "lon")
  nlon <- dim(lon)
  lonlist <- lon
  lat <- ncvar_get(check, "lat")
  nlat <- dim(lat)
  latlist <- lat
  
  ##Obtain time coverage information
  t <- ncvar_get(check, "time")
  nt <- dim(t)
  t.start <- nt - 119
  t.end <- nt - 60
  
  ##create 3d matrix
  tmp.array <- ncvar_get(check, "gpp")
  nc_close(check)
  
  ##extract lon lat according to time j
  j <- t.start
  
  #dim(tmp.array)
  sliceDF <- tmp.array[,,j]
  sliceDF <- t(sliceDF)
  
  DF <- im(sliceDF, xcol = lon, 
           yrow = lat) 
  
  DF <- as.data.frame(DF)
  outDF1 <- subset(DF, x <= 180)
  outDF2 <- subset(DF, x > 180)
  outDF2$x <- outDF2$x - 360
  DF <- rbind(outDF1, outDF2)
  
  for (j in (t.start+1):t.end)
  {
    sliceDF <- tmp.array[,,j]
    sliceDF <- t(sliceDF)
    
    iDF <- im(sliceDF, xcol = lon, 
              yrow = lat) 
    
    iDF <- as.data.frame(iDF)
    outDF1 <- subset(iDF, x <= 180)
    outDF2 <- subset(iDF, x > 180)
    outDF2$x <- outDF2$x - 360
    iDF <- rbind(outDF1, outDF2)
    
    DF[,(j-t.start+3)] <- iDF[,3]
  }
  
  colnames(DF) <- c("x","y", c(t.start:t.end))
  
  ##Extract for Australia, and for the period of 200101-200512
  outDF <- subset(DF, x <= lon.max & x >= lon.min & y <= lat.max & y >= lat.min)
  
  
  ##Save data as excel file
  outName <- DatFiles[i]
  outName <- sub("*\\.nc", "", outName)
  outNamef <- substr(outName, 1,(nchar(outName)-14))
  outName <- paste(outNamef, "_200101-200512.csv", sep="")
  
  write.table(outDF, paste(destDir, "/", outName, sep=""),
              row.names=F, col.names=T, sep=",")
}


#########################################################
##Downscale CMIP5 data spatial resolution to modis resolution
consist_reso <- function()
{
  for (i in 1:length(DatFiles))
  {
    inName <- file.path(cmipsourceDir, DatFiles[i], fsep = .Platform$file.sep)
    cmip <- read.table(inName, header=T, sep=",")
    
    for (j in 3:62)
    {
      cmipsub <- cbind(cmip[,1:2], cmip[,j])
      coordinates(cmipsub) <- ~x+y
      gridded(cmipsub) <- T
      r <- raster(cmipsub)
      
      outDF[,j] <- extract(r, corDF)
    }
    
    write.table(outDF, paste(destDir, DatFiles[i], sep="/"),
                col.names=T,row.names=F, sep=",")
    
  }
}

#########################################################
##Project MODIS data onto CMIP individual spaital resolution
modisToCmip <- function()
{
  for (i in 1:length(cmipFiles))
  {
    inName1 <- file.path(cmipDir, cmipFiles[i], fsep = .Platform$file.sep)
    outName <- cmipFiles[i]
    outNamef <- substr(outName, 10,(nchar(outName)-36))
    
    cmipDF <- read.table(inName1, header=T,sep=",")
    cmipDF <- as.data.frame(cmipDF[,1:3])
    colnames(cmipDF) <- c("x","y","value")
    
    coordinates(cmipDF) <- ~x+y
    gridded(cmipDF) <-T
    cmip_r <- raster(cmipDF)
    
    ##project modis onto cmip grids
    j <- 1
    inName2 <- file.path(modisDir, modisFiles[j], fsep = .Platform$file.sep)
    
    d <- raster(inName2, varname = "total")
    out <- resample(d, cmip_r, fun=mean)
    outDF <- rasterToPoints(out)
    outDF <- as.data.frame(outDF, stringsAsFactors=F)
    colnames(outDF) <- c("x","y","GPP")
    
    for (j in 2:length(modisFiles))
    {
      inName2 <- file.path(modisDir, modisFiles[j], fsep = .Platform$file.sep)
      
      ##aggregate data to coarser resolution
      d <- raster(inName2, varname = "total")
      
      out <- resample(d, cmip_r, fun=mean)
      outP <- rasterToPoints(out)
      outP <- as.data.frame(outP, stringsAsFactors=F)
      colnames(outP) <- c("x","y","GPP")
      
      outDF[,(j+2)] <- outP$GPP
    }
    
    write.table(outDF, paste(destDir, "/MODIS_", outNamef, ".csv", sep=""),
                col.names=T,row.names=F, sep=",")
    
  }
  
}
#########################################################
##Compute CMIP vs MODIS spatial comparison
##including absolute diff, % diff
##5-yr mean annual GPP CMIP
##5-yr mean annual GPP modis
##histogram of % diff

compareF <- function()
{
  for (i in 1:length(model))
  {
    cmipName <- paste(cmipsourceDir, "/gpp_Lmon_", model[i], "_historical_r1i1p1_200101-200512.csv",sep="")
    modisName <- paste(modisDir, "/MODIS_", model[i], ".csv",sep="")
    
    cmip <- read.table(cmipName, header=T, sep=",")
    modis <- read.table(modisName, header=T, sep=",")
    
    ##Convert cmip5 gpp unit from kg m-2 s-1 to g m-2 d-1
    cmip[,3:62] <- cmip[,3:62] * 1000 * 3600 * 24
    
    modis_ext <- as.data.frame(cmip[,1:2])
    
    ##Extract Modis based on cmip grid number
    for(j in 3:62)
    {
      corDF <- as.data.frame(cbind(cmip[,1:2]))
      coordinates(corDF) <- ~x+y
      gridded(corDF) <- T

      modisnew <- as.data.frame(cbind(modis[,1:2], modis[,j]))
      coordinates(modisnew) <- ~x+y
      gridded(modisnew) <- T
      modis_r <- raster(modisnew)
      
      out <- extract(modis_r, corDF)
      modis_ext[,j] <- out
    }
    
    ##Compute % diff between cmip and modis, at cmip grids
    outdat <- (cmip[,3:62]-modis_ext[,3:62])/modis_ext[,3:62] * 100
    out <- cbind(cmip[,1:2], outdat)
    
    write.table(out, paste(destDir1, "/gpp_Lmon_", model[i], "_historical_r1i1p1_200101-200512.csv",sep=""),
                col.names=T,row.names=F, sep=",")
    
    ##Compute monthly sum (g m-2 month-1) both modis and cmip
    cmip_sum <- out[,1:2]
    cmip_sum[3:62] <- cmip[,3:62] * 31
    modis_sum <- modis_ext
    modis_sum[3:62] <- modis_ext[,3:62] * 31
    
    m30 <- c(4,6,9,11,
             16,18,21,23,
             28,30,33,35,
             40,42,45,47,
             52,54,57,59)
    
    m28 <- c(2,14,26,50)
    m29 <- 38
    
    for (j in m30)
    {
      cmip_sum[,j+2] <- cmip_sum[,j+2]/31 * 30
      modis_sum[,j+2] <- modis_sum[,j+2]/31 *30
    }
    
    for (j in m28)
    {
      cmip_sum[,j+2] <- cmip_sum[,j+2]/31 * 28
      modis_sum[,j+2] <- modis_sum[,j+2]/31 *28
    }
    
    cmip_sum[,m29+2] <- cmip_sum[,m29+2]/31 * 29
    modis_sum[,m29+2] <- modis_sum[,m29+2]/31 *29
    
    cmip_sum <- as.data.frame(cmip_sum, stringsAsFactors=F)
    modis_sum <- as.data.frame(modis_sum, stringsAsFactors=F)
    
    ##Compute annual sum for both modis and gpp
    cmip_ann <- cmip_sum[,1:2]
    modis_ann <- modis_sum[,1:2]
    
    cmip_ann[,3] <- cmip_sum[,3] + cmip_sum[,4] + cmip_sum[,5] + cmip_sum[,6] + cmip_sum[,7] + cmip_sum[,8] +
      cmip_sum[,9] + cmip_sum[,10] + cmip_sum[,11] + cmip_sum[,12] + cmip_sum[,13] + cmip_sum[,14]
    
    cmip_ann[,4] <- cmip_sum[,15] + cmip_sum[,16] + cmip_sum[,17] + cmip_sum[,18] + cmip_sum[,19] + cmip_sum[,20] +
      cmip_sum[,21] + cmip_sum[,22] + cmip_sum[,23] + cmip_sum[,24] + cmip_sum[,25] + cmip_sum[,26]
    
    cmip_ann[,5] <- cmip_sum[,27] + cmip_sum[,28] + cmip_sum[,29] + cmip_sum[,30] + cmip_sum[,31] + cmip_sum[,32] +
      cmip_sum[,33] + cmip_sum[,34] + cmip_sum[,35] + cmip_sum[,36] + cmip_sum[,37] + cmip_sum[,38]
    
    cmip_ann[,6] <- cmip_sum[,39] + cmip_sum[,40] + cmip_sum[,41] + cmip_sum[,42] + cmip_sum[,43] + cmip_sum[,44] +
      cmip_sum[,45] + cmip_sum[,46] + cmip_sum[,47] + cmip_sum[,48] + cmip_sum[,49] + cmip_sum[,50]
    
    cmip_ann[,7] <- cmip_sum[,51] + cmip_sum[,52] + cmip_sum[,53] + cmip_sum[,54] + cmip_sum[,55] + cmip_sum[,56] +
      cmip_sum[,57] + cmip_sum[,58] + cmip_sum[,59] + cmip_sum[,60] + cmip_sum[,61] + cmip_sum[,62]
    
    cmip_ann <- as.data.frame(cmip_ann, stringsAsFactors=F)
    colnames(cmip_ann) <- c("x","y","y2001","y2002","y2003","y2004","y2005")
    
    modis_ann[,3] <- modis_sum[,3] + modis_sum[,4] + modis_sum[,5] + modis_sum[,6] + modis_sum[,7] + modis_sum[,8] +
      modis_sum[,9] + modis_sum[,10] + modis_sum[,11] + modis_sum[,12] + modis_sum[,13] + modis_sum[,14]
    
    modis_ann[,4] <- modis_sum[,15] + modis_sum[,16] + modis_sum[,17] + modis_sum[,18] + modis_sum[,19] + modis_sum[,20] +
      modis_sum[,21] + modis_sum[,22] + modis_sum[,23] + modis_sum[,24] + modis_sum[,25] + modis_sum[,26]
    
    modis_ann[,5] <- modis_sum[,27] + modis_sum[,28] + modis_sum[,29] + modis_sum[,30] + modis_sum[,31] + modis_sum[,32] +
      modis_sum[,33] + modis_sum[,34] + modis_sum[,35] + modis_sum[,36] + modis_sum[,37] + modis_sum[,38]
    
    modis_ann[,6] <- modis_sum[,39] + modis_sum[,40] + modis_sum[,41] + modis_sum[,42] + modis_sum[,43] + modis_sum[,44] +
      modis_sum[,45] + modis_sum[,46] + modis_sum[,47] + modis_sum[,48] + modis_sum[,49] + modis_sum[,50]
    
    modis_ann[,7] <- modis_sum[,51] + modis_sum[,52] + modis_sum[,53] + modis_sum[,54] + modis_sum[,55] + modis_sum[,56] +
      modis_sum[,57] + modis_sum[,58] + modis_sum[,59] + modis_sum[,60] + modis_sum[,61] + modis_sum[,62]
    
    modis_ann <- as.data.frame(modis_ann, stringsAsFactors=F)
    colnames(modis_ann) <- c("x","y","y2001","y2002","y2003","y2004","y2005")
    
    ##Compute 5-yr mean for modis and gpp
    cmip_ann[,"totmean"] <- (cmip_ann$y2001 + cmip_ann$y2002 + cmip_ann$y2003 +
      cmip_ann$y2004 + cmip_ann$y2005)/5
    
    modis_ann[,"totmean"] <- (modis_ann$y2001 + modis_ann$y2002 + modis_ann$y2003 +
                               modis_ann$y2004 + modis_ann$y2005)/5
    
    write.table(cmip_ann, paste(destDir2, "/gpp_Lmon_", model[i], "_historical_r1i1p1_200101-200512.csv",sep=""),
                col.names=T,row.names=F, sep=",")
    
    write.table(modis_ann, 
                paste("~/Documents/PostDoc/GDAY/MODIS_GPP/extracted/MODIS_", model[i], "_annual_sum.csv", sep=""),
                sep=",", col.names=T, row.names=F)
    
    
    ##Compute absolute diff for 5yr mean annual GPP
    ann_diff <- cmip_ann[,1:2]
    ann_diff[,"ann_abs"] <- cmip_ann[,"totmean"] - modis_ann[,"totmean"]
        
    ##Compute % diff for 5yr mean annual GPP
    ann_percent <- cmip_ann[,1:2]
    ann_percent[,"ann_percent"] <- (cmip_ann[,"totmean"] - modis_ann[,"totmean"]) / modis_ann[,"totmean"] * 100
    
    modis_plot <- cbind(modis_ann[,1:2], modis_ann[,7])
    cmip_plot <- cbind(cmip_ann[,1:2], cmip_ann[,7])
    
    ##Plotting figures
    pdf(paste(plotdest, "/", outNamef, ".pdf", sep=""))
    
    titleM <- expression("annual mean GPP (g C " * m^-2 * yr^-1 * ")")
    coordinates(modis_plot) <- ~x+y
    gridded(modis_plot) <- T
    r <- raster(modis_plot)
    plot(r)
    legend("topleft", "MODIS")
    title(titleM)
    
    titleM <- expression("annual mean GPP (g C " * m^-2 * yr^-1 * ")")
    coordinates(cmip_plot) <- ~x+y
    gridded(cmip_plot) <- T
    r <- raster(cmip_plot)
    plot(r)
    legend("topleft", model[i])
    title(titleM)
    
    titleM <- expression("GPP difference (g C " * m^-2 * yr^-1 * ")")
    coordinates(ann_diff) <- ~x+y
    gridded(ann_diff) <- T
    r <- raster(ann_diff)
    plot(r)
    legend("topleft", paste(model[i], " - MODIS", sep=""))
    title(titleM)
    
    titleM <- expression("GPP percent difference (%)")
    coordinates(ann_percent) <- ~x+y
    gridded(ann_percent) <- T
    r <- raster(ann_percent)
    plot(r)
    legend("topleft", model[i])
    title(titleM)
    
    hist(ann_percent$ann_percent, xlab = "GPP percent difference (%)",
         main = NA, col = "grey",prob = T)
    lines(density(ann_percent$ann_percent), col = "red", lwd = 2)
    
    ##Generating time series gpp of continental mean (averaged across all grids)
    ##For both modis and cmip
    modis_colmean <- colMeans(modis_ext, na.rm = T)
    cmip_colmean <- colMeans(cmip, na.rm=T)
    new <- as.data.frame(cbind(modis_colmean, cmip_colmean))
    new <- new[3:62,]
    
    y.min <- round(min(new) - 1, 0)
    y.max <- round(max(new) + 1, 0)
    y.lab <- expression("Monthly mean daily GPP (g C " * m^-2 * d^-1 * " )")
    
    new$date <- seq(as.Date("2001/01/01"), 
                    by = "month", length.out = 60)
    
    plot(modis_colmean~date, data = new, type = "b", ylim=c(y.min,y.max),
         lwd = 2, col = "black", pch = 16, ylab = y.lab)
    lines(cmip_colmean~date, data = new, type="b", lwd = 2, col = "red",
          pch = 16)
    legend("topright", c("MODIS", model[i]), col = c("black", "red"), pch = 16, lwd = 2)
    
    taylor.diagram(modis_ann$totmean, cmip_ann$totmean, main = "Normalized Taylor Diagram of Annual GPP",
                   pch = 16, xlab = "Standard deviation", normalize=T)
    legend("topright", c("MODIS", model[i]), col = c("black", "red"), pch = c(1,16))
    
    dev.off()

  }
}

#########################################################
##Convert from x, y, month 1, month 2, ... to
##x, y, date, value format
##Apply to gridded monthly gpp data of modis and cmip
##Save onto the same dataframe for taylor diagram and other statistical analyses

DataPrepForTaylor <- function()
{
  datelist <- seq(as.Date("2001/01/01"), by = "month", length.out = 60)
  outDF <- matrix(nrow=60, ncol=length(model)+2)
  outDF <- as.data.frame(outDF)
  colnames(outDF) <- c("date", "modis", model)
  outDF$date <- datelist
  
  for (i in 1:length(model))
  {
    cmipName <- paste(cmipsourceDir, "/gpp_Lmon_", model[i], "_historical_r1i1p1_200101-200512.csv",sep="")
    modisName <- paste(modisDir, "/MODIS_", model[i], ".csv",sep="")
    
    cmip <- read.table(cmipName, header=T, sep=",")
    modis <- read.table(modisName, header=T, sep=",")
    
    ##Convert cmip5 gpp unit from kg m-2 s-1 to g m-2 d-1
    cmip[,3:62] <- cmip[,3:62] * 1000 * 3600 * 24
    
    modis_ext <- as.data.frame(cmip[,1:2])
    
    gn <- length(cmip$x)
    
    ##Extract Modis based on cmip grid number
    for(j in 3:62)
    {
      corDF <- as.data.frame(cbind(cmip[,1:2]))
      coordinates(corDF) <- ~x+y
      gridded(corDF) <- T
      
      modisnew <- as.data.frame(cbind(modis[,1:2], modis[,j]))
      coordinates(modisnew) <- ~x+y
      gridded(modisnew) <- T
      modis_r <- raster(modisnew)
      
      out <- extract(modis_r, corDF)
      modis_ext[,j] <- out
    }
    
    ##Convert from x, y, month 1... to x, y, month, year, value format
    #xlist <- rep(cmip[1:length(cmip$x), 1], times = 60)
    #ylist <- rep(cmip[1:length(cmip$x), 2], times = 60)
    #m <- rep(c(1:12), each = gn, times = 5)
    #y <- rep(c(2001:2005), each = gn*12)
    #outDF <- as.data.frame(cbind(xlist, ylist, y, m))

    #for (j in 1:60)
    #{
    #  s <- 1 + (j-1) * gn
    #  e <- j * gn
    #  outDF[s:e,"modis"] <- modis_ext[1:gn,j+2]
    #  outDF[s:e, "cmip"] <- cmip[1:gn, j+2]
    #}
    
    cmipC <- colMeans(cmip, na.rm=T)
    modisC <- colMeans(modis_ext, na.rm=T)
    
    tempDF <- rbind(cmipC, modisC)
    tempDF <- as.data.frame(t(tempDF[,3:62]))
    colnames(tempDF) <- c("CMIP", "MODIS")
    outDF$modis <- tempDF$MODIS
    outDF[,i+2] <- tempDF$CMIP
    
  }
  
  write.table(outDF, paste(destDir, "/CMIP5_models_time_series_data.csv", sep=""),
              col.names=T, row.names=F, sep=",")
}