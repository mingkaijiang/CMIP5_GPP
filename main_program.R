#########################################################
##CMIP5 output assessment

#########################################################
source("code/Function_CMIP5_V1.2.R")


#########################################################
# note: modis datasets are projected onto cmip resolution
# cmip unit is in kg m-2 s-1, need to convert to g m-2 d-1

pdf("data/plot/Taylor_diamgram_GPP.pdf")

### read in modis (m) and cmip (c) datasets
m1DF <- read.csv("data/modis/historic/MODIS_CMCC-CESM.csv")
c1DF <- read.csv("data/extracted/historic/gpp_Lmon_CMCC-CESM_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m1DF$mean <- rowMeans(m1DF[,3:62])
m1DF <- m1DF[,c("x", "y", "mean")]

c1DF$mean <- rowMeans(c1DF[,3:62])
c1DF <- c1DF[,c("x", "y", "mean")]

# convert names
p1DF <- merge(m1DF, c1DF, by.x=c("x","y"), by.y=c("x","y"))
p1DF$c1 <- p1DF$mean.y * 1000 * 3600 * 24
colnames(p1DF) <- c("x", "y", "m1", "test", "c1")

# plot
p1 <- taylor.diagram(p1DF$m1, p1DF$c1, col="red", pch=19,  pcex=2, normalize=T,
                     xlab="Normalized GPP", main="")

### read in modis (m) and cmip (c) datasets
m2DF <- read.csv("data/modis/historic/MODIS_GFDL-ESM2G.csv")
c2DF <- read.csv("data/extracted/historic/gpp_Lmon_GFDL-ESM2G_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m2DF$mean <- rowMeans(m2DF[,3:62])
m2DF <- m2DF[,c("x", "y", "mean")]

c2DF$mean <- rowMeans(c2DF[,3:62])
c2DF <- c2DF[,c("x", "y", "mean")]

# convert names
p2DF <- merge(m2DF, c2DF, by.x=c("x","y"), by.y=c("x","y"))
p2DF$c2 <- p2DF$mean.y * 1000 * 3600 * 24
colnames(p2DF) <- c("x", "y", "m2", "test", "c2")

# plot
p2 <- taylor.diagram(p2DF$m2, p2DF$c2, add=T, col="blue", pch=19,  pcex=2, normalize=T)


### read in modis (m) and cmip (c) datasets
m3DF <- read.csv("data/modis/historic/MODIS_GFDL-ESM2M.csv")
c3DF <- read.csv("data/extracted/historic/gpp_Lmon_GFDL-ESM2M_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m3DF$mean <- rowMeans(m3DF[,3:62])
m3DF <- m3DF[,c("x", "y", "mean")]

c3DF$mean <- rowMeans(c3DF[,3:62])
c3DF <- c3DF[,c("x", "y", "mean")]

# convert names
p3DF <- merge(m3DF, c3DF, by.x=c("x","y"), by.y=c("x","y"))
p3DF$c3 <- p3DF$mean.y * 1000 * 3600 * 24
colnames(p3DF) <- c("x", "y", "m3", "test", "c3")

# plot
p3 <- taylor.diagram(p3DF$m3, p3DF$c3, add=T, col="orange", pch=19,  pcex=2, normalize=T)



### read in modis (m) and cmip (c) datasets
m3DF <- read.csv("data/modis/historic/MODIS_GISS-E2-H-CC.csv")
c3DF <- read.csv("data/extracted/historic/gpp_Lmon_GISS-E2-H-CC_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m3DF$mean <- rowMeans(m3DF[,3:62])
m3DF <- m3DF[,c("x", "y", "mean")]

c3DF$mean <- rowMeans(c3DF[,3:62])
c3DF <- c3DF[,c("x", "y", "mean")]

# convert names
p3DF <- merge(m3DF, c3DF, by.x=c("x","y"), by.y=c("x","y"))
p3DF$c3 <- p3DF$mean.y * 1000 * 3600 * 24
colnames(p3DF) <- c("x", "y", "m3", "test", "c3")

# plot
p4 <- taylor.diagram(p3DF$m3, p3DF$c3, add=T, col="green", pch=19,  pcex=2, normalize=T)


### read in modis (m) and cmip (c) datasets
m3DF <- read.csv("data/modis/historic/MODIS_GISS-E2-H.csv")
c3DF <- read.csv("data/extracted/historic/gpp_Lmon_GISS-E2-H_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m3DF$mean <- rowMeans(m3DF[,3:62])
m3DF <- m3DF[,c("x", "y", "mean")]

c3DF$mean <- rowMeans(c3DF[,3:62])
c3DF <- c3DF[,c("x", "y", "mean")]

# convert names
p3DF <- merge(m3DF, c3DF, by.x=c("x","y"), by.y=c("x","y"))
p3DF$c3 <- p3DF$mean.y * 1000 * 3600 * 24
colnames(p3DF) <- c("x", "y", "m3", "test", "c3")

# plot
p5 <- taylor.diagram(p3DF$m3, p3DF$c3, add=T, col="purple", pch=19,  pcex=2, normalize=T)

### read in modis (m) and cmip (c) datasets
m3DF <- read.csv("data/modis/historic/MODIS_GISS-E2-R-CC.csv")
c3DF <- read.csv("data/extracted/historic/gpp_Lmon_GISS-E2-R-CC_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m3DF$mean <- rowMeans(m3DF[,3:62])
m3DF <- m3DF[,c("x", "y", "mean")]

c3DF$mean <- rowMeans(c3DF[,3:62])
c3DF <- c3DF[,c("x", "y", "mean")]

# convert names
p3DF <- merge(m3DF, c3DF, by.x=c("x","y"), by.y=c("x","y"))
p3DF$c3 <- p3DF$mean.y * 1000 * 3600 * 24
colnames(p3DF) <- c("x", "y", "m3", "test", "c3")

# plot
p6 <- taylor.diagram(p3DF$m3, p3DF$c3, add=T, col="brown", pch=19,  pcex=2, normalize=T)

### read in modis (m) and cmip (c) datasets
m3DF <- read.csv("data/modis/historic/MODIS_GISS-E2-R.csv")
c3DF <- read.csv("data/extracted/historic/gpp_Lmon_GISS-E2-R_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m3DF$mean <- rowMeans(m3DF[,3:62])
m3DF <- m3DF[,c("x", "y", "mean")]

c3DF$mean <- rowMeans(c3DF[,3:62])
c3DF <- c3DF[,c("x", "y", "mean")]

# convert names
p3DF <- merge(m3DF, c3DF, by.x=c("x","y"), by.y=c("x","y"))
p3DF$c3 <- p3DF$mean.y * 1000 * 3600 * 24
colnames(p3DF) <- c("x", "y", "m3", "test", "c3")

# plot
p7 <- taylor.diagram(p3DF$m3, p3DF$c3, add=T, col="cyan", pch=19,  pcex=2, normalize=T)

### read in modis (m) and cmip (c) datasets
m3DF <- read.csv("data/modis/historic/MODIS_HadCM3.csv")
c3DF <- read.csv("data/extracted/historic/gpp_Lmon_HadCM3_historical_r1i1p1_200101-200512.csv")

# compute 5-year mean daily GPP for modis and cmip
m3DF$mean <- rowMeans(m3DF[,3:62])
m3DF <- m3DF[,c("x", "y", "mean")]

c3DF$mean <- rowMeans(c3DF[,3:62])
c3DF <- c3DF[,c("x", "y", "mean")]

# convert names
p3DF <- merge(m3DF, c3DF, by.x=c("x","y"), by.y=c("x","y"))
p3DF$c3 <- p3DF$mean.y * 1000 * 3600 * 24
colnames(p3DF) <- c("x", "y", "m3", "test", "c3")

# plot
p8 <- taylor.diagram(p3DF$m3, p3DF$c3, add=T, col="pink", pch=19,  pcex=2, normalize=T)

dev.off()