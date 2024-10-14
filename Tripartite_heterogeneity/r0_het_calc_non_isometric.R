####initial code####
rm(list=ls())

#set working directory
setwd("")

####Code for first time with new simulation data####
# #read in data to add header and remove one column
# rd <- read.csv("RD_non_isometric.csv", header = FALSE)
# dd <- read.csv("DD_non_isometric.csv", header = FALSE)
#  
# colnames(dd) <- c("Row", "sigma1", "iota1", "sigma2", "iota2", "sigma3", "iota3", "time",
#                            "S1", "S2", "S3", "I1", "I2", "I3",
#                            "I11", "I21", "I31", "I12", "I22", "I32", "I13", "I23", "I33")
# colnames(rd) <- c("Row", "sigma1", "iota1", "sigma2", "iota2", "sigma3", "iota3", "time",
#                            "S1", "S2", "S3", "I1", "I2", "I3")
#  
# #remove "Row" column
# rd <- subset(rd, select = -Row)
# dd <- subset(dd, select = -Row)
# # 
# write.csv(rd, file = "RD_non_isometric.csv", row.names = FALSE)
# write.csv(dd, file = "DD_non_isometric.csv", row.names = FALSE)

####Read in data####
rd <- read.csv("RD_non_isometric.csv", header = TRUE)
dd <- read.csv("DD_non_isometric.csv", header = TRUE)

#create values for i and j
i <- 1:nrow(rd)
j <- 1:nrow(dd)

####R0 calculation for RD####

#create function for the R0 equation for RD
r0.func.rd <- function(i){
  k <- 0.5
  S1 <- 49
  S2 <- 49
  S3 <- 49
  d <- 1
  a <- 0.7
  sigma1 <- rd$sigma1[[i]]
  sigma2 <- rd$sigma2[[i]]
  sigma3 <- rd$sigma3[[i]]
  iota1 <- rd$iota1[[i]]
  iota2 <- rd$iota2[[i]]
  iota3 <- rd$iota3[[i]]
  R0.rd <- (k*((S1*iota1*sigma1) + (S2*iota2*sigma2) + (S3*iota3*sigma3)))/(d+a)
  return(R0.rd)
}

#get R0 value for each simulation of RD
r0.rd <- lapply(i, r0.func.rd)

#ard R0 values to rd data frame
rd$R0 <- unlist(r0.rd)

#remove function and list
rm(r0.func.rd, r0.rd)

####R0 calculation for DD####

#create function for the R0 equation for DD
r0.func.dd <- function(j){
  k <- 0.5
  S1 <- 49
  S2 <- 49
  S3 <- 49
  d <- 1
  a <- 0.7
  sigma1 <- dd$sigma1[[j]]
  sigma2 <- dd$sigma2[[j]]
  sigma3 <- dd$sigma3[[j]]
  iota <- max(dd$iota1[[j]], dd$iota2[[j]], dd$iota3[[j]])
  R0.dd <- (k*iota*((S1*sigma1) + (S2*sigma2) + (S3*sigma3)))/(d+a)
  return(R0.dd)
}

#get R0 value for each simulation of DD
r0.dd <- lapply(j, r0.func.dd)

#add R0 values to dd data frame
dd$R0 <- unlist(r0.dd)

#remove function and list
rm(r0.func.dd, r0.dd)

####Calculate heterogeneity for starting RD####
#create function to calculate euclidean distances to weighted centroid (0.5, 05) and then calculate
#average weighted distance to the weighted centroid
dists.start <- function(i){
  dis1 <- sqrt((rd$sigma1[[i]] - 0.5)^2 + (rd$iota1[[i]] - 0.5)^2)
  dis2 <- sqrt((rd$sigma2[[i]] - 0.5)^2 + (rd$iota2[[i]] - 0.5)^2)
  dis3 <- sqrt((rd$sigma3[[i]] - 0.5)^2 + (rd$iota3[[i]] - 0.5)^2)
  het.start <- ((dis1 * (49 + 1)) + (dis2 * (49 + 1)) + (dis3 * (49 + 1)))/(49 + 1 + 49 + 1 + 49 + 1)
  return(het.start)
}

#run the functions with the data
het.s <- lapply(i, dists.start)

#put values in the original data frame
rd$het.start <- unlist(het.s)

#remove functions
rm(dists.start, het.s)

####Calculate heterogeneity for ending RD####
##create function for getting ending weighted centroid for sigma
sig.cent.end <- function(i){
  c.sig <- ((rd$sigma1[[i]] * (rd$S1[[i]] + rd$I1[[i]])) 
            + (rd$sigma2[[i]] * (rd$S2[[i]] + rd$I2[[i]]))
            + (rd$sigma3[[i]] * (rd$S3[[i]] + rd$I3[[i]]))
  )/(rd$S1[[i]] + rd$I1[[i]] + rd$S2[[i]] + rd$I2[[i]] + rd$S3[[i]] + rd$I3[[i]])
  return(c.sig)
}

##create function for getting ending weighted centroid for iota
iot.cent.end <- function(i){
  c.iot <- ((rd$iota1[[i]] * (rd$S1[[i]] + rd$I1[[i]])) 
            + (rd$iota2[[i]] * (rd$S2[[i]] + rd$I2[[i]]))
            + (rd$iota3[[i]] * (rd$S3[[i]] + rd$I3[[i]]))
  )/(rd$S1[[i]] + rd$I1[[i]] + rd$S2[[i]] + rd$I2[[i]] + rd$S3[[i]] + rd$I3[[i]])
  return(c.iot)
}

#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.end <- function(i){
  dis1 <- sqrt((rd$sigma1[[i]] - sig.cent.end(i))^2 + (rd$iota1[[i]] - iot.cent.end(i))^2)
  dis2 <- sqrt((rd$sigma2[[i]] - sig.cent.end(i))^2 + (rd$iota2[[i]] - iot.cent.end(i))^2)
  dis3 <- sqrt((rd$sigma3[[i]] - sig.cent.end(i))^2 + (rd$iota3[[i]] - iot.cent.end(i))^2)
  het.end <- ((dis1 * (rd$S1[[i]] + rd$I1[[i]])) +
                (dis2 * (rd$S2[[i]] + rd$I2[[i]])) +
                (dis3 * (rd$S3[[i]] + rd$I3[[i]])))/(rd$S1[[i]] + rd$I1[[i]] + rd$S2[[i]] + rd$I2[[i]] + rd$S3[[i]] + rd$I3[[i]])
  return(het.end)
}

#run the functions with the data
het.e <- lapply(i, dists.end)

#put values in the original data frame
rd$het.end <- unlist(het.e)

#remove functions
rm(sig.cent.end, iot.cent.end, dists.end, het.e)

####Calculate heterogeneity for starting DD####
#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.start <- function(j){
  disI11 <- sqrt((dd$sigma1[[j]] - 0.5)^2 + (dd$iota1[[j]] - 0.5)^2)
  disI21 <- sqrt((dd$sigma2[[j]] - 0.5)^2 + (dd$iota1[[j]] - 0.5)^2)
  disI31 <- sqrt((dd$sigma3[[j]] - 0.5)^2 + (dd$iota1[[j]] - 0.5)^2)
  disI12 <- sqrt((dd$sigma1[[j]] - 0.5)^2 + (dd$iota2[[j]] - 0.5)^2)
  disI22 <- sqrt((dd$sigma2[[j]] - 0.5)^2 + (dd$iota2[[j]] - 0.5)^2)
  disI32 <- sqrt((dd$sigma3[[j]] - 0.5)^2 + (dd$iota2[[j]] - 0.5)^2)
  disI13 <- sqrt((dd$sigma1[[j]] - 0.5)^2 + (dd$iota3[[j]] - 0.5)^2)
  disI23 <- sqrt((dd$sigma2[[j]] - 0.5)^2 + (dd$iota3[[j]] - 0.5)^2)
  disI33 <- sqrt((dd$sigma3[[j]] - 0.5)^2 + (dd$iota3[[j]] - 0.5)^2)
  disS1 <- abs(dd$sigma1[[j]] - 0.5)
  disS2 <- abs(dd$sigma2[[j]] - 0.5)
  disS3 <- abs(dd$sigma3[[j]] - 0.5)
  het.start <- ((disI11 * (1/3)) + (disI21 * (1/3)) + (disI31 * (1/3)) + (disS1 * 49) +  
                  (disI12 * (1/3)) + (disI22 * (1/3)) + (disI32 * (1/3)) + (disS2 * 49) +
                  (disI13 * (1/3)) + (disI23 * (1/3)) + (disI33 * (1/3)) + (disS3 * 49))/(49 + 49 + 49 + 1 + 1 + 1)
  return(het.start)
}

#run the functions with the data
het.s <- lapply(j, dists.start)

#put values in the original data frame
dd$het.start <- unlist(het.s)

#remove functions
rm(dists.start, het.s)

####Calculate heterogeneity for ending DD####
##create function for getting ending weighted centroid for sigma
sig.cent.end <- function(j){
  c.sig <- ((dd$sigma1[[j]] * (dd$S1[[j]] + dd$I11[[j]] + dd$I12[[j]] + dd$I13[[j]])) 
            + (dd$sigma2[[j]] * (dd$S2[[j]] + dd$I21[[j]] + dd$I22[[j]] + dd$I23[[j]]))
            + (dd$sigma3[[j]] * (dd$S3[[j]] + dd$I31[[j]] + dd$I32[[j]] + dd$I33[[j]]))
  )/(dd$S1[[j]] + dd$I11[[j]] + dd$I12[[j]] + dd$I13[[j]]
     + dd$S2[[j]] + dd$I21[[j]] + dd$I22[[j]] + dd$I23[[j]]
     + dd$S3[[j]] + dd$I31[[j]] + dd$I32[[j]] + dd$I33[[j]])
  return(c.sig)
}

##create function for getting ending weighted centroid for iota
iot.cent.end <- function(j){
  c.iot <- ((dd$iota1[[j]] * (dd$I11[[j]] + dd$I21[[j]] + dd$I31[[j]])) 
            + (dd$iota2[[j]] * (dd$I12[[j]] + dd$I22[[j]] + dd$I32[[j]]))
            + (dd$iota3[[j]] * (dd$I13[[j]] + dd$I23[[j]] + dd$I33[[j]]))
  )/(dd$I11[[j]] + dd$I21[[j]] + dd$I31[[j]]
     + dd$I12[[j]] + dd$I22[[j]] + dd$I32[[j]]
     + dd$I13[[j]] + dd$I23[[j]] + dd$I33[[j]])
  return(c.iot)
}

#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.end <- function(j){
  disI11 <- sqrt((dd$sigma1[[j]] - sig.cent.end(j))^2 + (dd$iota1[[j]] - iot.cent.end(j))^2)
  disI21 <- sqrt((dd$sigma2[[j]] - sig.cent.end(j))^2 + (dd$iota1[[j]] - iot.cent.end(j))^2)
  disI31 <- sqrt((dd$sigma3[[j]] - sig.cent.end(j))^2 + (dd$iota1[[j]] - iot.cent.end(j))^2)
  disI12 <- sqrt((dd$sigma1[[j]] - sig.cent.end(j))^2 + (dd$iota2[[j]] - iot.cent.end(j))^2)
  disI22 <- sqrt((dd$sigma2[[j]] - sig.cent.end(j))^2 + (dd$iota2[[j]] - iot.cent.end(j))^2)
  disI32 <- sqrt((dd$sigma3[[j]] - sig.cent.end(j))^2 + (dd$iota2[[j]] - iot.cent.end(j))^2)
  disI13 <- sqrt((dd$sigma1[[j]] - sig.cent.end(j))^2 + (dd$iota3[[j]] - iot.cent.end(j))^2)
  disI23 <- sqrt((dd$sigma2[[j]] - sig.cent.end(j))^2 + (dd$iota3[[j]] - iot.cent.end(j))^2)
  disI33 <- sqrt((dd$sigma3[[j]] - sig.cent.end(j))^2 + (dd$iota3[[j]] - iot.cent.end(j))^2)
  disS1 <- abs(dd$sigma1[[j]] - sig.cent.end(j))
  disS2 <- abs(dd$sigma2[[j]] - sig.cent.end(j))
  disS3 <- abs(dd$sigma3[[j]] - sig.cent.end(j))
  het.end <- ((disI11 * dd$I11[[j]]) + (disI21 * dd$I21[[j]]) + (disI31 * dd$I31[[j]]) + (disS1 * dd$S1[[j]])
              + (disI12 * dd$I12[[j]]) + (disI22 * dd$I22[[j]]) + (disI32 * dd$I32[[j]]) + (disS2 * dd$S2[[j]])
              + (disI13 * dd$I13[[j]]) + (disI23 * dd$I23[[j]]) + (disI33 * dd$I33[[j]]) + (disS3 * dd$S3[[j]])
  )/(dd$I11[[j]] + dd$I21[[j]] + dd$I31[[j]] + dd$S1[[j]]
     + dd$I12[[j]] + dd$I22[[j]] + dd$I32[[j]] + dd$S2[[j]]
     + dd$I13[[j]] + dd$I23[[j]] + dd$I33[[j]] + dd$S3[[j]])
  return(het.end)
}

#run the functions with the data
het.e <- lapply(j, dists.end)

#put values in the original data frame
dd$het.end <- unlist(het.e)

#remove functions
rm(sig.cent.end, iot.cent.end, dists.end, het.e)


####Save files as csv####
write.csv(rd, file = "RD_het_R0_non_isometric.csv", row.names = FALSE)
write.csv(dd, file = "DD_het_R0_non_isometric.csv", row.names = FALSE)
