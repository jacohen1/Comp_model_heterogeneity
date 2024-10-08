####initial code####
rm(list=ls())

#set working directory
setwd("")

#read in data
rd <- read.csv("RD_results.csv", header = TRUE)
dd <- read.csv("DD_results.csv", header = TRUE)

#create values for i and j
i <- 1:nrow(rd)
j <- 1:nrow(dd)

####R0 calculation for RD####

#create function for the R0 equation for RD
r0.func.rd <- function(i){
  k <- 0.5
  S1 <- 49
  S2 <- 49
  rho11 <- 1
  rho12 <- (1-rho11)
  rho22 <- 1
  rho21 <- (1-rho22)
  d <- 1
  alpha <- 0.7
  sigma1 <- rd$sigma1[[i]]
  sigma2 <- rd$sigma2[[i]]
  iota1 <- rd$iota1[[i]]
  iota2 <- rd$iota2[[i]]
  R0.rd <- (k*((S1*sigma1*((iota1*rho11)+(iota2*rho12)))+(S2*sigma2*((iota1*rho21)+(iota2*rho22)))))/(d+alpha)
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
  gamma11 <- 1
  gamma12 <- (1-gamma11)
  gamma22 <- 1
  gamma21 <- (1-gamma22)
  d <- 1
  alpha <- 0.7
  sigma1 <- dd$sigma1[[j]]
  sigma2 <- dd$sigma2[[j]]
  iota1 <- dd$iota1[[j]]
  iota2 <- dd$iota2[[j]]
  R0.dd <- ((k*((gamma11*iota1)+(gamma22*iota2))+
               sqrt(((gamma11^2)*(iota1^2))+(4*gamma12*gamma21*iota1*iota2)-
                      (2*gamma11*gamma22*iota1*iota2)+((gamma22^2)*(iota2^2))))*
              ((S1*sigma1)+(S2*sigma2)))/(2*(d+alpha))
  return(R0.dd)
}

#get R0 value for each simulation of DD
r0.dd <- lapply(j, r0.func.dd)

#add R0 values to dd data frame
dd$R0 <- unlist(r0.dd)

#remove function and list
rm(r0.func.dd, r0.dd)

####Heterogeneity calculation for RD####
#create function for getting starting weighted centroid for sigma
sig.cent.start <- function(i){
  c.sig <- ((rd$sigma1[[i]] * (49 + 1)) 
            + (rd$sigma2[[i]] * (49 + 1))
  )/(49 + 1 + 49 + 1)
  return(c.sig)
}

#create function for getting starting weighted centroid for iota
iot.cent.start <- function(i){
  c.iot <- ((rd$iota1[[i]] * (49 + 1)) 
            + (rd$iota2[[i]] * (49 + 1))
  )/(49 + 1 + 49 + 1)
  return(c.iot)
}

#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.start <- function(i){
  dis1 <- sqrt((rd$sigma1[[i]] - sig.cent.start(i))^2 + (rd$iota1[[i]] - iot.cent.start(i))^2)
  dis2 <- sqrt((rd$sigma2[[i]] - sig.cent.start(i))^2 + (rd$iota2[[i]] - iot.cent.start(i))^2)
  het.start <- ((dis1 * (49 + 1)) + (dis2 * (49 + 1)))/(49 + 1 + 49 + 1)
  return(het.start)
}

#run the functions with the data
het.s <- lapply(i, dists.start)

#put values in the original data frame
rd$het.start <- unlist(het.s)

#remove functions
rm(sig.cent.start, iot.cent.start, dists.start, het.s)

##create function for getting ending weighted centroid for sigma
sig.cent.end <- function(i){
  c.sig <- ((rd$sigma1[[i]] * (rd$S1[[i]] + rd$I1[[i]])) 
            + (rd$sigma2[[i]] * (rd$S2[[i]] + rd$I2[[i]]))
  )/(rd$S1[[i]] + rd$I1[[i]] + rd$S2[[i]] + rd$I2[[i]])
  return(c.sig)
}

##create function for getting ending weighted centroid for iota
iot.cent.end <- function(i){
  c.iot <- ((rd$iota1[[i]] * (rd$S1[[i]] + rd$I1[[i]])) 
            + (rd$iota2[[i]] * (rd$S2[[i]] + rd$I2[[i]]))
  )/(rd$S1[[i]] + rd$I1[[i]] + rd$S2[[i]] + rd$I2[[i]])
  return(c.iot)
}

#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.end <- function(i){
  dis1 <- sqrt((rd$sigma1[[i]] - sig.cent.end(i))^2 + (rd$iota1[[i]] - iot.cent.end(i))^2)
  dis2 <- sqrt((rd$sigma2[[i]] - sig.cent.end(i))^2 + (rd$iota2[[i]] - iot.cent.end(i))^2)
  het.end <- ((dis1 * (rd$S1[[i]] + rd$I1[[i]])) + (dis2 * (rd$S2[[i]] + rd$I2[[i]])))/(rd$S1[[i]] + rd$I1[[i]] + rd$S2[[i]] + rd$I2[[i]])
  return(het.end)
}

#run the functions with the data
het.e <- lapply(i, dists.end)

#put values in the original data frame
rd$het.end <- unlist(het.e)

#remove functions
rm(sig.cent.end, iot.cent.end, dists.end, het.e)

####Heterogeneity calculation for DD####
#create function for getting starting weighted centroid for sigma
sig.cent.start <- function(j){
  c.sig <- ((dd$sigma1[[j]] * (49 + 1)) 
            + (dd$sigma2[[j]] * (49 +1))
  )/(49 + 1 + 49 + 1)
  return(c.sig)
}

#create function for getting starting weighted centroid for iota
iot.cent.start <- function(j){
  c.iot <- ((dd$iota1[[j]] * (1)) 
            + (dd$iota2[[j]] * (1))
  )/(1 + 1)
  return(c.iot)
}

#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.start <- function(j){
  disI11 <- sqrt((dd$sigma1[[j]] - sig.cent.start(j))^2 + (dd$iota1[[j]] - iot.cent.start(j))^2)
  disI21 <- sqrt((dd$sigma2[[j]] - sig.cent.start(j))^2 + (dd$iota1[[j]] - iot.cent.start(j))^2)
  disI12 <- sqrt((dd$sigma1[[j]] - sig.cent.start(j))^2 + (dd$iota2[[j]] - iot.cent.start(j))^2)
  disI22 <- sqrt((dd$sigma2[[j]] - sig.cent.start(j))^2 + (dd$iota2[[j]] - iot.cent.start(j))^2)
  disS1 <- abs(dd$sigma1[[j]] - sig.cent.start(j))
  disS2 <- abs(dd$sigma2[[j]] - sig.cent.start(j))
  het.start <- ((disI11 * (0.5)) + (disI21 * (0.5)) + (disI12 * (0.5)) + (disI22 * (0.5)) + (disS1 * (49)) + (disS2 * (49)))/(49 + 1 + 1 + 49)
  return(het.start)
}

#run the functions with the data
het.s <- lapply(j, dists.start)

#put values in the original data frame
dd$het.start <- unlist(het.s)

#remove functions
rm(sig.cent.start, iot.cent.start, dists.start, het.s)

##create function for getting ending weighted centroid for sigma
sig.cent.end <- function(j){
  c.sig <- ((dd$sigma1[[j]] * (dd$S1[[j]] + dd$I11[[j]] + dd$I12[[j]])) 
            + (dd$sigma2[[j]] * (dd$S2[[j]] + dd$I21[[j]] + dd$I22[[j]]))
  )/(dd$S1[[j]] + dd$I11[[j]] + dd$I12[[j]] + dd$S2[[j]] + dd$I21[[j]] + dd$I22[[j]])
  return(c.sig)
}

##create function for getting ending weighted centroid for iota
iot.cent.end <- function(j){
  c.iot <- ((dd$iota1[[j]] * (dd$I11[[j]] + dd$I21[[j]])) 
            + (dd$iota2[[j]] * (dd$I12[[j]] + dd$I22[[j]]))
  )/(dd$I11[[j]] + dd$I21[[j]] + dd$I12[[j]] + dd$I22[[j]])
  return(c.iot)
}

#create function to calculate euclidean distances to weighted centroid and then calculate
#average weighted distance to the weighted centroid
dists.end <- function(j){
  disI11 <- sqrt((dd$sigma1[[j]] - sig.cent.end(j))^2 + (dd$iota1[[j]] - iot.cent.end(j))^2)
  disI21 <- sqrt((dd$sigma2[[j]] - sig.cent.end(j))^2 + (dd$iota1[[j]] - iot.cent.end(j))^2)
  disI12 <- sqrt((dd$sigma1[[j]] - sig.cent.end(j))^2 + (dd$iota2[[j]] - iot.cent.end(j))^2)
  disI22 <- sqrt((dd$sigma2[[j]] - sig.cent.end(j))^2 + (dd$iota2[[j]] - iot.cent.end(j))^2)
  disS1 <- abs(dd$sigma1[[j]] - sig.cent.end(j))
  disS2 <- abs(dd$sigma2[[j]] - sig.cent.end(j))
  het.end <- ((disI11 * dd$I11[[j]]) + (disI21 * dd$I21[[j]]) + (disI12 * dd$I12[[j]]) + (disI22 * dd$I22[[j]]) + (disS1 * dd$S1[[j]]) + (disS2 * dd$S2[[j]]))/(dd$I11[[j]] + dd$I21[[j]] + dd$I12[[j]] + dd$I22[[j]] + dd$S1[[j]] + dd$S2[[j]])
  return(het.end)
}

#run the functions with the data
het.e <- lapply(j, dists.end)

#put values in the original data frame
dd$het.end <- unlist(het.e)

#remove functions
rm(sig.cent.end, iot.cent.end, dists.end, het.e)

mycov <- function(i){
  correlation <- cov(
    c(rd$sigma1[i], rd$sigma2[i]),
    c(rd$iota1[i], rd$iota2[i]))
  return(correlation)
}

i <- 1:nrow(rd)
cov.res <- lapply(i, mycov)

rd$cov <- unlist(cov.res)

rm(cov.res, mycov)

####Find maximum iota value for simulation DD####
max.iot <- function(j){
  #if iota1 does not equal iota2, get maximum value
  if(dd$iota1[[j]] != dd$iota2[[j]]){
    maximum <- max(dd$iota1[[j]], dd$iota2[[j]])
  }
  #if they do equal each other, make the maximum value 0.5
  else{
    maximum <- 0.5
  }
  #return maximum value
  return(maximum)
}

#find the max
max.dd <- lapply(j, max.iot)
#add to the dataframe
dd$max.iota <- unlist(max.dd)

rm(max.iot, max.dd)

####Find euclidean distance of top point to sigma = 1, iota = 1 for simulation RD####
euc.dist.rd <- function(i){
  #get Equclidean distances for two groups
  dis1 <- sqrt((rd$sigma1[[i]] - 1)^2 + (rd$iota1[[i]] - 1)^2)
  dis2 <- sqrt((rd$sigma2[[i]] - 1)^2 + (rd$iota2[[i]] - 1)^2)
  
  #find which distance is shorter and get the corresponding sigma and iota value
  dis <- min(dis1, dis2)
  if(dis == dis1){
    result <- list(dis, rd$sigma1[[i]], rd$iota1[[i]])
  }
  if(dis == dis2){
    result <- list(dis, rd$sigma2[[i]], rd$iota2[[i]])
  }
  
  #return the smaller of the two distances
  return(result)
}

#find the minimum distance for each simulation
euc.dist <- lapply(i, euc.dist.rd)
#separate the three values into separate lists
min.dist <- sapply(euc.dist, "[[", 1)
sigma.min.dist <- sapply(euc.dist, "[[", 2)
iota.min.dist <- sapply(euc.dist, "[[", 3)
#add to the dataframe
rd$min.dist <- min.dist
rd$sigma.min.dist <- sigma.min.dist
rd$iota.min.dist <- iota.min.dist

rm(euc.dist.rd, euc.dist, min.dist, sigma.min.dist, iota.min.dist)

####save files as csv####
write.csv(dd, "DD_centroid_r0_het.csv", row.names = FALSE)
write.csv(rd, "RD_centroid_r0_het.csv", row.names = FALSE)
