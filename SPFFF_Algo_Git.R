####Approaches for Fast flexible filling split plot designs
library(cluster)
library(acebayes)
library(DiceKriging)
library(lhs)
library(AlgDesign)
library(minimaxdesign)
library(keras)
library(e1071)
library(MOLHD)
library(rpart)
library(flexclust)


SplitPlotFFFDesign <- function(dWP, nWP, dSP, noverall, nsim = noverall*nWP, printsteps = 100, seed = NULL){
  if(is.null(seed) == FALSE){set.seed(seed)}
  X <- matrix(runif(nsim * (dWP + dSP)), ncol= dWP + dSP)
  ncl <- nsim#
  Clusterassignment <- 1:nrow(X)
  Tab <- length(table(Clusterassignment))
# a first run without loop:
  Xsplit <- split(X, Clusterassignment)
  Xsplit <- lapply(Xsplit, helpfct1, d = dWP + dSP)
  ClusterNos <- as.numeric(names(Xsplit))
  LocalMeans <- t(sapply(Xsplit, colMeans))
  LocalN <- unlist(lapply(Xsplit, nrow))
  nCl <- length(LocalN)
  Divisor <- matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow = FALSE) + matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow= TRUE)
  DistLocalM <- as.matrix(dist(LocalMeans))^2
  DistX <- DistLocalM/Divisor
  diag(DistX) <- Inf
  MIN <- min(DistX)
  MIN.ind <- sort(which(DistX == MIN,arr.ind = TRUE)[1,])
  ### here I need to get the indices of above
  Clusterassignment[Clusterassignment == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
  Tab <- Tab - 1
### now iterate by updating DistX and LocalN  
  
  while(Tab >= noverall){
    NewMean <- LocalMeans[MIN.ind[2],]*LocalN[MIN.ind[2]] / sum(LocalN[MIN.ind]) + LocalMeans[MIN.ind[1],] * LocalN[MIN.ind[1]] / sum(LocalN[MIN.ind])
    LocalN[MIN.ind[1]] <- sum(LocalN[MIN.ind])
    RowNamesLocalMeans <- as.numeric(row.names(LocalMeans))
    LocalMeans[MIN.ind[1],] <- NewMean
    LocalN <- LocalN[-MIN.ind[2]]  
    LocalMeans <- LocalMeans[-MIN.ind[2],]
    Disttemp <- (dist2(matrix(NewMean, ncol = dSP + dWP), LocalMeans) / (1/LocalN + 1/LocalN[MIN.ind[1]]))^2
    DistX <- DistX[-MIN.ind[2],-MIN.ind[2]]
    DistX[MIN.ind[1],] <- Disttemp
    DistX[,MIN.ind[1]] <- Disttemp
    DistX[MIN.ind[1],MIN.ind[1]] <- Inf
    MIN <- min(DistX)
    MIN.ind <- sort(which(DistX == MIN,arr.ind = TRUE)[1,])
    Tab <- Tab - 1
    Clusterassignment[RowNamesLocalMeans[MIN.ind[2]]] <- Clusterassignment[RowNamesLocalMeans[MIN.ind[1]]]
    if(Tab%%printsteps == 0){ 
      print(c("Stage 1, " , Tab))
    }
  }

  SPMeans <- LocalMeans
  ClusterassignmentWP <- 1:noverall  
  Tab <- length(table(ClusterassignmentWP))
  Candidates2 <- SPMeans[,1:dWP]
  while(Tab > nWP){
    Xsplit <- split(Candidates2, ClusterassignmentWP)
    Xsplit <- lapply(Xsplit, helpfct1, d = dWP)
    ClusterNos <- as.numeric(names(Xsplit))
    LocalMeans <- t(sapply(Xsplit, colMeans))
    if(dWP == 1){
      LocalMeans <- t(LocalMeans)
    }
    LocalN <- unlist(lapply(Xsplit, nrow))
    nCl <- length(LocalN)
    Divisor <- matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow= FALSE) +   matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow= TRUE)
    DistLocalM <- as.matrix(dist(LocalMeans))^2
    DistX <- DistLocalM/Divisor
    diag(DistX) <- Inf
    MIN <- min(DistX)
    MIN.ind <- which(DistX == MIN,arr.ind = TRUE)[1,]
    ### here I need to get the indices of above
    ClusterNos[MIN.ind]
    ClusterassignmentWP[ClusterassignmentWP == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    Tab <- length(table(ClusterassignmentWP))
    if(Tab%%printsteps == 0){ print(c("Stage 2, " , Tab))}
  }
  ### Create final Design  
  XsplitWP <- split(Candidates2, ClusterassignmentWP)
  XsplitWP <- lapply(XsplitWP, helpfct1, d = dWP)
  WPMeans <- t(sapply(XsplitWP, colMeans))
  if(dWP == 1){
    WPMeans <- t(WPMeans)
  }
  SPMeansBU <- SPMeans  
  for(i in 1:noverall){
    #    i <- 1
    SPMeans[i, 1:dWP] <- WPMeans[ClusterassignmentWP[i] == row.names(WPMeans),]
  }
  DoE <- data.frame(ClusterassignmentWP, SPMeans)
  DoE = DoE[sort(ClusterassignmentWP, index.return = TRUE)$ix,]
  row.names(DoE) = 1:noverall
  DoE <- DoE[sort(DoE[,1], index.return = TRUE)$ix,]
  DoE[,2:(dWP + dSP + 1)] <- DoE[,2:(dWP + dSP + 1)] * 2 - 1
  DoE[,1] <- c(1, cumsum((DoE[-1,1] - DoE[-n,1]) > 0) + 1)
  return(DoE)
}





