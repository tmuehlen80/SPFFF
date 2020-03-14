####Approaches for Fast flexible filling split plot designs
library(cluster)
# library(mvtnorm)
library(acebayes)
library(DiceKriging)
library(lhs)
library(AlgDesign)
#library(ICAOD)
library(minimaxdesign)
library(keras)
library(e1071)
library(MOLHD)
library(rpart)
library(flexclust)
#citation(package = "e1071")
#citation(package = "base")
#citation(package = "DiceKriging")
#citation(package = "keras")
#citation(package = "lhs")
#citation(package = "acebayes")
#ModelOutputSPFFF1[[i]] <- ModelPreds(
# DoE<- SPFFF[[DoEindicesSPFFF[i]]][,-1]
# WP <- SPFFF[[DoEindicesSPFFF[i]]][,1]
# YReg <- YsSPFFF1[[i]][,2]
# YClass <- YsSPFFF1[[i]][,1]
# DoEVal <- DoEVal
# YRegVal <- YRegVal1
# YClassVal <- YClassVal1
# nepochs <- 10
# nugget <- 0.0001
# costSVMReg <- 1000
# gammaSVMReg <- 0.0001
# costSVMClass <- 100
# gammaSVMclass <- 1



ModelPreds <- function(DoE, WP, YReg, YClass, DoEVal, YRegVal, YClassVal, nepochs = 2, nugget = 0.00001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1){  
  DataReg <- data.frame(Y = YReg, DoE)
  DataClass <- data.frame(Y=as.factor(YClass), DoE)
  LM <- lm(Y ~ X1 + X2 + X3 + X4 + X1*X2 + X1*X3 + X1*X4 + X2*X3 + X2*X4 + X3*X4 + I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2), data=DataReg)
  PredValRegLM <- predict(LM, newdata = DoEVal)
  #  str(PredValLM)
  #plot(Y1.Val[,2], PredVal1.1.2)
  Uniques <- duplicated(DoE) == FALSE
  KM <- km(~1, design = DoE[Uniques,], response = YReg[Uniques], nugget = nugget)
  PredValRegKM <- predict( KM, newdata = DoEVal, type = "UK")$mean
  svm.modelReg <- svm( Y ~ ., data = DataReg, cost = costSVMReg, gamma = gammaSVMReg)
  PredValRegSVM <- predict(svm.modelReg, DoEVal)
  
  KerasModelReg <- keras_model_sequential()
  KerasModelReg %>% layer_dense(units = 4, activation = "relu", input_shape = c(d)) %>% layer_dense(units = 6, activation = "relu")%>% layer_dense(units = 6, activation = "relu")%>% layer_dense(units = 12, activation = "relu")%>% layer_dense(units = 6, activation = "relu")%>% layer_dense(units = 6, activation = "relu") %>% layer_dense(units = 1, activation = "linear")
  KerasModelReg %>% compile(loss = 'mean_squared_error', optimizer = "adam", metrics = "mean_squared_error")
  k_set_value(KerasModelReg$optimizer$lr, 0.03)
  history <- KerasModelReg %>% keras::fit(as.matrix(DoE), YReg, epochs = nepochs, batch_size = 10, validation_split = 0)
  #KerasModel %>% predict
  PredValRegKeras <- predict(KerasModelReg, as.matrix(DoEVal))
  
  GLM <- glm(Y ~ X1 + X2 + X3 + X4 + X1*X2 + X1*X3 + X1*X4 + X2*X3 + X2*X4 + X3*X4, data=DataClass, family = binomial())
  PredValClassGLM <- round(predict(GLM, newdata = DoEVal, type = "response"))
  svm.modelClass <- svm( Y ~ ., data = DataClass, cost = costSVMClass, gamma = gammaSVMclass)
  PredValClassSVM <- predict(svm.modelClass, DoEVal)
  KerasModelClass <- keras_model_sequential()
  KerasModelClass %>% layer_dense(units = 4, activation = "relu", input_shape = c(d)) %>% layer_dense(units = 4, activation = "relu") %>% layer_dense(units = 2, activation = "softmax")
  KerasModelClass %>% compile(loss = 'categorical_crossentropy', optimizer = "adam", metrics = "accuracy")
  k_set_value(KerasModelClass$optimizer$lr, 0.01)
  yforkeras <- to_categorical(YClass)
  history <- KerasModelClass %>% keras::fit(as.matrix(DoE), yforkeras, epochs = nepochs, batch_size = 50, validation_split = 0)
  PredValClassKeras <- apply(predict(KerasModelClass, as.matrix(DoEVal)), 1, which.max) - 1 
  ResultsVal <- data.frame("YRegVal" = YRegVal, "YClassVal"=YClassVal, "LM"=PredValRegLM, "KM"=PredValRegKM, "KerasReg"=PredValRegKeras, "SVMReg"=PredValRegSVM, "GLM"=PredValClassGLM, "KerasClass"=PredValClassKeras, "SVMClass"=as.numeric(PredValClassSVM) - 1)
  FitObjects <- list(LM = LM, KM = KM, KerasModelReg = KerasModelReg, SVMReg = svm.modelReg, GLM = GLM, KerasModelClass = KerasModelClass, SVMClass = svm.modelClass)
  return(list(ResultsVal = ResultsVal, FitObjects = FitObjects))
}


Testbed2 <- function(DoE, beta_pf, beta_num){  # assuming lin + 2fi
  # scaled on [-1,1]^d
  DoE <- (DoE + 1) / 2
  MM <- model.matrix(~.^2, data = DoE)
  MM <- cbind(MM, DoE^2 )
  names(MM)[(2+ d + choose(d,2)):(1 + 2*d + choose(d, 2))] <- paste(names(DoE), "sq", sep = "")
  y1_1 <- round(sigmoid(as.matrix(MM) %*% matrix(beta_pf, ncol = 1)))
  y1_2 <- as.matrix(MM) %*% matrix(beta_num, ncol = 1)
  Y <- data.frame(y1 = y1_1, y2 = y1_2)
  return(Y)  
}

Testbed1 <- function(DoE, w=4, t=2){
  ##########################################################################
  #
  # CANTILEVER BEAM FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUTS AND INPUTS:
  #
  # y  = c(D, S), where D = displacement and S = stress
  # xx = c(R, E, X, Y)
  # w  = width (optional)
  # t  = thickness (optional)
  #
  #########################################################################
  xx <- (DoE + 1)/2
  R <- xx[,1]*8000 + 36000
  E <- xx[,2] * 1.45e6*4 + (2.9e7 - 2*1.45e6)
  X <- xx[,3] *400 + 300
  Y <- xx[,4] * 400 + 800
  
  0*1.45e6*4 + (2.9e7 - 2*1.45e6)
  
  L   <- 100
  D_0 <- 2.2535
  
  Sterm1 <- 600*Y / (w*(t^2))
  Sterm2 <- 600*X / ((w^2)*t)
  
  S <- Sterm1 + Sterm2 
  
  Dfact1 <- 4*(L^3) / (E*w*t)
  Dfact2 <- sqrt((Y/(t^2))^2 + (X/(w^2))^2)
  
  D <- Dfact1 * Dfact2
  
  y2 <- c(D, S)
  y2 <- D
  y1 <- round(sigmoid(y2 - 4.3))
  Y <- data.frame(y1 = y1, y2 = y2)
  return(Y)
}

Testbed1old <- function(DoE){
  DoE <- (DoE + 1)/2
  y1_1 <- round(sigmoid((branin(DoE[,c(1,2)])*branin(sqrt(DoE[,c(3,4)])) - 1500)/900))$X2
  y1_2 <- (branin(DoE[,c(1,2)])*branin(DoE[,c(3,4)]))$X2
  Y <- data.frame(y1 = y1_1, y2 = y1_2)
  return(Y)  
}

sigmoid <- function(x){
  result <- 1/(1 + exp(-x))
  return(result)
}

helpfct1 <- function(x, d){
  if(is.null(dim(x))){
    x <- matrix(x, ncol = d)
  }
  return(x)
}


### implement hierarchical clustering?he
#help.search("cluster")
#DoE1 <- rbind(c(0.5, 0.5), c(0.5, -0.5), c(-0.5, 0.5), c(-0.5, -0.5))
# Minimax2(DoE1, nsim = 50000)
# 
# DoE1 <- (SPFFFscaled[[40]][,-1] + 1)/2
# plot(DoE1, xlim = c(0,1), ylim = c(0,1))
# time()
# t1 <- Sys.time()
# Minimax2(DoE1, nsim = 1e6)
# t2 <- Sys.time()
# t2 - t1
# t1 <- Sys.time()
# mM <- mMdist(DoE1, neval = 1e6)
# mM$dist
# t2 <- Sys.time()
# t2 - t1
# 
# mM$far.pt
# points(mM$far.pt[1],mM$far.pt[2], col = "red")
# as.matrix(dist(rbind(c(0.1742064, 1), DoE1)))[1,]
# DoE1
# sqrt(miM(as.matrix((SPFFF[[1]][,-1] + 1)/2), num = 1000))
# #Minimax(DoE1, nsim = 50000)
# miM(DoE1, num = 100)
# sqrt(0.5^2 + 0.5^2)
Minimax2 <- function(Design, lower = NULL, upper = NULL, nsim, seed = NULL){
  ### important: So far assuming Unit square as design space!  
  ### still to be optimized.... in time and precision
  ### assumed design space: [-1, 1]^d
  if(!is.null(seed)){set.seed(seed)}
#Design <- matrix(runif(20*2), ncol = 2)
#nsim <- 1000
  Design <- unique(Design)
  n <- nrow(Design)
  d <- ncol(Design)
  Random <- matrix(runif(nsim*d), ncol = d) * 2 - 1
#  Random <- matrix(runif(nsim*d), ncol = d)
  Dists <- dist2(Design, Random)
  ArgMins <- apply(Dists, 2, which.min)
  Mins <- apply(Dists, 2, min)
  Distsi <- numeric(n)
  for(i in 1:n){
    Distsi[i] <- max(Dists[i, which(ArgMins == i)], na.rm = TRUE)
  }  
#min(Distsi)
  results <- c(mean(Distsi), sd(Distsi), median(Distsi), min(Distsi), max(Distsi))
  names(results) <- c("mean", "sd", "median", "min", "max")
  return(results)  
#  return(c(mean(Distsi), sd(Distsi), median(Distsi), min(Distsi)))  
}

Minimax <- function(Design, lower = NULL, upper = NULL, nsim, seed = NULL){
  ### important: So far assuming Unit square as design space!  
  ### still to be optimized.... in time and precision
  ### assumed design space: [-1, 1]^d
  if(!is.null(seed)){set.seed(seed)}
  n <- nrow(Design)
  d <- ncol(Design)
  Random <- matrix(runif(nsim*d), ncol = d) * 2 - 1
#  Diststemp <- as.matrix(dist(rbind(Design, Random)))[1:n, (n + 1):(n + nsim)]
#  max(apply(Diststemp, 2, min))
#  str(Diststemp)
  MIN <- Inf 
  for(i in 1:nsim){
    Mintemp <- min(sqrt(apply((Design - matrix(Random[i,], ncol = d, nrow = n, byrow = TRUE))^2, 1, sum)))
    if(Mintemp < MIN){
      MIN <- Mintemp
    }
  }
#  print(MIN)
  return(MIN)  
}

Doptimal <- function(DoEtoEval){
  MM <- as.matrix(cbind(model.matrix(~.^2 ,DoEtoEval), DoEtoEval^2))
  Doptimal <- det(t(MM) %*% MM )
  return(Doptimal)
}


DoptimalWP <- function(DoEtoEval, WP, sigmaWP = 1, sigmae = 1){
  MM <- as.matrix(cbind(model.matrix(~.^2 ,DoEtoEval), DoEtoEval^2))
  WPCov <- diag(rep(sigmae, length(WP))) + sigmaWP * model.matrix(~. -1, data.frame(WP = as.factor(WP))) %*% t(model.matrix(~. -1, data.frame(WP = as.factor(WP))))
  #  solve(WPCov)
  Doptimal <- det(t(MM) %*% solve(WPCov) %*% MM )
  return(Doptimal)
}

Ioptimal <- function(DoEtoEval){
#DoEtoEval <- DoEs[[1]][,-1]
  MM <- as.matrix(cbind(model.matrix(~.^2 ,DoEtoEval), DoEtoEval^2))
  d <- ncol(DoEtoEval)
  dstar <- d*(d - 1)/2
  B <- cbind(c(1, rep(0, d + dstar), rep(1/3, d)), 
             rbind(rep(0, d), diag(rep(1/3, d)), matrix(0, ncol = d, nrow = dstar), matrix(0, ncol = d, nrow = d)), 
             rbind(rep(0, dstar), matrix(0, nrow = d, ncol = dstar), if(dstar == 1){1/9}else{diag(rep(1/9, dstar))}, matrix(0, ncol = dstar, nrow = d)),
             rbind(rep(1/3, d), matrix(0, ncol = d, nrow = d), matrix(0, nrow = dstar, ncol = d), (1/45)* (diag(rep(4, d)) + matrix(5, ncol = d, nrow = d))))
  Iopt <- sum(diag(solve(t(MM) %*% MM ) %*% B))
  return(Iopt)
}

IoptimalWP <- function(DoEtoEval, WP, sigmaWP = 1, sigmae = 1){
  MM <- as.matrix(cbind(model.matrix(~.^2 ,DoEtoEval), DoEtoEval^2))
  d <- ncol(DoEtoEval)
  dstar <- d*(d - 1)/2
  WPCov <- diag(rep(sigmae, length(WP))) + sigmaWP * model.matrix(~. -1, data.frame(WP = as.factor(WP))) %*% t(model.matrix(~. -1, data.frame(WP = as.factor(WP))))
  B <- cbind(c(1, rep(0, d + dstar), rep(1/3, d)), 
             rbind(rep(0, d), diag(rep(1/3, d)), matrix(0, ncol = d, nrow = dstar), matrix(0, ncol = d, nrow = d)), 
             rbind(rep(0, dstar), matrix(0, nrow = d, ncol = dstar), if(dstar == 1){1/9}else{diag(rep(1/9, dstar))}, matrix(0, ncol = dstar, nrow = d)),
             rbind(rep(1/3, d), matrix(0, ncol = d, nrow = d), matrix(0, nrow = dstar, ncol = d), (1/45)* (diag(rep(4, d)) + matrix(5, ncol = d, nrow = d))))
  IoptWP <- sum(diag(solve(t(MM) %*% solve(WPCov) %*% MM ) %*% B))
  return(IoptWP)
}

#################
### first version: Start constructing the Whole Plot FFF design, 
### then fit in the SP FFF design
#################

SplitPlotFFFDesign <- function(dWP, nWP, dSP, noverall, nsim = noverall*nWP, printsteps = 100, seed = NULL){
#dWP <- 1
#nWP <- 10
#dSP <- 1
#noverall <- 20
#printsteps <- 10
#nsim <- 40
   if(is.null(seed) == FALSE){set.seed(seed)}
    ### Step1: Create WP DoE:
  candidates <- matrix(runif(nsim*dWP), ncol= dWP)
  X <- candidates
  ncl <- nWP
  Clusterassignment <- 1:nrow(X)
  Tab <- length(table(Clusterassignment))
  ### loop to merge clusters:
  while(Tab > nWP){
#Tab > nWP
    Xsplit <- split(X, Clusterassignment)
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
#str(DistX)
#str(DistLocalM)
    diag(DistX) <- Inf
    #    str(DistX)
    MIN <- min(DistX)
    MIN.ind <- which(DistX == MIN,arr.ind = TRUE)[1,]
    ### here I need to get the indices of above
    ClusterNos[MIN.ind]
    Clusterassignment[Clusterassignment == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    Tab <- length(table(Clusterassignment))
#Tab
        if(Tab%%printsteps == 0){ 
        print(c("Stage 1, " , Tab))
      }
  }
#str(ClusterNos)
#str(unique(Clusterassignment))
    ### now programm the next stage: Do clustering for subplot level:
  Xsplit <- split(X, Clusterassignment)
  Xsplit <- lapply(Xsplit, helpfct1, d = dWP)
  WPMeans <- t(sapply(Xsplit, colMeans))
#  str(WPMeans)
  ### generate new candidates, now with WP as given:
  if(dWP == 1){
    Candidates2 <- cbind(kronecker(matrix(WPMeans, ncol = 1),matrix(1, ncol = 1, nrow = nsim/nWP)), matrix(runif(nsim*dSP), ncol = dSP))
  }else{
    Candidates2 <- cbind(kronecker(WPMeans,matrix(1, ncol = 1, nrow = nsim/nWP)), matrix(runif(nsim*dSP), ncol = dSP))
  }
#  str(Candidates2)
  ClusterassignmentWP <- kronecker(1:nWP, rep(1, nsim/nWP))
  ClusterassignmentSP <- 1:nsim
  Tab <- length(table(ClusterassignmentSP))
  ### loop to merge clusters:
    while(Tab > noverall){
#Tab > noverall
    Xsplit <- split(Candidates2, ClusterassignmentSP)
    Xsplit <- lapply(Xsplit, helpfct1, d = dSP + dWP)
    ClusterNos <- as.numeric(names(Xsplit))
    LocalMeans <- t(sapply(Xsplit, colMeans))
#str(LocalMeans)
#    if(dSP == 1){
#      LocalMeans <- t(LocalMeans)
#    }
    LocalN <- unlist(lapply(Xsplit, nrow))
    nCl <- length(LocalN)
#nCl
    Divisor <- matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow= FALSE) +   matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow= TRUE)
#str(Divisor)
    DistLocalM <- as.matrix(dist(LocalMeans))^2
#str(DistLocalM)
#str(Divisor)
    DistX <- DistLocalM/Divisor
#str(DistLocalM)
#str(DistX)
    diag(DistX) <- Inf
    #    str(DistX)
    MIN <- min(DistX)
    MIN.ind <- which(DistX == MIN,arr.ind = TRUE)[1,]
    ### here I need to get the indices of above
    ClusterNos[MIN.ind]
    ClusterassignmentSP[ClusterassignmentSP == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    Tab <- length(table(ClusterassignmentSP))
#str(Tab)
    if(Tab%%printsteps == 0){ print(c("Stage 2, " , Tab))}
  }
#ClusterNos
#length(ClusterassignmentSP)    
  ### Create final Design  
  XsplitSP <- split(Candidates2, ClusterassignmentSP)
  XsplitSP <- lapply(XsplitSP, helpfct1, d = dSP + dWP)
  SPMeans <- t(sapply(XsplitSP, colMeans))
#  pairs(SPMeans)
  ### assign SP to WPs:
  WPlist <- numeric(noverall)
  for(i in 1:noverall){
    if(dWP == 1){
      MIN.i <- which.min(as.matrix(dist(c(SPMeans[i,1:dWP], WPMeans)))[1,-1])
      SPMeans[i,dWP] <- WPMeans[MIN.i]
    }else{    
      MIN.i <- which.min(as.matrix(dist(rbind(SPMeans[i,1:dWP], WPMeans)))[1,-1])
      SPMeans[i,1:dWP] <- WPMeans[MIN.i,]
    }
    WPlist[i] <- MIN.i
  }
  DoE <- data.frame(WPlist, SPMeans)
#  pairs(SPMeans, col = WPlist)
  DoE <- DoE[sort(DoE[,1], index.return = TRUE)$ix,]
  return(DoE)
}

########################
### second option for an algorithm:
### start with SP design points and then
### for forward to refine WP Factors in a second step
########################



SplitPlotFFFDesign2 <- function(dWP, nWP, dSP, noverall, nsim = noverall*nWP, printsteps = 100, seed = NULL){
  if(is.null(seed) == FALSE){set.seed(seed)}
#    dWP <- 1
#    nWP <- 10
#    dSP <- 2
#    noverall <- 20
#    nsim <- 60
#    seed <- 1001
#    printsteps <- 5
  ### Step1: Create WP DoE:
  X <- matrix(runif(nsim*(dWP + dSP)), ncol= dWP + dSP)
  ncl <- nsim#
  Clusterassignment <- 1:nrow(X)
  Tab <- length(table(Clusterassignment))
  while(Tab > noverall){
#Tab > noverall
    Xsplit <- split(X, Clusterassignment)
    Xsplit <- lapply(Xsplit, helpfct1, d = dWP + dSP)
    ClusterNos <- as.numeric(names(Xsplit))
    LocalMeans <- t(sapply(Xsplit, colMeans))
#    if(dWP == 1){
#      LocalMeans <- t(LocalMeans)
#    }
    LocalN <- unlist(lapply(Xsplit, nrow))
    nCl <- length(LocalN)
    #print(c("ncl", ncl))
    Divisor <- matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow = FALSE) +   matrix(1 / LocalN, ncol = nCl, nrow = nCl, byrow= TRUE)
    #    str(Divisor)
    DistLocalM <- as.matrix(dist(LocalMeans))^2
#str(Divisor)
#str(DistLocalM)
    DistX <- DistLocalM/Divisor
    diag(DistX) <- Inf
    #    str(DistX)
    MIN <- min(DistX)
    MIN.ind <- which(DistX == MIN,arr.ind = TRUE)[1,]
    ### here I need to get the indices of above
    ClusterNos[MIN.ind]
    Clusterassignment[Clusterassignment == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    #    length(Clusterassignment)
    #  Clusterassignment[MIN.ind] <- min(MIN.ind)
    Tab <- length(table(Clusterassignment))
    if(Tab%%printsteps == 0){ 
      print(c("Stage 1, " , Tab))
    }
  }
  Xsplit <- split(X, Clusterassignment)
  Xsplit <- lapply(Xsplit, helpfct1, d = dWP + dSP)
  SPMeans <- t(sapply(Xsplit, colMeans))
#  str(SPMeans)
  ### second stage: go further for WP factors:
  ### continue from here
  
  ClusterassignmentWP <- 1:noverall  
  Tab <- length(table(ClusterassignmentWP))
  Candidates2 <- SPMeans[,1:dWP]
#  str(Candidates2)
  while(Tab > nWP){
#Tab
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
    #    str(Divisor)
    DistLocalM <- as.matrix(dist(LocalMeans))^2
#str(Divisor)
#str(DistLocalM)
    DistX <- DistLocalM/Divisor
    diag(DistX) <- Inf
    #    str(DistX)
    MIN <- min(DistX)
    MIN.ind <- which(DistX == MIN,arr.ind = TRUE)[1,]
    ### here I need to get the indices of above
    ClusterNos[MIN.ind]
    ClusterassignmentWP[ClusterassignmentWP == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    #    length(Clusterassignment)
    #  Clusterassignment[MIN.ind] <- min(MIN.ind)
    Tab <- length(table(ClusterassignmentWP))
    if(Tab%%printsteps == 0){ print(c("Stage 2, " , Tab))}
  }
  
#  str(Candidates2  )
  #ClusterassignmentSP
#  table(ClusterassignmentWP)
  #  pairs(Candidates2, col = ClusterassignmentWP)
  #  table(ClusterassignmentWP, ClusterassignmentSP)
  
  ### Create final Design  
  XsplitWP <- split(Candidates2, ClusterassignmentWP)
  XsplitWP <- lapply(XsplitWP, helpfct1, d = dWP)
  WPMeans <- t(sapply(XsplitWP, colMeans))
  if(dWP == 1){
    WPMeans <- t(WPMeans)
  }
#  pairs(WPMeans)
#  WPMeans  
#  row.names(WPMeans) 
#  sort(unique(ClusterassignmentWP))
  ### assign SP to WPs:
  
  #ClusterassignmetSP
  SPMeansBU <- SPMeans  
#  WPlist <- numeric(noverall)
  for(i in 1:noverall){
#    i <- 1
    SPMeans[i, 1:dWP] <- WPMeans[ClusterassignmentWP[i] == row.names(WPMeans),]
  }
  DoE <- data.frame(ClusterassignmentWP, SPMeans)
  DoE = DoE[sort(ClusterassignmentWP, index.return = TRUE)$ix,]
  row.names(DoE) = 1:noverall
#  DoE
#  pairs(SPMeans, col = ClusterassignmentWP)
  DoE <- DoE[sort(DoE[,1], index.return = TRUE)$ix,]
  return(DoE)
}


SplitPlotFFFDesign2_fast <- function(dWP, nWP, dSP, noverall, nsim = noverall*nWP, printsteps = 100, seed = NULL){
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
  
#DistX
#MIN.ind
#Clusterassignment
#LocalMeans
  while(Tab >= noverall){
    NewMean <- LocalMeans[MIN.ind[2],]*LocalN[MIN.ind[2]] / sum(LocalN[MIN.ind]) + LocalMeans[MIN.ind[1],] * LocalN[MIN.ind[1]] / sum(LocalN[MIN.ind])
    LocalN[MIN.ind[1]] <- sum(LocalN[MIN.ind])
#str(LocalN)
    RowNamesLocalMeans <- as.numeric(row.names(LocalMeans))
#LocalMeans
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
#    MIN.ind
#    LocalMeans
#row.names(LocalMeans)
    ### here I need to get the indices of above
    Clusterassignment[RowNamesLocalMeans[MIN.ind[2]]] <- Clusterassignment[RowNamesLocalMeans[MIN.ind[1]]]
#    Clusterassignment[Clusterassignment == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    if(Tab%%printsteps == 0){ 
      print(c("Stage 1, " , Tab))
    }
  }
#table(Clusterassignment)
#print(LocalMeans)
#print(dim(LocalMeans))
#print(MIN.ind)
#print(LocalN)
#length(table(Clusterassignment))
#length(LocalN)
#dim(LocalMeans)
#table(Clusterassignment)
#  Xsplit <- split(X, Clusterassignment)
#  Xsplit <- lapply(Xsplit, helpfct1, d = dWP + dSP)
#  SPMeans <- t(sapply(Xsplit, colMeans))
  #  str(SPMeans)
  ### second stage: go further for WP factors:
  ### continue from here
  SPMeans <- LocalMeans
  ClusterassignmentWP <- 1:noverall  
  Tab <- length(table(ClusterassignmentWP))
  Candidates2 <- SPMeans[,1:dWP]
  #  str(Candidates2)
  while(Tab > nWP){
    #Tab
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
    #    str(Divisor)
    DistLocalM <- as.matrix(dist(LocalMeans))^2
    #str(Divisor)
    #str(DistLocalM)
    DistX <- DistLocalM/Divisor
    diag(DistX) <- Inf
    #    str(DistX)
    MIN <- min(DistX)
    MIN.ind <- which(DistX == MIN,arr.ind = TRUE)[1,]
    ### here I need to get the indices of above
    ClusterNos[MIN.ind]
    ClusterassignmentWP[ClusterassignmentWP == max(ClusterNos[MIN.ind])] <- min(ClusterNos[MIN.ind]) 
    #    length(Clusterassignment)
    #  Clusterassignment[MIN.ind] <- min(MIN.ind)
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
  return(DoE)
}

#######################
### Kriging models, first step: pointwise corr
#######################

gaussk <- function(h, theta){
  result <- exp(-(1/2) * (h / theta)^2)
  return(result)
}
exponentialk <- function(h, theta){
  result <- exp(-h/theta)
  return(result)
}
matern3_2k <- function(h, theta){
  result <- (1+sqrt(3)*h/theta)*exp(-sqrt(3)*h/theta)
  return(result)
}
matern5_2k <- function(h, theta){
  result <- (1+sqrt(5)*h/theta+(1/3)*5*(h/theta)^2) * exp(-sqrt(5)*h/theta)
  return(result)
}


Rfunc1 <- function(H, theta, COR){
  #   Rfunc soll die Elementweise Cor Matrix enhalten. -> Nachher ?ber Prod zusammenf?gen.
  #H <- as.matrix(dist(LHD[[1]][,1]))
  #theta <- 1
  #COR <- "matern5_2k"
  ARGS <- list(H, theta)
  names(ARGS) <- c("h", "theta")
  R <- do.call(COR, args = ARGS, quote = FALSE)
  return(R)
}

LLH <- function(Param, DoE, WP, y, COR){
  d <- ncol(DoE)
  #Param <- Paramintial
  Parameters <- list(Param[1:(d+1)], Param[(d+ 2):(2*d+1)], Param[2*d+2], Param[2*d+3], Param[2*d+4])
  beta <- Parameters[[1]]
  theta <- Parameters[[2]]
  sigma_Spatial <- Parameters[[3]]
  sigma_WP <- Parameters[[4]]
  sigma <- Parameters[[5]]
  #  COR <- "gaussk"
  n <- nrow(DoE)
  trend <- model.matrix(~., data = DoE)
  CorSpatial <- matrix(1, ncol = n, nrow = n)
  for(k in 1:d){
    H <- as.matrix(dist(DoE[,k]))
    CorSpatial <- CorSpatial* Rfunc1(H, theta[k], COR)
  }
  WPF <- as.factor(WP)
  MM <- model.matrix(~WPF-1)
  CovWP <- sigma_WP^2 * MM%*% t(MM)
  C <- CorSpatial*sigma_Spatial^2 + CovWP + sigma^2*diag(n)
  #  LLH <- log((1/(sqrt(2*pi*det(R)))) * exp(t(trend%*%beta - y) %*% solve(R)%*% (trend %*% beta - y)))
  #  LLH1 <- - 0.5*log((sqrt(2*pi*det(C))))  -  t(trend%*%beta - y) %*% solve(C)%*% (trend %*% beta - y) 
  #  LLH2 <-  - 0.5 * n * log(2*pi) - 0.5 * n * det(C) - 0.5* t((y - trend%*% beta)) %*% solve(C) %*% (y - trend%*% beta)# important to note: det(C) calculated the log of the det....
  LLH <- dmvnorm(as.numeric(y),mean = as.numeric(trend%*% beta), sigma = C, log = TRUE)
  #  return(c(LLH1, LLH2, det(C)))
  return(as.numeric(LLH))
}


KrigingParamEst <- function(Cor, DoE, y, WP, ninitial){
  d <- ncol(DoE)
  n <- nrow(DoE)
  Paramintial <- runif(d + 1 + d + 3)
  #Parameters <- list(beta, theta, sigma_Spatial, sigma_WP, sigma)
  #LLH(Paramintial, DoE, WP, y,COR)
  ### alternativer for llh:
  ### Optimizing Likelihood:
  betahatOLS <- as.numeric(lm(as.numeric(y)~as.matrix(DoE))$coefficients)
  ParamintialMatrix <- cbind(matrix(betahatOLS, ncol = d + 1, nrow = n_initial, byrow = TRUE), matrix(runif((d+3)*n_initial), nrow = n_initial) + 0.01)
  
  #ParamintialMatrix <- matrix(runif((d+1+d+3)*n_initial), nrow = n_initial) + 0.01
  #LLstartValues <- apply(ParamintialMatrix, 1, LLH2, DoE = DoE, WP = WP, y = y, COR = COR)
  str(ParamintialMatrix)
  LLstartValues <- apply(ParamintialMatrix, 1, LLH, DoE = DoE, WP = WP, y = y, COR = COR)
  
  imax <- which.max(LLstartValues)
  thetastart <- ParamintialMatrix[imax,]
  lower <- c(rep(-1000, d + 1), rep(0.0001, d + 3) )
  upper <- c(rep(100000, d + 1), rep(10, d), rep(100000, 3) )
  Optim <- optim(thetastart, LLH, DoE = DoE, y = y, WP = WP, COR = Cor, method = "L-BFGS-B", control = list(fnscale = -1, trace = 2, REPORT = 5, maxit = 1000), lower = lower, upper = upper)
  #Optim
  print(Optim$convergence)
  ParOptim <- Optim$par
  return(ParOptim)
}
#ParOptim <- KrigingParamEst(COR, DoE, y, WP, ninitial)


### construct a prediction:
#newdata <- DoEVal[1:10,]
#str(DoEVal)
#str(beta)
PredFct <- function(newdata, DoE, y, Paravec){
  d <- ncol(DoE)  
  ParTemp <- list(Paravec[1:(d+1)], Paravec[(d+ 2):(2*d+1)], Paravec[2*d+2], Paravec[2*d+3], Paravec[2*d+4])
  beta <- ParTemp[[1]]
  theta <- ParTemp[[2]]
  sigma_Spatial <- ParTemp[[3]]
  sigma_WP <- ParTemp[[4]]
  sigma <- ParTemp[[5]]
  nDoE <- nrow(DoE)
  npred <- nrow(newdata)
  n <- nrow(DoE)
  CorSpatial <- matrix(1, ncol = nDoE, nrow = nDoE)
  for(k in 1:d){
    H <- as.matrix(dist(DoE[,k]))
    CorSpatial <- CorSpatial* Rfunc1(H, theta[k], COR)
  }
  WPF <- as.factor(WP)
  MM <- model.matrix(~WPF-1)
  CovWP <- sigma_WP^2 * MM%*% t(MM)
  C <- CorSpatial*sigma_Spatial^2 + CovWP + sigma^2*diag(n)
  X <- model.matrix(~., data = DoE)
  betaII <- as.numeric(solve(t(X) %*% solve(C) %*% X) %*%  t(X) %*% solve(C) %*% y)
  beta
  betaII
  
  LinearPart <- model.matrix(~., data =newdata) %*% betaII
  Corpart <- numeric(nrow(newdata))
  for(i in 1:nrow(newdata)){
    kleinr <- kleinrfunc(newdata[i,], DoE, theta, COR)
    Corpart[i] <- as.numeric(t(kleinr) %*% solve(C) %*% (y - X %*% betaII))
  }
  result <- LinearPart + Corpart
  return(result)
}


# 


kleinrfunc <- function(x, DoE, theta, cor){   #
  n <- nrow(DoE)
  d <- ncol(DoE)
  r <- rep(1, d)
  for(k in 1:d){
    h <- abs(as.numeric(x[k]) - DoE[,k])
    r <- r * Rfunc1(h, theta[k], cor)
  }
  return(r)
}

###################
### Intuitive version, just round for WP factors:
###################

#FFFRounded <- function(dWP, nWP, dSP, noverall, nsim = noverall*nWP, printsteps = 100, seed = NULL)
  
###################
### compare designs:
###################


#DoE1 <- SplitPlotFFFDesign(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)
#DoE2 <- SplitPlotFFFDesign2(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)
#DoE2



###########################
### 
###########################








#########################
### old code base:
##########################


### Constructing Likelihood, potentially faster implementation:

# LLH2 <- function(Paramintial, DoE, WP, y, COR){
#   Parameters <- list(Paramintial[1:(d+1)], Paramintial[(d+ 2):(2*d+1)], Paramintial[2*d+2], Paramintial[2*d+3], Paramintial[2*d+4])
#   #  print(Parameters)
#   beta <- Parameters[[1]]
#   theta <- Parameters[[2]]
#   #  print(theta)
#   sigma_Spatial <- Parameters[[3]]
#   sigma_WP <- Parameters[[4]]
#   sigma <- Parameters[[5]]
#   COR <- "gaussk"
#   n <- nrow(DoE)
#   d <- ncol(DoE)
#   trend <- model.matrix(~., data = DoE)
#   CorSpatial <- matrix(1, ncol = n, nrow = n)
#   for(k in 1:d){
#     #    print(k)
#     H <- as.matrix(dist(DoE[,k]))
#     #    print(H)
#     #    print(theta[k])
#     #    print(COR)
#     CorSpatial <- CorSpatial* Rfunc1(H, theta[k], COR)
#   }
#   
#   WPF <- as.factor(WP)
#   MM <- model.matrix(~WPF-1)
#   CovWP <- sigma_WP * MM%*% t(MM)
#   C <- CorSpatial*sigma_Spatial + CovWP + sigma*diag(n)
#   CholC <- chol(C)
#   M <- backsolve(t(CholC), trend)
#   Tinv_y <- backsolve(t(CholC), y, upper.tri = FALSE)
#   Q <- qr.Q(qr(M))
#   H <- Q %*% t(Q)
#   z <- Tinv_y - H %*% Tinv_y
#   v <- t(z) %*% z/n
#   loglik <- -0.5*(n * log(2 * pi * v) + 2 * sum(log(diag(CholC))) + n)
#   loglik2 <- -0.5 * (n * log(2*pi) + 2 * sum(log(diag(CholC)))) + t(Z) %*% z
#   Results <- list(loglik = loglik, Chol = CholC, C = C, trend = trend)
#   #  LLH <- log((1/(sqrt(2*pi*det(R)))) * exp(t(trend%*%beta - y) %*% solve(R)%*% (trend %*% beta - y)))
#   #  LLH1 <- - 0.5*log((sqrt(2*pi*det(R))))  -  t(trend%*%beta - y) %*% solve(R)%*% (trend %*% beta - y) 
#   #  return(c(LLH, LLH1, det(R)))
#   #  return(Results)
#   return(Results)
# }
# LLHmanual(Parameters, DoE, WP, y,COR)
# LLH2(Paramintial, DoE, WP, y,COR)$loglik
# Parameters
# 
# 
# 
# kleinRfunc2bsb.w1 <- function(x, VP, theta, COR, basisints, argmax, weighting = TRUE, wtype = "dbeta", limits){   # x kann auch wieder etwas mit mehreren Elementern sein? nein
#   ### here norder is not necessary
#   LHD <- VP$lhd
#   if(is.null(LHD)){
#     k1 <- 0
#   }else{
#     xLHD <- x[[2]]
#     k1 <- ncol(LHD)
#     thetaLHD <- theta[[1]]
#   }
#   
#   n <- nrow(VP$Coefs[[1]])
#   #  k1 <- ncol(LHD)
#   k2 <- length(VP$Coefs)
#   thetafunc <- theta[[2]]
#   RLHD <- rep(1, n)
#   ### bis hier!
#   if(k1 > 0){
#     for(k in 1:k1){
#       H <- abs(xLHD[k] - LHD[,k])
#       RLHD <- RLHD * Rfunc1(H, thetaLHD[k], COR)
#     }
#   }
#   Rfuncs <- rep(1, n)
#   for(k in 1:k2){
#     #k <- 1
#     if(weighting == TRUE){
#       if(wtype == "dbeta"){
#         Weights <- diag(dbeta(argmax, theta[[3]][k, 1], theta[[3]][k, 2]))
#         #      Weights <- sqrt(Weights / sum(Weights))
#         Weights <-   nbasis * Weights / sum(Weights)
#       }
#       if(wtype == "dlin"){
#         Weights <- diag(dlin(argmax, theta[[3]][k], limits = limits))
#         #      Weights <- sqrt(Weights / sum(Weights))
#         Weights <-   nbasis * Weights / sum(Weights)
#       }
#     }else{
#       Weights <- diag(rep(1 , ncol(basisints)))
#     }
#     H <- numeric(n)
#     for(i in 1:n){
#       H[i] <- sqrt(t(t(VP$Coefs[[k]])[,i] - t(t(x[[1]][[k]]))) %*% Weights %*%
#                      basisints %*% Weights %*% (t(VP$Coefs[[k]])[,i] - t(t(x[[1]][[k]]))))
#       #DIST <- sqrt(apply((t(Diffbeta) %*% Weights2 %*% EIGEN[[2]] %*% sqrt(diag(EIGEN[[1]])))^2, 1, sum)) 
#       
#     }
#     Rfuncs <- Rfuncs * Rfunc1(H, thetafunc[[k]], COR)
#   }
#   Rcombined <- RLHD * Rfuncs
#   return(Rcombined)
# }
# 
# 
# 
# yhatfuncbsb.w1 <- function(x, y, theta, COR = "gaussk", VP, nbasis, norder, eps.R = 0.000000000001, weighting = TRUE, wtype = "dbeta", limits){
#   # x sind die zu vorhersagenden Werte
#   ### theta wie vorher auch
#   #COR <- "matern5_2k"
#   #VP <- Best
#   #eps.R <- 0.000000001
#   #nbasis <- 5
#   #norder <- 3
#   bsb <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
#   T <- seq(0, 1, length.out = 50000)
#   Evals <- eval.basis(T, bsb)
#   basisints <- matrix(0, nbasis, nbasis)
#   for(i in 1:nbasis){
#     for(j in 1:nbasis){
#       basisints[i, j] <- mean(Evals[,i] * Evals[,j])
#     }
#   }
#   EIGEN <- eigen(basisints)
#   MAX <- apply(Evals, 2, max)
#   argmax <- numeric(nbasis)
#   for(i in 1: nbasis){
#     argmax[i] <- mean(T[(t(Evals) == MAX)[i,]])
#   }
#   if(argmax[1] == 0){
#     argmax[1] <- 0.01
#   }
#   if(argmax[nbasis] == 1){
#     argmax[nbasis] <- 1 - 0.01
#   }
#   LHD <- VP$lhd
#   if(is.null(LHD)){
#     k1 <- 0
#   }else{
#     k1 <- ncol(LHD)
#   }
#   k2 <- length(VP$Coefs)
#   n <- nrow(VP$Coefs[[1]])
#   npred <- nrow(x[[1]][[1]])
#   DM <- matrix(1, ncol = 1, nrow = n)
#   p <- 1    #  number of trend parameters
#   ### basisints berechnen
#   
#   R <- Rfunc2bsb.w1(VP, theta, COR,  EIGEN, argmax, weighting = weighting, wtype = wtype, limits = limits) + diag(eps.R, ncol = n, nrow = n)
#   Rinvs <- solve(R)
#   betadach <- solve(t(DM) %*% Rinvs %*% DM) %*% t(DM) %*% Rinvs %*% y
#   Faktor2 <- Rinvs %*% (y - DM %*% betadach)
#   sigmadachhoch2 <- (1 / (n - p)) * t((y - DM %*% betadach)) %*% Rinvs %*% (y - DM %*% betadach)
#   ydach <- numeric(npred)
#   sigmadachx <- numeric(npred)
#   for(i in (1:npred)){
#     xpred <- vector(2, mode = "list")
#     xpred[[2]] <- x[[2]][i,]
#     xpredfunc <- vector(k2, mode = "list")
#     for(k in 1:k2){
#       xpredfunc[[k]] <- x[[1]][[k]][i,]
#     }
#     xpred[[1]] <- xpredfunc
#     kleinrx <- kleinRfunc2bsb.w1(xpred, VP, theta, COR, basisints, argmax, weighting = weighting, wtype = wtype, limits = limits)
#     ydach[i] <- betadach + t(kleinrx) %*% Faktor2
#     sigmadachx[i] <- sigmadachhoch2 *  (1 - t(kleinrx) %*% Rinvs %*% kleinrx + (1 - t(kleinrx) %*% Rinvs %*% DM) * (1 / (t(DM) %*% Rinvs %*% DM)) * (1 - t(kleinrx) %*% Rinvs %*% DM))
#   }
#   result <- cbind(ydach, sigmadachx)
#   return(result)
# }
# 
# 
# 
# 
# 
# yvalpred <- as.numeric(predict.km(KM, newdata = xval, type = "UK")$mean)
# 
# str(yvalpred)
# str(yval)
# plot(yval, yvalpred)
# abline(0,1)
# ### construct an ML estimator based  on code for FANOVA paper?
# library(DiceKriging)
# ?km
# d <- 2
# covStruct.create("gauss", d, known.covparam = "None", var.names = c("X1", "X2"))
# covFFF <- function(X, WP, sigmaWP, sigmaSP, theta){
#   d <- ncol(X)
#   nWP <- length(WP)
#   n <- nrow(X)
#   for(k in 1:d){
#     
#   }
# }
# 
# LogLfuncbSPFFF <- function(theta, y, VP, WPvec, COR, argmax, eps.R = 0.000000000001, limits){
#   #### Likelihood according to Oliver.
#   #theta <- thetastart[[1]]
#   # theta <- thetastart[[1]]
#   # y
#   # VP
#   # COR
#   # EIGEN
#   # argmax
#   # eps.R
#   # weighting
#   # wtype
#   # limits
#   
#   n <- length(y)
#   R <- Rfunc2bsb.w1(VP, theta, COR, EIGEN, argmax, weighting = weighting, wtype = wtype, limits = limits) + diag(eps.R, ncol = n, nrow = n)
#   F <- matrix(rep(1, n), ncol = 1)
#   T <- chol(R)
#   M <- backsolve(t(T), F, upper.tri = FALSE)
#   Tinv_y <- backsolve(t(T), matrix(y, ncol = 1), upper.tri = FALSE)
#   Q <- qr.Q(qr(M))
#   H <- Q %*% t(Q)
#   z <- Tinv_y - H %*% Tinv_y			# linear regression properties
#   v <- t(z) %*% z / n
#   loglik <- -0.5 * (n * log(2 * pi * v) + 2 * sum (log(diag(T))) + n)
#   return(loglik)
# }
# 
# ##############################
# ### Reuse old Code base from B-spline paper:
# ##############################
# 
# system.time(thetaML.wTl <- MLoptimfuncbsb.w1(seed = SEED + 2, DoE, y, n.initialtries = nini,
#                                              limits = LIMITS, eps.R = eps.R,  COR = COR, nbasis, norder = norder,
#                                              weighting = TRUE, wtype = WTYPE, MAXIT = Nit))
# 
# MLoptimBS <- function(x, ysim, nsim = nsim, wtype = wtype, limits = limits, SEED = SEED, weighting = weighting,
#                       DoE = DoE, n.initialtries = n.initialtries,  eps.R = eps.R, COR = COR, nbasis = nbasis, norder = norder, MAXIT = MAXIT, TRACE = TRACE){ 
#   #x <- 1
#   SEED1 <- SEED
#   repeat{
#     thetatemp <- MLoptimfuncbsb.w1(seed = SEED1 + x + 3, DoE, ysim[x,], n.initialtries = n.initialtries,
#                                    limits = limits, eps.R = eps.R,  COR = COR, nbasis = nbasis, norder = norder,
#                                    weighting = weighting, wtype = wtype, MAXIT = MAXIT, TRACE = TRACE)
#     if(weighting == TRUE){
#       if(thetatemp[[4]] == 0) break
#     }
#     if(weighting == FALSE){
#       if(thetatemp[[3]] == 0) break
#     }
#     SEED1 <- SEED1 + nsim
#   }
#   return(thetatemp)
# }
# 
# #############################################
# ### GP functions
# #############################################
# 
# gaussk <- function(h, theta){
#   result <- exp(-(1/2) * (h / theta)^2)
#   return(result)
# }
# exponentialk <- function(h, theta){
#   result <- exp(-h/theta)
#   return(result)
# }
# matern3_2k <- function(h, theta){
#   result <- (1+sqrt(3)*h/theta)*exp(-sqrt(3)*h/theta)
#   return(result)
# }
# matern5_2k <- function(h, theta){
#   result <- (1+sqrt(5)*h/theta+(1/3)*5*(h/theta)^2) * exp(-sqrt(5)*h/theta)
#   return(result)
# }
# 
# 
# Rfunc1 <- function(H, theta, COR){
#   #   Rfunc soll die Elementweise Cor Matrix enhalten. -> Nachher ?ber Prod zusammenf?gen.
#   #H <- as.matrix(dist(LHD[[1]][,1]))
#   #theta <- 1
#   #COR <- "matern5_2k"
#   ARGS <- list(H, theta)
#   names(ARGS) <- c("h", "theta")
#   R <- do.call(COR, args = ARGS, quote = FALSE)
#   return(R)
# }
# 
# 
# loobsb.w1 <- function(y, theta, COR = "gaussk", VP, nbasis, norder, eps.R = 0.000000000001, weighting = TRUE, wtype = "dbeta", limits, var = FALSE){
#   # x sind die zu vorhersagenden Werte
#   ### theta wie vorher auch
#   # theta <- thetaML.w1
#   # COR = "gaussk"
#   # VP <- DoE
#   # weighting = TRUE
#   bsb <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
#   T <- seq(0, 1, length.out = 50000)
#   Evals <- eval.basis(T, bsb)
#   basisints <- matrix(0, nbasis, nbasis)
#   for(i in 1:nbasis){
#     for(j in 1:nbasis){
#       basisints[i, j] <- mean(Evals[,i] * Evals[,j])
#     }
#   }
#   EIGEN <- eigen(basisints)
#   MAX <- apply(Evals, 2, max)
#   argmax <- numeric(nbasis)
#   for(i in 1: nbasis){
#     argmax[i] <- mean(T[(t(Evals) == MAX)[i,]])
#   }
#   if(argmax[1] == 0){
#     argmax[1] <- 0.01
#   }
#   if(argmax[nbasis] == 1){
#     argmax[nbasis] <- 1 - 0.01
#   }
#   LHD <- VP$lhd
#   if(is.null(LHD)){
#     k1 <- 0
#   }else{
#     k1 <- ncol(LHD)
#   }
#   k2 <- length(VP$Coefs)
#   n <- nrow(VP$Coefs[[1]])
#   DM <- matrix(1, ncol = 1, nrow = n - 1)
#   p <- 1    #  number of trend parameters
#   ydach <- numeric(n)
#   sigmadachx <- numeric(n)
#   ### basisints berechnen
#   for(i in 1:n){
#     VPtemp <- VP
#     VPtemp$lhd <- VP$lhd[-i,]
#     for(k in 1:k2){
#       VPtemp$Coefs[[k]] <- VPtemp$Coefs[[k]][-i,]
#     }
#     R <- Rfunc2bsb.w1(VPtemp, theta, COR,  EIGEN, argmax, weighting = weighting, wtype = wtype, limits) + diag(eps.R, ncol = n - 1, nrow = n - 1)
#     Rinvs <- solve(R)
#     betadach <- solve(t(DM) %*% Rinvs %*% DM) %*% t(DM) %*% Rinvs %*% y[-i]
#     Faktor2 <- Rinvs %*% (y[-i] - DM %*% betadach)
#     sigmadachhoch2 <- (1 / (n - p)) * t((y[-i] - DM %*% betadach)) %*% Rinvs %*% (y[-i] - DM %*% betadach)
#     #ydach <- numeric(npred)
#     # sigmadachx <- numeric(npred)
#     #for(i in (1:npred)){
#     xpred <- vector(2, mode = "list")
#     xpred[[2]] <- VP$lhd[i,]
#     xpredfunc <- vector(k2, mode = "list")
#     for(k in 1:k2){
#       xpredfunc[[k]] <- VP$Coefs[[k]][i,]
#     }
#     xpred[[1]] <- xpredfunc
#     kleinrx <- kleinRfunc2bsb.w1(xpred, VPtemp, theta, COR, basisints, argmax, weighting = weighting, wtype = wtype, limits = limits)
#     ydach[i] <- betadach + t(kleinrx) %*% Faktor2
#     sigmadachx[i] <- sigmadachhoch2 *  (1 - t(kleinrx) %*% Rinvs %*% kleinrx + (1 - t(kleinrx) %*% Rinvs %*% DM) * (1 / (t(DM) %*% Rinvs %*% DM)) * (1 - t(kleinrx) %*% Rinvs %*% DM))
#     #}
#   }
#   if(var == TRUE){
#     result <- cbind(ydach, sigmadachx)
#   }
#   if(var == FALSE){
#     result <- ydach
#   }
#   return(result)
# }
# 
# Rfunc2bsb.w1 <- function(VP, theta, COR, EIGEN, argmax, weighting = TRUE, wtype = "dbeta", limits){
#   ### here, norder is not necessary
#   ### here theta is a list with 3 elements, the 3rd element describing the weighting parameters
#   ### assumption: Same weighting process for all functional inputs, but potentially different parameters
#   ### theta[[3]] is a k2 x 2 matrix, first column cosisting mean, second std dev for a normal density
#   
#   ### inconsistency: in DoE: lhd at second place, at theta, real inputs first...
#   #  VP
#   #  theta
#   #  COR
#   #  EIGEN
#   #  argmax
#   #  weighting
#   #  wtype
#   #  limits
#   
#   LHD <- VP$lhd
#   if(is.null(LHD)){
#     k1 <- 0
#   }else{
#     k1 <- ncol(LHD)
#     thetaLHD <- theta[[1]]
#   }
#   n <- nrow(VP$Coefs[[1]])
#   #  EIGEN <- eigen(basisints)
#   Comb <- combn(n, 2)
#   k2 <- length(VP$Coefs)
#   thetafunc <- theta[[2]]
#   RLHD <- matrix(1, ncol = n, nrow = n)
#   if(k1 > 0){
#     for(k in 1:k1){
#       RLHD <- RLHD * Rfunc1(as.matrix(dist(LHD[,k])), thetaLHD[k], COR)
#     }
#   }
#   Rfuncs <- matrix(1, ncol = n, nrow = n)
#   for(k in 1:k2){
#     #    k <- 1
#     if(weighting == TRUE){
#       if(wtype == "dbeta"){
#         Weights2 <- diag(dbeta(argmax, theta[[3]][k, 1], theta[[3]][k, 2]))
#         #      Weights2 <- sqrt(Weights2 / sum(Weights2))
#         Weights2 <-  nbasis *Weights2 / sum(Weights2)
#       }
#       if(wtype == "dlin"){
#         Weights2 <- diag(dlin(argmax, theta[[3]][k], limits = limits))
#         #      Weights2 <- sqrt(Weights2 / sum(Weights2))
#         Weights2 <-  nbasis *Weights2 / sum(Weights2)
#       }
#     }else{
#       Weights2 <- diag(rep(1, ncol(EIGEN[[2]])))
#     }
#     Coefstemp <- t(VP$Coefs[[k]])
#     Diffbeta <- (Coefstemp[,Comb[1,]] - Coefstemp[,Comb[2,]])
#     DIST <- sqrt(apply((t(Diffbeta) %*% Weights2 %*% EIGEN[[2]] %*% 
#                           sqrt(diag(EIGEN[[1]])))^2, 1, sum)) 
#     Dist <- matrix(0, n, n)
#     Dist[t(Comb)] <- DIST
#     Dist <- Dist + t(Dist)
#     Rfuncs <- Rfuncs * Rfunc1(Dist, thetafunc[[k]], COR)
#   }
#   Rcombined <- RLHD * Rfuncs
#   return(Rcombined)
# }
# #Rfunc2bsb.w1 <- cmpfun(Rfunc2bsb.w1)
# 
# 
# 
# kleinRfunc2bsb.w1 <- function(x, VP, theta, COR, basisints, argmax, weighting = TRUE, wtype = "dbeta", limits){   # x kann auch wieder etwas mit mehreren Elementern sein? nein
#   ### here norder is not necessary
#   LHD <- VP$lhd
#   if(is.null(LHD)){
#     k1 <- 0
#   }else{
#     xLHD <- x[[2]]
#     k1 <- ncol(LHD)
#     thetaLHD <- theta[[1]]
#   }
#   
#   n <- nrow(VP$Coefs[[1]])
#   #  k1 <- ncol(LHD)
#   k2 <- length(VP$Coefs)
#   thetafunc <- theta[[2]]
#   RLHD <- rep(1, n)
#   ### bis hier!
#   if(k1 > 0){
#     for(k in 1:k1){
#       H <- abs(xLHD[k] - LHD[,k])
#       RLHD <- RLHD * Rfunc1(H, thetaLHD[k], COR)
#     }
#   }
#   Rfuncs <- rep(1, n)
#   for(k in 1:k2){
#     #k <- 1
#     if(weighting == TRUE){
#       if(wtype == "dbeta"){
#         Weights <- diag(dbeta(argmax, theta[[3]][k, 1], theta[[3]][k, 2]))
#         #      Weights <- sqrt(Weights / sum(Weights))
#         Weights <-   nbasis * Weights / sum(Weights)
#       }
#       if(wtype == "dlin"){
#         Weights <- diag(dlin(argmax, theta[[3]][k], limits = limits))
#         #      Weights <- sqrt(Weights / sum(Weights))
#         Weights <-   nbasis * Weights / sum(Weights)
#       }
#     }else{
#       Weights <- diag(rep(1 , ncol(basisints)))
#     }
#     H <- numeric(n)
#     for(i in 1:n){
#       H[i] <- sqrt(t(t(VP$Coefs[[k]])[,i] - t(t(x[[1]][[k]]))) %*% Weights %*%
#                      basisints %*% Weights %*% (t(VP$Coefs[[k]])[,i] - t(t(x[[1]][[k]]))))
#       #DIST <- sqrt(apply((t(Diffbeta) %*% Weights2 %*% EIGEN[[2]] %*% sqrt(diag(EIGEN[[1]])))^2, 1, sum)) 
#       
#     }
#     Rfuncs <- Rfuncs * Rfunc1(H, thetafunc[[k]], COR)
#   }
#   Rcombined <- RLHD * Rfuncs
#   return(Rcombined)
# }
# 
# yhatfuncbsb.w1 <- function(x, y, theta, COR = "gaussk", VP, nbasis, norder, eps.R = 0.000000000001, weighting = TRUE, wtype = "dbeta", limits){
#   # x sind die zu vorhersagenden Werte
#   ### theta wie vorher auch
#   #COR <- "matern5_2k"
#   #VP <- Best
#   #eps.R <- 0.000000001
#   #nbasis <- 5
#   #norder <- 3
#   bsb <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
#   T <- seq(0, 1, length.out = 50000)
#   Evals <- eval.basis(T, bsb)
#   basisints <- matrix(0, nbasis, nbasis)
#   for(i in 1:nbasis){
#     for(j in 1:nbasis){
#       basisints[i, j] <- mean(Evals[,i] * Evals[,j])
#     }
#   }
#   EIGEN <- eigen(basisints)
#   MAX <- apply(Evals, 2, max)
#   argmax <- numeric(nbasis)
#   for(i in 1: nbasis){
#     argmax[i] <- mean(T[(t(Evals) == MAX)[i,]])
#   }
#   if(argmax[1] == 0){
#     argmax[1] <- 0.01
#   }
#   if(argmax[nbasis] == 1){
#     argmax[nbasis] <- 1 - 0.01
#   }
#   LHD <- VP$lhd
#   if(is.null(LHD)){
#     k1 <- 0
#   }else{
#     k1 <- ncol(LHD)
#   }
#   k2 <- length(VP$Coefs)
#   n <- nrow(VP$Coefs[[1]])
#   npred <- nrow(x[[1]][[1]])
#   DM <- matrix(1, ncol = 1, nrow = n)
#   p <- 1    #  number of trend parameters
#   ### basisints berechnen
#   
#   R <- Rfunc2bsb.w1(VP, theta, COR,  EIGEN, argmax, weighting = weighting, wtype = wtype, limits = limits) + diag(eps.R, ncol = n, nrow = n)
#   Rinvs <- solve(R)
#   betadach <- solve(t(DM) %*% Rinvs %*% DM) %*% t(DM) %*% Rinvs %*% y
#   Faktor2 <- Rinvs %*% (y - DM %*% betadach)
#   sigmadachhoch2 <- (1 / (n - p)) * t((y - DM %*% betadach)) %*% Rinvs %*% (y - DM %*% betadach)
#   ydach <- numeric(npred)
#   sigmadachx <- numeric(npred)
#   for(i in (1:npred)){
#     xpred <- vector(2, mode = "list")
#     xpred[[2]] <- x[[2]][i,]
#     xpredfunc <- vector(k2, mode = "list")
#     for(k in 1:k2){
#       xpredfunc[[k]] <- x[[1]][[k]][i,]
#     }
#     xpred[[1]] <- xpredfunc
#     kleinrx <- kleinRfunc2bsb.w1(xpred, VP, theta, COR, basisints, argmax, weighting = weighting, wtype = wtype, limits = limits)
#     ydach[i] <- betadach + t(kleinrx) %*% Faktor2
#     sigmadachx[i] <- sigmadachhoch2 *  (1 - t(kleinrx) %*% Rinvs %*% kleinrx + (1 - t(kleinrx) %*% Rinvs %*% DM) * (1 / (t(DM) %*% Rinvs %*% DM)) * (1 - t(kleinrx) %*% Rinvs %*% DM))
#   }
#   result <- cbind(ydach, sigmadachx)
#   return(result)
# }
# 
# LogLfuncbsb.w1 <- function(theta, y, VP, COR, EIGEN, argmax, eps.R = 0.000000000001, weighting, wtype = "dbeta", limits){
#   #### Likelihood according to Oliver.
#   #theta <- thetastart[[1]]
#   # theta <- thetastart[[1]]
#   # y
#   # VP
#   # COR
#   # EIGEN
#   # argmax
#   # eps.R
#   # weighting
#   # wtype
#   # limits
#   
#   n <- length(y)
#   R <- Rfunc2bsb.w1(VP, theta, COR, EIGEN, argmax, weighting = weighting, wtype = wtype, limits = limits) + diag(eps.R, ncol = n, nrow = n)
#   F <- matrix(rep(1, n), ncol = 1)
#   T <- chol(R)
#   M <- backsolve(t(T), F, upper.tri = FALSE)
#   Tinv_y <- backsolve(t(T), matrix(y, ncol = 1), upper.tri = FALSE)
#   Q <- qr.Q(qr(M))
#   H <- Q %*% t(Q)
#   z <- Tinv_y - H %*% Tinv_y			# linear regression properties
#   v <- t(z) %*% z / n
#   loglik <- -0.5 * (n * log(2 * pi * v) + 2 * sum (log(diag(T))) + n)
#   return(loglik)
# }
# 
# #LogLfuncbsb.w1(thetastart[[1]], basisints = basisints, y = y, VP = VP, COR = COR, eps.R = eps.R)
# 
# 
# LogLfunc2bsb.w1 <- function(theta, y, VP, COR, EIGEN, argmax, eps.R = 0.000000000001, weighting, wtype = "dbeta", limits){
#   #### LogLfunc2 hat als input unlist(theta)
#   ### struktur theta: zuerst die thetas f?r den LHD part, dann die
#   ### die f?r den func part, dann die f?r die Weight function
#   ### here theta has k1 + k2 + k2*2 values for dbeta and k1 + k2 + k2 for dlin
#   ### upperdlin and lowerdlin for range of a parameter for dlin
#   if(is.null(VP$lhd)){
#     k1 <- 0
#   }else{
#     k1 <- ncol(VP$lhd)
#   }
#   k2 <- length(VP$Coefs)
#   if(weighting == TRUE){
#     if(wtype == "dbeta"){  
#       if(k1 > 0){
#         thetalist <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)], matrix(theta[(k1 + k2 + 1):(k1 + k2*3)], ncol = 2, byrow = TRUE))
#       }else{
#         thetalist <- list(NULL, theta[(k1 + 1):(k1 + k2)], matrix(theta[(k1 + k2 + 1):(k1 + k2*3)], ncol = 2, byrow = TRUE))
#       }}
#     if(wtype == "dlin"){  
#       if(k1 > 0){
#         thetalist <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)], theta[(k1 + k2 + 1):(k1 + k2*2)])
#       }else{
#         thetalist <- list(NULL, theta[(k1 + 1):(k1 + k2)], theta[(k1 + k2 + 1):(k1 + k2*2)])
#       }
#     }
#   }else{
#     if(k1 > 0){
#       thetalist <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)])
#     }else{
#       thetalist <- list(NULL, theta[(k1 + 1):(k1 + k2)])
#     }
#   }
#   result <- LogLfuncbsb.w1(thetalist, y, VP, COR, EIGEN, argmax, eps.R, weighting = weighting, wtype = wtype, limits = limits)
#   return(result)
# }
# 
# MLoptimfuncbsb.w1 <- function(seed = 1, VP, y, n.initialtries = 50,
#                               limits = c(0.001, 4), eps.R = 0.00000001,  COR = "gaussk", nbasis, norder, weighting = TRUE, wtype = "dbeta", MAXIT = 100, TRACE = 2){
#   
#   # seed <- SEED + 2
#   # VP <- DoE
#   # y
#   # n.initialtries <- nini
#   # limits <- LIMITS
#   # eps.R
#   # COR
#   # nbasis
#   # norder
#   # weighting <- TRUE
#   # wtype <- "dlin"
#   # MAXIT <- Nit
#   # TRACE <- 2
#   
#   bsb <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
#   T <- seq(0, 1, length.out = 50000)
#   Evals <- eval.basis(T, bsb)
#   basisints <- matrix(0, nbasis, nbasis)
#   for(i in 1:nbasis){
#     for(j in 1:nbasis){
#       basisints[i, j] <- mean(Evals[,i] * Evals[,j])
#     }
#   }
#   EIGEN <- eigen(basisints)
#   MAX <- apply(Evals, 2, max)
#   argmax <- numeric(nbasis)
#   for(i in 1: nbasis){
#     argmax[i] <- mean(T[(t(Evals) == MAX)[i,]])
#   }
#   if(argmax[1] == 0){
#     argmax[1] <- 0.01
#   }
#   if(argmax[nbasis] == 1){
#     argmax[nbasis] <- 1 - 0.01
#   }
#   
#   set.seed(seed)
#   n <- length(y)
#   #  DM <- matrix(1, ncol = 1, nrow = n)
#   if(is.null(VP$lhd)){
#     k1 <- 0
#   }else{
#     k1 <- ncol(VP$lhd)
#   }
#   
#   k2 <- length(VP$Coefs)
#   thetastart <- vector(n.initialtries, mode = "list")
#   if(weighting == TRUE){
#     if(k1 > 0){
#       for(i in 1:n.initialtries){
#         if(wtype == "dbeta"){
#           thetastart[[i]] <- list(runif(k1, limits[1], limits[2]),
#                                   runif(k2, limits[1], limits[2]),
#                                   #              matrix(runif(k2*2, limits[1], limits[2]), ncol = 2))
#                                   matrix(runif(k2*2, limits[1], 8), ncol = 2))
#         }
#         if(wtype == "dlin"){
#           thetastart[[i]] <- list(runif(k1, limits[1], limits[2]),
#                                   runif(k2, limits[1], limits[2]),
#                                   runif(k2, 0, 2))
#         }
#       }
#     }else{
#       for(i in 1:n.initialtries){
#         thetastart[[i]] <- list(NULL,
#                                 runif(k2, limits[1], limits[2]),
#                                 matrix(runif(k2*2, limits[1], 8), ncol = 2))
#       }
#     }
#   }else{
#     if(k1 > 0){
#       for(i in 1:n.initialtries){
#         thetastart[[i]] <- list(runif(k1, limits[1], limits[2]),
#                                 runif(k2, limits[1], limits[2]))
#       }
#     }else{
#       for(i in 1:n.initialtries){
#         thetastart[[i]] <- list(NULL,
#                                 runif(k2, limits[1], limits[2]))
#       }
#     }
#   }
#   ## choosing initial points
#   #LogLfuncbsb.w1(theta[[1]], y, VP, COR, EIGEN, argmax, eps.R = eps.R, weighting, wtype, limits)
#   thetastart[[1]]
#   
#   LLinitial <- unlist(lapply(thetastart, LogLfuncbsb.w1, EIGEN = EIGEN, argmax = argmax, y = y, VP = VP, COR = COR, eps.R = eps.R, weighting = weighting, wtype = wtype, limits = limits))
#   #theta <- unlist(thetastart[[which.max(LLinitial)]])
#   thetastart <- thetastart[[which.max(LLinitial)]]
#   if(weighting == TRUE){
#     if(k1 > 0){
#       #temptheta
#       #thetastart <- temptheta
#       theta <- c(thetastart[[1]], thetastart[[2]], t(thetastart[[3]]))
#       #thetalist <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)], matrix(theta[(k1 + k2 + 1):(k1 + k2*3)], ncol = 2, byrow = TRUE))
#       
#     }else{
#       theta <- c(thetastart[[2]], t(thetastart[[3]]))
#     }
#   }else{
#     if(k1 > 0){
#       theta <- c(thetastart[[1]], thetastart[[2]])
#     }else{
#       theta <- c(thetastart[[2]])
#     }
#   }
#   theta
#   
#   ###
#   Optim <- optim(theta, LogLfunc2bsb.w1, y = y, VP = VP, COR = COR, EIGEN = EIGEN, argmax = argmax, eps.R = eps.R, weighting = weighting, wtype = wtype, limits = limits, method = "L-BFGS-B", control = list(fnscale = -1, trace = TRACE, REPORT = 5, maxit = MAXIT),
#                  lower = limits[1], upper = limits[2])
#   print(Optim$convergence)
#   theta <- Optim$par
#   if(weighting == TRUE){
#     if(k1 > 0){
#       if(wtype == "dbeta"){
#         theta <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)], matrix(theta[(k1 + k2 + 1):(k1 + k2*3)], ncol = 2, byrow = TRUE), Optim$convergence)
#       }
#       if(wtype == "dlin"){
#         theta <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)], theta[(k1 + k2 + 1):(k1 + k2*2)], Optim$convergence)
#       }
#     }else{
#       if(wtype == "dbeta"){
#         theta <- list(NULL, theta[(k1 + 1):(k1 + k2)], matrix(theta[(k1 + k2 + 1):(k1 + k2*3)], ncol = 2, byrow = TRUE),  Optim$convergence)
#       }
#       if(wtype == "dlin"){
#         theta <- list(NULL, theta[(k1 + 1):(k1 + k2)], theta[(k1 + k2 + 1):(k1 + k2*2)],  Optim$convergence)
#       }
#     }
#   }else{
#     if(k1 > 0){
#       theta <- list(theta[1:k1], theta[(k1 + 1):(k1 + k2)], Optim$convergence)
#     }else{
#       theta <- list(NULL, theta[(k1 + 1):(k1 + k2)], Optim$convergence)
#     }
#   }
#   return(theta)
# }
# 
# 
# 
# d
# CorSpatial <- matrix(1, ncol = n, nrow = n)
# for(k in 1:d){
#   H <- as.matrix(dist(DoE[,k]))
#   CorSpatial <- CorSpatial* Rfunc1(H, theta, COR)
# }
# str(CorSpatial)
# #CorM <- Rfunc1(H, theta, COR)
# str(CorSpatial)  
# WPF <- as.factor(WP)
# MM <- model.matrix(~WPF-1)
# str(MM)
# CovWP <- sigma_WP * MM%*% t(MM)
# R <- CorSpatial*sigma_Spatial + CovWP + sigma*diag(n)
# str(R)
# 
# ### next step: create trend function
# trend <- model.matrix(~., data = DoE)
# str(trend)
# ###############
# ### put together likelihood:
# F <- trend
# T <- chol(R)
# d
# beta<- matrix(runif(d+1), ncol = 1)
# theta <- runif(d)
