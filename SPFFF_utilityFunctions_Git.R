#### Utility functions to evaluate designs, create test functions and make predictions based on different methods:
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


ModelPreds <- function(DoE, WP, YReg, YClass, DoEVal, YRegVal, YClassVal, nepochs = 2, nugget = 0.00001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1){  
### utitlity function to collect different predictions
  d <- ncol(DoE)
  DataReg <- data.frame(Y = YReg, DoE)
  DataClass <- data.frame(Y=as.factor(YClass), DoE)
  LM <- lm(Y ~ X1 + X2 + X3 + X4 + X1*X2 + X1*X3 + X1*X4 + X2*X3 + X2*X4 + X3*X4 + I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2), data=DataReg)
  PredValRegLM <- predict(LM, newdata = DoEVal)
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


Minimax <- function(Design, lower = NULL, upper = NULL, nsim, seed = NULL){
  ### assumed design space: [-1, 1]^d
  if(!is.null(seed)){set.seed(seed)}
  Design <- unique(Design)
  n <- nrow(Design)
  d <- ncol(Design)
  Random <- matrix(runif(nsim*d), ncol = d) * 2 - 1
  Dists <- dist2(Design, Random)
  ArgMins <- apply(Dists, 2, which.min)
  Mins <- apply(Dists, 2, min)
  Distsi <- numeric(n)
  for(i in 1:n){
    Distsi[i] <- max(Dists[i, which(ArgMins == i)], na.rm = TRUE)
  }  
  results <- c(mean(Distsi), sd(Distsi), median(Distsi), min(Distsi), max(Distsi))
  names(results) <- c("mean", "sd", "median", "min", "max")
  return(results)  
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