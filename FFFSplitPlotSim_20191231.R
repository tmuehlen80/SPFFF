#?source
#setwd("/home/muehlenth/Documents/08_Learning/SplitPlotFastFlexibleFillingDesign")
#setwd("/home/rstudio/SPFFF")
#install.packages('mvtnorm')
#install.packages('acebayes')
#ssetwd("~/08_Learning/Split Plot Fast Flexible Filling Design")

### test 2

getwd()
source(file = "FFFSplitplotFunctions_20200103.R")
##########################
### Define the number of competitors:
### The 2 algorithms
### D optimal RSM Split Plot DoE
### Maximin LHS
### unconstrained FFF Design
### acebayes design
##########################

n <- rep(c(20, 20, sort(rep(seq(25, 50, 5), 3))), 2)
dWP <- c(rep(1, length(n)/2), rep(2, length(n)/2))
dSP <- c(rep(1, length(n)/2), rep(2, length(n)/2))
nWP <- rep(c(8, 12, rep(c(8, 12, 16), 6)), 2)

nsim <- n*3
seed <- 1001
SPFFF <- vector(mode = "list", length = length(n))
SPFFFscaled <- vector(mode = "list", length = length(n))
timing <- numeric(length(n))
for(i in (1:length(n))){
#i <- 1
  print(i)
  time1 <- Sys.time()
  temp <- SplitPlotFFFDesign2_fast(dWP = dWP[i], nWP = nWP[i], dSP = dSP[i], noverall = n[i], nsim = nsim[i], seed = seed, printsteps = 40)
  temp[,2:(dWP[i] + dSP[i] + 1)] <- temp[,2:(dWP[i] + dSP[i] + 1)] * 2 - 1
  SPFFF[[i]] <- temp
  Mins <- apply(temp[,-1], 2, min)
  Maxs <- apply(temp[,-1], 2, max)
  MINS <- matrix(Mins, nrow = n[i], ncol = dWP[i] + dSP[i], byrow = TRUE)
  MAXS <- matrix(Maxs, nrow = n[i], ncol = dWP[i] + dSP[i], byrow = TRUE)
  temp2 <- data.frame(ClusterassignmentWP = temp[,1],(temp[,-1] - MINS) * 2/ (MAXS - MINS) - 1)
  SPFFFscaled[[i]] <- temp2
  time2 <- Sys.time()
  timing[i] <- time2 - time1
  print(i)
  print(time2 - time1)
}
### save created designs to file:

###  read I optimal SP designs:

SPiOpt <- vector(mode = "list", length = length(n))
for(i in (1:length(n))){
  temp <- read.table(paste("Custom Design_I_d", dWP[i] + dSP[i],"WP",nWP[i],"n", n[i], "_20191229.csv", sep = ""), sep = ",", header = TRUE, dec = ",")[,1:(dWP[i] + dSP[i] + 1)]
  SPiOpt[[i]] <- temp
}

### create LHS and read FFFs:

nLHS <- sort(rep(seq(20, 50, 5), 2))
dLHS <- rep(c(2,4), 7)
LHSs <- vector(length(nLHS), mode = "list")
FFFs <- vector(length(nLHS), mode = "list")
Iopts <- vector(length(nLHS), mode = "list")
set.seed(1258126312)
for(i in (1:length(nLHS))){
print(i)
  temp <- data.frame(maximinLHS(nLHS[i], dLHS[i], method = "iterative", maxIter = 10000) * 2 - 1)
  temp <- data.frame(WP = 1, temp)
  LHSs[[i]] <- temp
  temp <- read.table(paste("Fast Flexible Filling Design_d", dLHS[i], "n", nLHS[i], "_20191229.csv", sep = ""), sep = ",", header = TRUE, dec = ",")[,1:dLHS[i]]
  temp <- data.frame(WP = 1, temp)
  FFFs[[i]] <- temp
}

#save.image("DoEs_20200104.Rdata")
load("DoEs_20200105.Rdata")

ls()
n

#####################
### check some things
#####################
str(SPFFF[[1]])
apply(SPFFF[[3]], 2, summary)
str(SPiOpt[[3]])
apply(SPiOpt[[3]], 2, summary)
length(SPFFF)
length(SPiOpt)
str(n)
dWP
dSP
n
nWP

str(LHSs[[1]])
str(FFFs[[1]])
apply(LHSs[[3]], 2, summary)
apply(FFFs[[3]], 2, summary)
length(FFFs)
str(nLHS)
nLHS
dLHS


#####################
# calculate optimality criteria
#####################

source(file = "FFFSplitplotFunctions_20200103.R")

nmM <- 500000
miMnum <- 20
IoptimalsWP <- Ioptimals <- PhiMMs <- DoptimalsWP <- Doptimals <- MaxiMinsUnique <- Maximins <- Minimaxs <- Minimaxs3 <- numeric(length(n))
IoptimalsWPScaled <- IoptimalsScaled <- PhiMMsScaled <- DoptimalsWPScaled <- DoptimalsScaled <- MaxiMinsUniqueScaled <- MaximinsScaled <- MinimaxsScaled <- Minimaxs3Scaled <- numeric(length(n))
IoptimalsWPSP <- IoptimalsSP <- PhiMMsSP <- DoptimalsWPSP <- DoptimalsSP <- MaxiMinsUniqueSP <- MaximinsSP <- MinimaxsSP <- Minimaxs3SP <- numeric(length(n))
Minimaxs2 <- matrix(0, nrow = 5, ncol = length(n))
Minimaxs2Scaled <- matrix(0, nrow = 5, ncol = length(n))
Minimaxs2SP <- matrix(0, nrow = 5, ncol = length(n))
for(i in 1:(length(n))){
  print(i)
  Minimaxs3[i] <- miM(as.matrix(SPFFF[[i]][,-1]), num = miMnum)
  Maximins[i] <- min(dist(SPFFF[[i]][,-1]))
  MaxiMinsUnique[i] <- min(dist(unique(SPFFF[[i]][,-1])))
  Doptimals[i] <- Doptimal(SPFFF[[i]][,-1])
  PhiMMs[i] <- sum(1/dist(SPFFF[[i]][,-1])^2)
  Ioptimals[i] <- Ioptimal(SPFFF[[i]][,-1])
  Minimaxs2[, i] <- Minimax2(SPFFF[[i]][,-1], nsim = nmM)
  DoptimalsWP[i] <- DoptimalWP(SPFFF[[i]][,-1], SPFFF[[i]][,1])
  IoptimalsWP[i] <- IoptimalWP(SPFFF[[i]][,-1],SPFFF[[i]][,1])
  Minimaxs[i] <- Minimaxs2[5,i]
  
#  MinimaxsScaled[i] <- Minimax(SPFFFscaled[[i]][,-1], nsim = 200, seed = 1)
  Minimaxs3Scaled[i] <- miM(as.matrix(SPFFFscaled[[i]][,-1]), num = miMnum)
  MaximinsScaled[i] <- min(dist(SPFFFscaled[[i]][,-1]))
  MaxiMinsUniqueScaled[i] <- min(dist(unique(SPFFFscaled[[i]][,-1])))
  DoptimalsScaled[i] <- Doptimal(SPFFFscaled[[i]][,-1])
  PhiMMsScaled[i] <- sum(1/dist(SPFFFscaled[[i]][,-1])^2)
  IoptimalsScaled[i] <- Ioptimal(SPFFFscaled[[i]][,-1])
  Minimaxs2Scaled[, i] <- Minimax2(SPFFFscaled[[i]][,-1], nsim = nmM)
  DoptimalsWPScaled[i] <- DoptimalWP(SPFFFscaled[[i]][,-1], SPFFFscaled[[i]][,1])
  IoptimalsWPScaled[i] <- IoptimalWP(SPFFFscaled[[i]][,-1],SPFFFscaled[[i]][,1])
  MinimaxsScaled[i] <- Minimaxs2Scaled[5,i]
  
#  MinimaxsSP[i] <- Minimax(SPiOpt[[i]][,-1], nsim = 200, seed = 1)
  Minimaxs3SP[i] <- miM(as.matrix(SPiOpt[[i]][,-1]), num = miMnum)
  MaximinsSP[i] <- min(dist(SPiOpt[[i]][,-1]))
  MaxiMinsUniqueSP[i] <- min(dist(unique(SPiOpt[[i]][,-1])))
  DoptimalsSP[i] <- Doptimal(SPiOpt[[i]][,-1])
  PhiMMsSP[i] <- sum(1/dist(SPiOpt[[i]][,-1])^2)
  IoptimalsSP[i] <- Ioptimal(SPiOpt[[i]][,-1])
  Minimaxs2SP[, i] <- Minimax2(SPiOpt[[i]][,-1], nsim = nmM)
  DoptimalsWPSP[i] <- DoptimalWP(SPiOpt[[i]][,-1], SPiOpt[[i]][,1])
  IoptimalsWPSP[i] <- IoptimalWP(SPiOpt[[i]][,-1],SPiOpt[[i]][,1])
  MinimaxsSP[i] <- Minimaxs2SP[5,i]
}

#Minimax2(SPFFF[[i]][,-1], nsim = nmM)
#miM(as.matrix(SPFFF[[i]][,-1]))
#miM(as.matrix(SPFFFscaled[[i]][,-1]))
#mMdist(as.matrix(SPFFF[[i]][,-1]))

IoptimalsLHS <- PhiMMsLHS <- DoptimalsLHS <- MaxiMinsUniqueLHS <- MaximinsLHS <- MinimaxsLHS <- Minimaxs3LHS <- numeric(length(nLHS))
IoptimalsFFF <- PhiMMsFFF <- DoptimalsFFF <- MaxiMinsUniqueFFF <- MaximinsFFF <- MinimaxsFFF <- Minimaxs3FFF <- numeric(length(nLHS))
Minimaxs2LHS <- matrix(0, nrow = 5, ncol = length(nLHS))
Minimaxs2FFF <- matrix(0, nrow = 5, ncol = length(nLHS))
for(i in 1:(length(nLHS))){
  print(i)
#  MinimaxsLHS[i] <- Minimax(LHSs[[i]][,-1], nsim = 200, seed = 1)
  Minimaxs3LHS[i] <- miM(as.matrix(LHSs[[i]][,-1]), num = miMnum)
  MaximinsLHS[i] <- min(dist(LHSs[[i]][,-1]))
  MaxiMinsUniqueLHS[i] <- min(dist(unique(LHSs[[i]][,-1])))
  DoptimalsLHS[i] <- Doptimal(LHSs[[i]][,-1])
  PhiMMsLHS[i] <- sum(1/dist(LHSs[[i]][,-1])^2)
  IoptimalsLHS[i] <- Ioptimal(LHSs[[i]][,-1])
  Minimaxs2LHS[, i] <- Minimax2(LHSs[[i]][,-1], nsim = nmM)
  MinimaxsLHS[i] <- Minimaxs2LHS[5,i]
  
#  MinimaxsFFF[i] <- Minimax(FFFs[[i]][,-1], nsim = 200, seed = 1)
  Minimaxs3FFF[i] <- miM(as.matrix(FFFs[[i]][,-1]), num = miMnum)
  MaximinsFFF[i] <- min(dist(FFFs[[i]][,-1]))
  MaxiMinsUniqueFFF[i] <- min(dist(unique(FFFs[[i]][,-1])))
  DoptimalsFFF[i] <- Doptimal(FFFs[[i]][,-1])
  PhiMMsFFF[i] <- sum(1/dist(FFFs[[i]][,-1])^2)
  IoptimalsFFF[i] <- Ioptimal(FFFs[[i]][,-1])
  Minimaxs2FFF[, i] <- Minimax2(FFFs[[i]][,-1], nsim = nmM)
  MinimaxsFFF[i] <- Minimaxs2FFF[5,i]
}

##########################
### plot optimality criteria
##########################

LHS_criteria <- data.frame(n = nLHS, d = dLHS, Phi = PhiMMsLHS, Minimax = MinimaxsLHS, Minimax2 = Minimaxs2LHS[4,], Minimax3 = Minimaxs3LHS, Maximin = MaximinsLHS, MaximinUnique = MaxiMinsUniqueLHS, D = DoptimalsLHS, I = IoptimalsLHS)
FFF_criteria <- data.frame(n = nLHS, d = dLHS, Phi = PhiMMsFFF, Minimax = MinimaxsFFF, Minimax2 = Minimaxs2FFF[4,], Minimax3 = Minimaxs3FFF, Maximin = MaximinsFFF, MaximinUnique = MaxiMinsUniqueFFF, D = DoptimalsFFF, I = IoptimalsFFF)
SPFFF_criteria <- data.frame(n, nWP, dWP, dSP, Phi = PhiMMs, Minimax = Minimaxs, Minimax2 = Minimaxs2[4,], Minimax3 = Minimaxs3, Maximin = Maximins, MaximinUnique = MaxiMinsUnique, D = Doptimals, I = Ioptimals, ISP = IoptimalsWP, DSP = DoptimalsWP)
SPFFF_criteria['d'] <- SPFFF_criteria$dWP + SPFFF_criteria$dSP
SPFFFscaled_criteria <- data.frame(n, nWP, dWP, dSP, Phi = PhiMMsScaled, Minimax = MinimaxsScaled, Minimax2 = Minimaxs2Scaled[4,], Minimax3 = Minimaxs3Scaled, Maximin = MaximinsScaled, MaximinUnique = MaxiMinsUniqueScaled, D = DoptimalsScaled, I = IoptimalsScaled, ISP = IoptimalsWPScaled, DSP = DoptimalsWPScaled)
SPFFFscaled_criteria['d'] <- SPFFFscaled_criteria$dWP + SPFFFscaled_criteria$dSP
I_criteria <- data.frame(n, nWP, dWP, dSP, Phi = PhiMMsSP, Minimax = MinimaxsSP, Minimax2 = Minimaxs2SP[4,], Minimax3 = Minimaxs3SP, Maximin = MaximinsSP, MaximinUnique = MaxiMinsUniqueSP, D = DoptimalsSP, I = IoptimalsSP, ISP = IoptimalsWPSP, DSP = DoptimalsWPSP)
I_criteria['d'] <- I_criteria$dWP + I_criteria$dSP

#save.image("OptCriteria_20200106.Rdata")



### Check Minimax and Minimax2

plot(SPFFF_criteria$Minimax2, SPFFF_criteria$Minimax3)
### Maximin:

#dtemp <- 4
for(dtemp in c(2, 4)){ 
  YLIM <- range(c(FFF_criteria$Maximin[FFF_criteria$d == dtemp], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp], I_criteria$Maximin[I_criteria$d == dtemp], LHS_criteria$Maximin[LHS_criteria$d == dtemp]))
  YLIM[2] <- YLIM[2] * 1.3
  png(paste("Maximin_d", dtemp, "_20200102.png", sep = ""), 500, 500)
  plot(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$Maximin[FFF_criteria$d == dtemp], ylim = YLIM, xlab = "n", ylab = "Maximin")
  lines(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$Maximin[FFF_criteria$d == dtemp])
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$Maximin[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$Maximin[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  points(I_criteria$n[I_criteria$d == dtemp], I_criteria$Maximin[I_criteria$d == dtemp], col = "blue")
  lines(I_criteria$n[I_criteria$d == dtemp], I_criteria$Maximin[I_criteria$d == dtemp], col = "blue")
  points(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$Maximin[LHS_criteria$d == dtemp], col = "green")
  lines(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$Maximin[LHS_criteria$d == dtemp], col = "green")
  title(main = paste("Maximin Criterion, d =", dtemp))
  legend("topright", col = c("black", rep("red", 3), rep("orange", 3), "blue", "green"), legend = c("FFF", "SPFFF8", "SPFFFF12", "SPFFF16", "SPFFFs8", "SPFFFs12", "SPFFFs16", "I", "LHS"), pch = c(1, 8, 12, 16, 8, 12, 16, 1, 1))
  dev.off()
}

### Phi:
for(dtemp in c(2, 4)){ 
  YLIM <- range(c(FFF_criteria$Phi[FFF_criteria$d == dtemp], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp], LHS_criteria$Phi[LHS_criteria$d == dtemp]))
  YLIM[2] <- YLIM[2] * 1.3
  png(paste("Phi_d", dtemp, "_20200102.png", sep = ""), 500, 500)
  plot(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$Phi[FFF_criteria$d == dtemp], ylim = YLIM, xlab = "n", ylab = "Phi")
  lines(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$Phi[FFF_criteria$d == dtemp])
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$Phi[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$Phi[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  #points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$Phi[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  #lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$Phi[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  #points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$Phi[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  #lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$Phi[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  #points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$Phi[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  #lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$Phi[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  points(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$Phi[LHS_criteria$d == dtemp], col = "green")
  lines(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$Phi[LHS_criteria$d == dtemp], col = "green")
  title(main = paste("Phi Criterion, d =", dtemp))
  legend("topright", col = c("black", rep("red", 3), rep("orange", 3), "green"), legend = c("FFF", "SPFFF8", "SPFFFF12", "SPFFF16", "SPFFFs8", "SPFFFs12", "SPFFFs16", "LHS"), pch = c(1, 8, 12, 16, 8, 12, 16, 1))
  dev.off()
}
### Minimax:
#dtemp <- 4
for(dtemp in c(2, 4)){ 
  YLIM <- range(c(FFF_criteria$Minimax[FFF_criteria$d == dtemp], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp], I_criteria$Minimax[I_criteria$d == dtemp], LHS_criteria$Minimax[LHS_criteria$d == dtemp]))
  YLIM[2] <- YLIM[2] * 1.3
  png(paste("Minimax_d", dtemp, "_20200102.png", sep = ""), 500, 500)
  plot(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$Minimax[FFF_criteria$d == dtemp], ylim = YLIM, xlab = "n", ylab = "Minimax")
  lines(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$Minimax[FFF_criteria$d == dtemp])
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$Minimax[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$Minimax[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$Minimax[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$Minimax[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$Minimax[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$Minimax[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$Minimax[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$Minimax[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  points(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$Minimax[LHS_criteria$d == dtemp], col = "green")
  lines(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$Minimax[LHS_criteria$d == dtemp], col = "green")
  title(main = paste("Minimax Criterion, d =", dtemp))
  legend("topright", col = c("black", rep("red", 3), rep("orange", 3), rep("blue", 3), "green"), legend = c("FFF", "SPFFF8", "SPFFFF12", "SPFFF16", "SPFFFs8", "SPFFFs12", "SPFFFs16", "I8", "I12", "I16", "LHS"), pch = c(1, 8, 12, 16, 8, 12, 16, 8, 12, 16, 1))
  dev.off()
}

### I:
#dtemp <- 4
for(dtemp in c(2, 4)){ 
  YLIM <- range(c(FFF_criteria$I[FFF_criteria$d == dtemp], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp], SPFFF_criteria$I[SPFFF_criteria$d == dtemp], I_criteria$I[I_criteria$d == dtemp], LHS_criteria$I[LHS_criteria$d == dtemp]))
  YLIM[2] <- YLIM[2] * 1.3
  png(paste("I_d", dtemp, "_20200102.png", sep = ""), 500, 500)
  plot(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$I[FFF_criteria$d == dtemp], ylim = YLIM, xlab = "n", ylab = "I")
  lines(FFF_criteria$n[FFF_criteria$d == dtemp], FFF_criteria$I[FFF_criteria$d == dtemp])
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$I[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$I[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$I[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$I[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$I[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$I[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$I[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$I[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$I[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$I[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$I[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$I[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$I[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  points(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$I[LHS_criteria$d == dtemp], col = "green")
  lines(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$I[LHS_criteria$d == dtemp], col = "green")
  title(main = paste("I Criterion, d =", dtemp))
  legend("topright", col = c("black", rep("red", 3), rep("orange", 3), rep("blue", 3), "green"), legend = c("FFF", "SPFFF8", "SPFFFF12", "SPFFF16", "SPFFFs8", "SPFFFs12", "SPFFFs16", "I8", "I12", "I16", "LHS"), pch = c(1, 8, 12, 16, 8, 12, 16, 8, 12, 16, 1))
  dev.off()
}

### I SP:
dtemp <- 4
for(dtemp in c(2, 4)){ 
  YLIM <- range(c(FFF_criteria$ISP[FFF_criteria$d == dtemp], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp], I_criteria$ISP[I_criteria$d == dtemp], LHS_criteria$ISP[LHS_criteria$d == dtemp]))
  YLIM[2] <- YLIM[2] * 1.3
  png(paste("ISP_d", dtemp, "_20200102.png", sep = ""), 500, 500)
  plot(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], ylim = YLIM, xlab = "n", ylab = "ISP", col = "red", pch = 8)
  #points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 8], col = "red", pch = 8)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 12], col = "red", pch = 12)
  points(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  lines(SPFFF_criteria$n[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], SPFFF_criteria$ISP[SPFFF_criteria$d == dtemp & SPFFF_criteria$nWP == 16], col = "red", pch = 16)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 8], col = "orange", pch = 8)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 12], col = "orange", pch = 12)
  points(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  lines(SPFFFscaled_criteria$n[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], SPFFFscaled_criteria$ISP[SPFFFscaled_criteria$d == dtemp & SPFFFscaled_criteria$nWP == 16], col = "orange", pch = 16)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$ISP[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 8], I_criteria$ISP[I_criteria$d == dtemp & I_criteria$nWP == 8], col = "blue", pch = 8)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$ISP[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 12], I_criteria$ISP[I_criteria$d == dtemp & I_criteria$nWP == 12], col = "blue", pch = 12)
  points(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$ISP[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  lines(I_criteria$n[I_criteria$d == dtemp & I_criteria$nWP == 16], I_criteria$ISP[I_criteria$d == dtemp & I_criteria$nWP == 16], col = "blue", pch = 16)
  points(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$ISP[LHS_criteria$d == dtemp], col = "green")
  lines(LHS_criteria$n[LHS_criteria$d == dtemp], LHS_criteria$ISP[LHS_criteria$d == dtemp], col = "green")
  title(main = paste("ISP Criterion, d =", dtemp))
  legend("topright", col = c(rep("red", 3), rep("orange", 3), rep("blue", 3)), legend = c( "SPFFF8", "SPFFF12", "SPFFF16", "SPFFFs8", "SPFFFs12", "SPFFFs16", "I8", "I12", "I16"), pch = c(8, 12, 16, 8, 12, 16, 8, 12, 16, 1))
  dev.off()
}


#########################
### make prediction example comparison for SPFFF and FFF
#########################
source(file = "FFFSplitplotFunctions_20200103.R")
#Testbed1(SPFFF[[40]][,-1])

d <- 4
nval <- 100000
set.seed(7)
nepochs = 100
beta_pf <- round(runif(2*d + 1 + choose(d,2))*4-2)
beta_num <- round(runif(2*d + 1 + choose(d,2))*4-2)
beta_pf
beta_num

DoEVal <- data.frame(matrix(runif(nval*d), ncol = d) * 2 - 1)
Y1.Val <- Testbed1(DoEVal)
YRegVal1 <- Y1.Val[,2]
YClassVal1 <- Y1.Val[,1]  
Y2.Val <- Testbed2(DoEVal, beta_pf, beta_num)
YRegVal2 <- Y2.Val[,2]
YClassVal2 <- Y2.Val[,1]  
str(DoEVal)
DoEindicesSPFFF <- which(dWP == 2)
DoEindicesFFF <- which(FFF_criteria$d == 4)
YsSPFFF1 <- YsSPFFF2 <- vector(mode = "list", length = length(DoEindicesSPFFF))

ModelOutputFFF1 <- ModelOutputFFF2 <- YsFFF1 <- YsFFF2 <- vector(mode = "list", length = length(DoEindicesFFF))
ModelOutputSPFFF1 <- ModelOutputSPFFF2 <- ModelOutputSPFFFScaled1 <- ModelOutputSPFFFScaled2 <- YsSPFFFScaled1 <- YsSPFFFScaled2 <- vector(mode = "list", length = length(DoEindicesSPFFF))
Accs2SPFFFScaled <- Accs2SPFFF <- Accs1SPFFFScaled <- Accs1SPFFF <- matrix(0, ncol = 3, nrow = length(DoEindicesSPFFF))
Accs1FFF <- Accs2FFF <- matrix(0, ncol = 3, nrow = length(DoEindicesFFF))
Rsq2SPFFFScaled <- Rsq2SPFFF <- Rsq1SPFFFScaled <- Rsq1SPFFF <- matrix(0, ncol = 4, nrow = length(DoEindicesSPFFF))
Rsq1FFF <- Rsq2FFF <- matrix(0, ncol = 4, nrow = length(DoEindicesFFF))
### only use the 4d DoEs:

for(i in (1:length(DoEindicesSPFFF))){
#i <- 1
  YsSPFFF1[[i]] <- Testbed1(SPFFF[[DoEindicesSPFFF[i]]][,-1])   
  YsSPFFFScaled1[[i]] <- Testbed1(SPFFFscaled[[DoEindicesSPFFF[i]]][,-1])   
  YsSPFFF2[[i]] <- Testbed2(SPFFF[[DoEindicesSPFFF[i]]][,-1], beta_pf, beta_num)   
  YsSPFFFScaled2[[i]] <- Testbed2(SPFFFscaled[[DoEindicesSPFFF[i]]][,-1], beta_pf, beta_num)   
  ModelOutputSPFFF1[[i]] <- ModelPreds(DoE=SPFFF[[DoEindicesSPFFF[i]]][,-1], WP=SPFFF[[DoEindicesSPFFF[i]]][,1], YReg=YsSPFFF1[[i]][,2], YsSPFFF1[[i]][,1], DoEVal, YRegVal1, YClassVal1, nepochs = nepochs, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)
  ModelOutputSPFFF2[[i]] <- ModelPreds(DoE=SPFFF[[DoEindicesSPFFF[i]]][,-1], WP=SPFFF[[DoEindicesSPFFF[i]]][,1], YReg=YsSPFFF2[[i]][,2], YsSPFFF2[[i]][,1], DoEVal, YRegVal2, YClassVal2, nepochs = nepochs, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)
  ModelOutputSPFFFScaled1[[i]] <- ModelPreds(DoE=SPFFF[[DoEindicesSPFFF[i]]][,-1], WP=SPFFFscaled[[DoEindicesSPFFF[i]]][,1], YReg=YsSPFFFScaled1[[i]][,2], YsSPFFFScaled1[[i]][,1], DoEVal, YRegVal1, YClassVal1, nepochs = nepochs, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)
  ModelOutputSPFFFScaled2[[i]] <- ModelPreds(DoE=SPFFF[[DoEindicesSPFFF[i]]][,-1], WP=SPFFFscaled[[DoEindicesSPFFF[i]]][,1], YReg=YsSPFFFScaled2[[i]][,2], YsSPFFFScaled2[[i]][,1], DoEVal, YRegVal2, YClassVal2, nepochs = nepochs, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)
  Rsq1SPFFF[i, 1] <- (cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$LM))^2
  Rsq1SPFFF[i, 2] <- (cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$KM))^2
  Rsq1SPFFF[i, 3] <- (cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasReg))^2
  Rsq1SPFFF[i, 4] <- (cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMReg))^2
  Accs1SPFFF[i, 1] <- sum(diag(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$GLM))
  Accs1SPFFF[i, 2] <- sum(diag(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasClass))
  Accs1SPFFF[i, 3] <- sum(diag(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMClass))
  Rsq2SPFFF[i, 1] <- (cor(ModelOutputSPFFF2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF2[[i]]$ResultsVal$LM))^2
  Rsq2SPFFF[i, 2] <- (cor(ModelOutputSPFFF2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF2[[i]]$ResultsVal$KM))^2
  Rsq2SPFFF[i, 3] <- (cor(ModelOutputSPFFF2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF2[[i]]$ResultsVal$KerasReg))^2
  Rsq2SPFFF[i, 4] <- (cor(ModelOutputSPFFF2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF2[[i]]$ResultsVal$SVMReg))^2
  Accs2SPFFF[i, 1] <- sum(diag(table(ModelOutputSPFFF2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF2[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputSPFFF2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF2[[i]]$ResultsVal$GLM))
  Accs2SPFFF[i, 2] <- sum(diag(table(ModelOutputSPFFF2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF2[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputSPFFF2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF2[[i]]$ResultsVal$KerasClass))
  Accs2SPFFF[i, 3] <- sum(diag(table(ModelOutputSPFFF2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF2[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputSPFFF2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF2[[i]]$ResultsVal$SVMClass))

  Rsq1SPFFFScaled[i, 1] <- (cor(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$LM))^2
  Rsq1SPFFFScaled[i, 2] <- (cor(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$KM))^2
  Rsq1SPFFFScaled[i, 3] <- (cor(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$KerasReg))^2
  Rsq1SPFFFScaled[i, 4] <- (cor(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$SVMReg))^2
  Accs1SPFFFScaled[i, 1] <- sum(diag(table(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$GLM))
  Accs1SPFFFScaled[i, 2] <- sum(diag(table(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$KerasClass))
  Accs1SPFFFScaled[i, 3] <- sum(diag(table(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputSPFFFScaled1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled1[[i]]$ResultsVal$SVMClass))
  Rsq2SPFFFScaled[i, 1] <- (cor(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$LM))^2
  Rsq2SPFFFScaled[i, 2] <- (cor(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$KM))^2
  Rsq2SPFFFScaled[i, 3] <- (cor(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$KerasReg))^2
  Rsq2SPFFFScaled[i, 4] <- (cor(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YRegVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$SVMReg))^2
  Accs2SPFFFScaled[i, 1] <- sum(diag(table(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$GLM))
  Accs2SPFFFScaled[i, 2] <- sum(diag(table(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$KerasClass))
  Accs2SPFFFScaled[i, 3] <- sum(diag(table(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputSPFFFScaled2[[i]]$ResultsVal$YClassVal, ModelOutputSPFFFScaled2[[i]]$ResultsVal$SVMClass))
}


for(i in (1:length(DoEindicesFFF))){
  YsFFF1[[i]] <- Testbed1(FFFs[[DoEindicesFFF[i]]][,-1])   
  YsFFF2[[i]] <- Testbed2(FFFs[[DoEindicesFFF[i]]][,-1], beta_pf, beta_num)   
  ModelOutputFFF1[[i]] <- ModelPreds(DoE=FFFs[[DoEindicesFFF[i]]][,-1], WP=FFFs[[DoEindicesFFF[i]]][,1], YReg=YsFFF1[[i]][,2], YsFFF1[[i]][,1], DoEVal, YRegVal1, YClassVal1, nepochs = nepochs, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)
  ModelOutputFFF2[[i]] <- ModelPreds(DoE=FFFs[[DoEindicesFFF[i]]][,-1], WP=FFFs[[DoEindicesFFF[i]]][,1], YReg=YsFFF2[[i]][,2], YsFFF2[[i]][,1], DoEVal, YRegVal2, YClassVal2, nepochs = nepochs, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)

  Rsq1FFF[i, 1] <- (cor(ModelOutputFFF1[[i]]$ResultsVal$YRegVal, ModelOutputFFF1[[i]]$ResultsVal$LM))^2
  Rsq1FFF[i, 2] <- (cor(ModelOutputFFF1[[i]]$ResultsVal$YRegVal, ModelOutputFFF1[[i]]$ResultsVal$KM))^2
  Rsq1FFF[i, 3] <- (cor(ModelOutputFFF1[[i]]$ResultsVal$YRegVal, ModelOutputFFF1[[i]]$ResultsVal$KerasReg))^2
  Rsq1FFF[i, 4] <- (cor(ModelOutputFFF1[[i]]$ResultsVal$YRegVal, ModelOutputFFF1[[i]]$ResultsVal$SVMReg))^2
  Accs1FFF[i, 1] <- sum(diag(table(ModelOutputFFF1[[i]]$ResultsVal$YClassVal, ModelOutputFFF1[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputFFF1[[i]]$ResultsVal$YClassVal, ModelOutputFFF1[[i]]$ResultsVal$GLM))
  Accs1FFF[i, 2] <- sum(diag(table(ModelOutputFFF1[[i]]$ResultsVal$YClassVal, ModelOutputFFF1[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputFFF1[[i]]$ResultsVal$YClassVal, ModelOutputFFF1[[i]]$ResultsVal$KerasClass))
  Accs1FFF[i, 3] <- sum(diag(table(ModelOutputFFF1[[i]]$ResultsVal$YClassVal, ModelOutputFFF1[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputFFF1[[i]]$ResultsVal$YClassVal, ModelOutputFFF1[[i]]$ResultsVal$SVMClass))
  Rsq2FFF[i, 1] <- (cor(ModelOutputFFF2[[i]]$ResultsVal$YRegVal, ModelOutputFFF2[[i]]$ResultsVal$LM))^2
  Rsq2FFF[i, 2] <- (cor(ModelOutputFFF2[[i]]$ResultsVal$YRegVal, ModelOutputFFF2[[i]]$ResultsVal$KM))^2
  Rsq2FFF[i, 3] <- (cor(ModelOutputFFF2[[i]]$ResultsVal$YRegVal, ModelOutputFFF2[[i]]$ResultsVal$KerasReg))^2
  Rsq2FFF[i, 4] <- (cor(ModelOutputFFF2[[i]]$ResultsVal$YRegVal, ModelOutputFFF2[[i]]$ResultsVal$SVMReg))^2
  Accs2FFF[i, 1] <- sum(diag(table(ModelOutputFFF2[[i]]$ResultsVal$YClassVal, ModelOutputFFF2[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputFFF2[[i]]$ResultsVal$YClassVal, ModelOutputFFF2[[i]]$ResultsVal$GLM))
  Accs2FFF[i, 2] <- sum(diag(table(ModelOutputFFF2[[i]]$ResultsVal$YClassVal, ModelOutputFFF2[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputFFF2[[i]]$ResultsVal$YClassVal, ModelOutputFFF2[[i]]$ResultsVal$KerasClass))
  Accs2FFF[i, 3] <- sum(diag(table(ModelOutputFFF2[[i]]$ResultsVal$YClassVal, ModelOutputFFF2[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputFFF2[[i]]$ResultsVal$YClassVal, ModelOutputFFF2[[i]]$ResultsVal$SVMClass))
}

save.image("FittedModels_20200106.Rdata")

# For classification, i = 1: GLM, i = 2: NN, i = 3: SVM
# For Regression: i = 1: LM, i = 2: KM, i = 3: NN, i = 4: SVM 
#####################
# plot prediction results
#####################

Accs2SPFFFScaled
Accs2SPFFF
Accs2FFF
Rsq2SPFFFScaled
Rsq2SPFFF
Rsq2FFF
DoEindicesSPFFF
DoEindicesFFF
SPFFF_criteria
FFF_criteria

SPFFF_criteria[DoEindicesSPFFF, ]$nWP

SPFFF_criteria[DoEindicesSPFFF, ]$n
SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8
ClassType <- c("GLM", "NN", "SVM")
png("PredictionComp1Acc_20200106.png", width = 480*3/2, height = 480/2)
YLIM <- range(c(range(Accs1FFF), range(Accs1SPFFF), range(Accs1SPFFFScaled)))
YLIM[1] <- YLIM[1] *0.95
# For classification, i = 1: GLM, i = 2: NN, i = 3: SVM
# For Regression: i = 1: LM, i = 2: KM, i = 3: NN, i = 4: SVM 
par(mfrow = c(1,3), oma = c(0,0,2,0), mar = c(5, 4, 3, 2) + 0.1)
for(i in (1:3)){
  plot(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], xlab = "n", ylab = "Accuracy", ylim = YLIM, col = 1, pch = 15)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 1, pch = 15)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 1)
#  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i, pch = 16)
#  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 1, pch = 17)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 1)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 2, pch = 15)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 2)
#  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i, pch = 16)
#  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 2, pch = 17)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 2)
  points(FFF_criteria[DoEindicesFFF, ]$n,Accs1FFF[,i], col = 3, pch = 18)
  lines(FFF_criteria[DoEindicesFFF, ]$n,Accs1FFF[,i], col = 3)
  title(paste("Validation data set accuracy for", ClassType[i]))
}
title("Accuracy for Test bed 1 (Cantilever function), classification",outer = TRUE)
legend("bottomright", col = c(1, 1, 2, 2, 3), pch = c(15, 17, 15, 17, 18), legend = c("SPFFF8", "SPFFF16", "SPFFFs8", "SPFFFs16", "FFF"))
dev.off()

YLIM <- range(c(range(Accs2FFF), range(Accs2SPFFF), range(Accs2SPFFFScaled)))
plot(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs2SPFFF[,1][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], xlab = "n", ylab = "Accuracy", ylim = YLIM, col = 1)
for(i in (1:3)){
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8],]$n, Accs2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12],]$n, Accs2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16],]$n, Accs2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  points(FFF_criteria[DoEindicesFFF, ]$n,Accs2FFF[,i], col = i, pch = 15)
  lines(FFF_criteria[DoEindicesFFF, ]$n,Accs2FFF[,i], col = i)
}

YLIM <- range(c(range(Rsq1FFF, na.rm = TRUE), range(Rsq1SPFFF, na.rm = TRUE), range(Rsq1SPFFFScaled, na.rm = TRUE)))
YLIM <- c(0.8, 1)
png("PredictionComp1Rsq_20200106.png", width = 480*3/2, height = 480*3/2)
par(mfrow = c(2,2), oma = c(0,0,2,0), mar = c(5, 4, 3, 2) + 0.1)
RegType <- c("LM", "Kriging", "NN", "SVM")
# For classification, i = 1: GLM, i = 2: NN, i = 3: SVM
# For Regression: i = 1: LM, i = 2: KM, i = 3: NN, i = 4: SVM 
for(i in (1:4)){
#i <- 3
  plot(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], xlab = "n", ylab = "Rsquare", ylim = YLIM, pch = 15)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 1, pch = 15)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 1)
#  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
#  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 1, pch = 17)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq1SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 1)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 2, pch = 15)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = 2)
#  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], pch = 14, col = i)
#  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 2, pch = 17)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq1SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = 2)
  points(FFF_criteria[DoEindicesFFF, ]$n, Rsq1FFF[,i], col = 3, pch = 18)
  lines(FFF_criteria[DoEindicesFFF, ]$n, Rsq1FFF[,i], col = 3)
  title(paste("Validation data set Rsquare for", RegType[i]))
}
legend("bottomright", col = c(1, 1, 2, 2, 3), pch = c(15, 17, 15, 17, 18), legend = c("SPFFF8", "SPFFF16", "SPFFFs8", "SPFFFs16", "FFF"))
title("Rsquare for Test bed 1 (Cantilever function), regression", outer = TRUE)
dev.off()


YLIM <- range(c(range(Rsq2FFF, na.rm = TRUE), range(Rsq2SPFFF, na.rm = TRUE), range(Rsq2SPFFFScaled, na.rm = TRUE)))
plot(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq2SPFFF[,1][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], xlab = "n", ylab = "Rsquare", ylim = YLIM)
for(i in (1:4)){
  #i <- 3
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq2SPFFF[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], pch = 14, col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], ]$n, Rsq2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 8], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], pch = 14, col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], ]$n, Rsq2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 12], col = i)
  points(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], pch = 14, col = i)
  lines(SPFFF_criteria[DoEindicesSPFFF[SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], ]$n, Rsq2SPFFFScaled[,i][SPFFF_criteria[DoEindicesSPFFF, ]$nWP == 16], col = i)
  points(FFF_criteria[DoEindicesFFF, ]$n, Rsq2FFF[,i], pch = 15, col = i)
  lines(FFF_criteria[DoEindicesFFF, ]$n, Rsq2FFF[,i], col = i)
}



########################
# old
########################
i <- 1
plot(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$LM)
plot(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$KM)
plot(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasReg)
plot(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMReg)
(cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$LM))^2
(cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$KM))^2
(cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasReg))^2
(cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMReg))^2
table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$GLM)
table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasClass)
table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMClass)
sum(diag(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$GLM))) / sum(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$GLM))
sum(diag(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasClass))) / sum(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$KerasClass))
sum(diag(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMClass))) / sum(table(ModelOutputSPFFF1[[i]]$ResultsVal$YClassVal, ModelOutputSPFFF1[[i]]$ResultsVal$SVMClass))

str(ModelOutputSPFFF1[[1]]$ResultsVal$YClassVal, )

cor(ModelOutputSPFFF1[[i]]$ResultsVal$YRegVal, ModelOutputSPFFF1[[i]]$ResultsVal$LM)^2


par(mfrow = c(1,2))
plot(GLMrs, ylim = c(0, 1))
points(1:5, KerasClassrs, col = "red")
points(1:5, SVMClassrs, col = "blue")

plot(LMrs, ylim = c(0, 1))
points(1:5, KMrs, col = "red")
points(1:5, KerasRegrs, col = "blue")
points(1:5, SVMRegrs, col = "green")






#########################
### old
########################
### unify all designs to the same design space: by definition her [-1, 1]

apply(DoE1, 2, summary)
apply(DoE2, 2, summary)
apply(DoELHS, 2, summary)
apply(DoEFFF, 2, summary)
apply(DoESP, 2, summary)
### trying out acebayes functions:
### acebayes is unstable, thowing out an error
#start.d<-matrix(2 * randomLHS(n = noverall,k = d) - 1,nrow = noverall,ncol = d,
#                dimnames = list(as.character(1:noverall), c("x1", "x2", "x3", "x4")))
## Generate an initial design of appropriate dimension. The initial design is a 
## Latin hypercube sample.
#start.d
#low<-c(-3, 4, 5, -6, -2.5)
#upp<-c(3, 10, 11, 0, 3.5)
## Lower and upper limits of the uniform prior distributions.

#prior <- function(B){t(t(6*matrix(runif(n = 5 * B),ncol = 5)) + low)}
## Create a function which specifies the prior. This function will return a 
## B by 5 matrix where each row gives a value generated from the prior 
## distribution for the model parameters.
#example1 <- aceglm(formula=~x1+x2+x3+x4, start.d = start.d, family = binomial, prior = prior, method = "MC", N1 = 1, N2 = 0, B = c(1000, 1000))


########################
### calculate test bed responses based on the constructed DoE's
########################
nval <- 5000
DoEVal <- data.frame(matrix(runif(nval*d), ncol = d) * 2 - 1)
str(DoEVal)
str(DoE1)
Y1.Val <- Testbed1(DoEVal)
Y1.1 <- Testbed1(DoE1[,-1])
Y1.2 <- Testbed1(DoE2[,-1])
Y1.LHS <- Testbed1(DoELHS[,-1])
Y1.FFF <- Testbed1(DoEFFF[,-1])
Y1.SP <- Testbed1(DoESP[,-1]) 

set.seed(7)
beta_pf <- round(runif(2*d + 1 + choose(d,2))*4-2)
beta_num <- round(runif(2*d + 1 + choose(d,2))*4-2)
beta_pf
beta_num

Y2.Val <- Testbed2(DoEVal, beta_pf, beta_num)
Y2.1 <- Testbed2(DoE1[,-1], beta_pf, beta_num)
Y2.2 <- Testbed2(DoE2[,-1], beta_pf, beta_num)
Y2.LHS <- Testbed2(DoELHS[,-1], beta_pf, beta_num)
Y2.FFF <- Testbed2(DoEFFF[,-1], beta_pf, beta_num)
Y2.SP <- Testbed2(DoESP[,-1], beta_pf, beta_num)

########################
### set up models
########################


########################
### set up models all together:
########################

#DoE <- DoE1[,-1]
# DoE is DoE without WP
#WP <- DoE1[,1]
#YReg <- Y1.1[,2]
#YClass <- Y1.1[,1]
#DoEVal
#source(file = "FFFSplitPLotFunctions_20190812.R")

#######################
# Evaluate designs with 
#######################


#######################
# Evaluate designs  
#######################
### still to do: 
### - get the bayesian D-optimal design criteria from acebayes
### - Comparison D efficiency to a JMP DoE
### - D-optimal design
### - potentially also using the package AlgDesign
DoEs <- list(DoE1 = DoE1, DoE2 = DoE2, DoESP = DoESP, DoELHS = DoELHS, DoEFFF = DoEFFF)
DoEnames <- c('DoE1', 'DoE2', 'SP', 'LHS', 'FFF')
nmM <- 50000
IoptimalsWP <- Ioptimals <- PhiMMs <- DoptimalsWP <- Doptimals <- MaxiMinsUnique <- Maximins <- Minimaxs <- numeric(5)
Minimaxs2 <- matrix(0, nrow = 4, ncol = 5)
for(i in 1:5){
  print(i)
  Minimaxs[i] <- Minimax(DoEs[[i]][,-1], nsim = 1000, seed = 1)
  Maximins[i] <- min(dist(DoEs[[i]][,-1]))
  MaxiMinsUnique[i] <- min(dist(unique(DoEs[[i]][,-1])))
  Doptimals[i] <- Doptimal(DoEs[[i]][,-1])
  PhiMMs[i] <- sum(1/dist(DoEs[[i]][,-1])^2)
  Ioptimals[i] <- Ioptimal(DoEs[[i]][,-1])
  Minimaxs2[, i] <- Minimax2(DoEs[[i]][,-1], nsim = 50000)
  if(i < 4){
    DoptimalsWP[i] <- DoptimalWP(DoEs[[i]][,-1], DoEs[[i]][,1])
    IoptimalsWP[i] <- IoptimalWP(DoEs[[i]][,-1],DoEs[[i]][,1])
  }
}

Minimaxs2

DoEnames
OptimizationDirection <- c("na", "min", "min", "min", "max", "max", "max", "max", "min")
OptimalityCriteria <- data.frame(DoEnames, Ioptimals, IoptimalsWP, PhiMMs, Doptimals, DoptimalsWP, MaxiMinsUnique, Maximins, Minimaxs)
OptimalityCriteria
t(Minimaxs2)
#min(dist(unique(DoE1[,2:(dWP+ 1 )])))
#min(dist(unique(DoE2[,2:(dWP+ 1 )])))

#min(dist(DoE1[,(dWP + 2):(dWP + dSP + 1)]))
#min(dist(DoE2[,(dWP + 2):(dWP + dSP + 1)]))

############################
# Compare prediction power
############################

YRegVal <- Y1.Val[,2]
YClassVal <- Y1.Val[,1]  

Y1s <- list(Y1.1, Y1.2, Y1.SP, Y1.LHS, Y1.FFF)
Outputs <- vector("list", 5)
for(l in (1:5)){
  print(l)
  Outputs[[l]] <- ModelPreds(DoE=DoEs[[l]][,-1], WP=DoEs[[l]][,1], YReg=Y1s[[l]][,2], Y1s[[l]][,1], DoEVal, YRegVal, YClassVal, nepochs = 10, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)  
}  

str(Outputs[[l]][[1]])

ndesign <- length(DoEs)
GLMrs <- numeric(ndesign)
KerasClassrs <- numeric(ndesign)
SVMClassrs <- numeric(ndesign)
LMrs <- numeric(ndesign)
KMrs <- numeric(ndesign)
KerasRegrs <- numeric(ndesign)
SVMRegrs <- numeric(ndesign)

for(l in 1:5){
  GLMrs[l] <- sum(diag(table(Outputs[[l]][[1]]$GLM, Outputs[[l]][[1]]$YClassVal)))/sum(table(Outputs[[l]][[1]]$GLM, Outputs[[l]][[1]]$YClassVal))
  KerasClassrs[l] <- sum(diag(table(Outputs[[l]][[1]]$KerasClass, Outputs[[l]][[1]]$YClassVal)))/sum(table(Outputs[[l]][[1]]$KerasClass, Outputs[[l]][[1]]$YClassVal))
  SVMClassrs[l] <- sum(diag(table(Outputs[[l]][[1]]$SVMClass, Outputs[[l]][[1]]$YClassVal)))/sum(table(Outputs[[l]][[1]]$SVMClass, Outputs[[l]][[1]]$YClassVal))
  LMrs[l] <- (cor(Outputs[[l]][[1]]$LM, Outputs[[l]][[1]]$YRegVal))^2
  KMrs[l] <- (cor(Outputs[[l]][[1]]$KM, Outputs[[l]][[1]]$YRegVal))^2
  KerasRegrs[l] <- (cor(Outputs[[l]][[1]]$KerasReg, Outputs[[l]][[1]]$YRegVal))^2
  SVMRegrs[l] <- (cor(Outputs[[l]][[1]]$SVMReg, Outputs[[l]][[1]]$YRegVal))^2
}

par(mfrow = c(1,2))
plot(GLMrs, ylim = c(0, 1))
points(1:5, KerasClassrs, col = "red")
points(1:5, SVMClassrs, col = "blue")

plot(LMrs, ylim = c(0, 1))
points(1:5, KMrs, col = "red")
points(1:5, KerasRegrs, col = "blue")
points(1:5, SVMRegrs, col = "green")

#####################
### second examples:
#####################
YRegVal <- Y2.Val[,2]
YClassVal <- Y2.Val[,1]  
DoEs <- list(DoE1 = DoE1, DoE2 = DoE2, DoESP = DoESP, DoELHS = DoELHS, DoEFFF = DoEFFF)
str(DoEs)
Y2s <- list(Y2.1, Y2.2, Y2.SP, Y2.LHS, Y2.FFF)
Outputs <- vector("list", 5)
for(l in (1:5)){
  print(l)
  Outputs[[l]] <- ModelPreds(DoE=DoEs[[l]][,-1], WP=DoEs[[l]][,1], YReg=Y2s[[l]][,2], Y2s[[l]][,1], DoEVal, YRegVal, YClassVal, nepochs = 10, nugget = 0.0001, costSVMReg = 1000, gammaSVMReg = 0.0001, costSVMClass = 100, gammaSVMclass = 1)  
}  

str(Outputs[[l]][[1]])

ndesign <- length(DoEs)
GLMrs <- numeric(ndesign)
KerasClassrs <- numeric(ndesign)
SVMClassrs <- numeric(ndesign)
LMrs <- numeric(ndesign)
KMrs <- numeric(ndesign)
KerasRegrs <- numeric(ndesign)
SVMRegrs <- numeric(ndesign)

for(l in 1:5){
  GLMrs[l] <- sum(diag(table(Outputs[[l]][[1]]$GLM, Outputs[[l]][[1]]$YClassVal)))/sum(table(Outputs[[l]][[1]]$GLM, Outputs[[l]][[1]]$YClassVal))
  KerasClassrs[l] <- sum(diag(table(Outputs[[l]][[1]]$KerasClass, Outputs[[l]][[1]]$YClassVal)))/sum(table(Outputs[[l]][[1]]$KerasClass, Outputs[[l]][[1]]$YClassVal))
  SVMClassrs[l] <- sum(diag(table(Outputs[[l]][[1]]$SVMClass, Outputs[[l]][[1]]$YClassVal)))/sum(table(Outputs[[l]][[1]]$SVMClass, Outputs[[l]][[1]]$YClassVal))
  LMrs[l] <- (cor(Outputs[[l]][[1]]$LM, Outputs[[l]][[1]]$YRegVal))^2
  KMrs[l] <- (cor(Outputs[[l]][[1]]$KM, Outputs[[l]][[1]]$YRegVal))^2
  KerasRegrs[l] <- (cor(Outputs[[l]][[1]]$KerasReg, Outputs[[l]][[1]]$YRegVal))^2
  SVMRegrs[l] <- (cor(Outputs[[l]][[1]]$SVMReg, Outputs[[l]][[1]]$YRegVal))^2
}

par(mfrow = c(1,2))
plot(GLMrs, ylim = c(0, 1))
points(1:5, KerasClassrs, col = "red")
points(1:5, SVMClassrs, col = "blue")

plot(LMrs, ylim = c(0, 1))
points(1:5, KMrs, col = "red")
points(1:5, KerasRegrs, col = "blue")
points(1:5, SVMRegrs, col = "green")


