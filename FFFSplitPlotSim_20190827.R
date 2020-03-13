#?source
source(file = "FFFSplitPLotFunctions_20190812.R")

dWP <- 1
nWP <- 10
dSP <- 3
noverall <- 30
nsim <- 200
seed <- 1002
printsteps <- 50
DoE1 <- SplitPlotFFFDesign(dWP = dWP , nWP = nWP, dSP = dSP , noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)


min(dist(DoE1[,-1]))
min(dist(unique(DoE1[,2:(dWP+ 1 )])))
min(dist(DoE1[,(dWP + 2):(dWP + dSP + 1)]))

### potentielle Verbesserungen: Auch kategorialle Variable, Constraint design spaces....  

dWP <- 2
nWP <- 20
dSP <- 2
noverall <- 100
nsim <- 500
seed <- 1003
printsteps <- 5

DoE1 <- SplitPlotFFFDesign(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)
DoE2 <- SplitPlotFFFDesign2(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)

##########################
### compare designs:
##########################
Minimax(DoE1[,-1], nsim = 5000, seed = 1)
Minimax(DoE2[,-1], nsim = 5000, seed = 1)

min(dist(DoE1[,-1]))
min(dist(DoE2[,-1]))

min(dist(unique(DoE1[,2:(dWP+ 1 )])))
min(dist(unique(DoE2[,2:(dWP+ 1 )])))

min(dist(DoE1[,(dWP + 2):(dWP + dSP + 1)]))
min(dist(DoE2[,(dWP + 2):(dWP + dSP + 1)]))

sum(1/dist(DoE1[,-1])^2)
sum(1/dist(DoE2[,-1])^2)

pairs(DoE1[,-1], col = DoE1[,1])
pairs(DoE2[,-1], col = DoE2[,1])

### Comparison D efficiency to a JMP DoE
#######################
### Ein Model aufstellen
#######################

### zuerst ein normales GP Model, juhuu!;-)
#library(geoR)
#library(geoRglm)
### use branin*branin as example:

dWP <- 2
nWP <- 20
dSP <- 2
noverall <- 100
nsim <- 500
seed <- 1003
printsteps <- 5
### use DoE's from previos step

DoE1 <- SplitPlotFFFDesign(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)
DoE2 <- SplitPlotFFFDesign2(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)

source(file = "FFFSplitPLotFunctions_20190812.R")

DoE <- DoE2[,-1]
WP <- DoE2[,1]
nval <- 500
d <- dWP + dSP 
xval <- data.frame(matrix(runif(nval*d), ncol = d))
colnames(xval) <- colnames(DoE)

#testbed1:
Testbed1(DoE)
Y <- Testbed1(DoE)
Yval <- Testbed1(xval)
str(Yval)

#str(y1_1)
#table(round(y1_1))
plot(Y$y1, Y$y2, log = "y")
pairs(DoE, pch = Y$y1 + 1, cex = 2, col = Y$y2)


#i <- i + 1
set.seed(7)
beta_pf <- round(runif(2*d + 1 + choose(d,2))*4-2)
beta_num <- round(runif(2*d + 1 + choose(d,2))*4-2)
beta_pf
beta_num

TB2 <- Testbed2(DoE, beta_pf, beta_num)
str(TB2)
plot(TB2)

#KM <- km(~., design = data.frame(DoE), response = y, covtype = "gauss")
#plot(KM)

##############################
### construct competitor designs: 
### one from acebayes package for logistic regression
### one split plot RSM from JMP
### one FFF design from JMP
### one maximinLHS from package lhs
### one doe from package acebayes

names(DoE)
LHS <- maximinLHS(noverall, d, method = "iterative", maxIter = 100)
?read.table
FFF <- read.table("Fast Flexible Filling Design_4factors_n30_20190827.csv", sep = ",", header = TRUE, dec = ",")
FFF <- FFF[,1:4]
names(FFF) <- names(DoE)
str(FFF)
SPDoE <- read.table("RSM2etcfactors2htcFactors_10WPn30_20190827.csv", sep = ",", header = TRUE, dec = ",")
SPDoE <- SPDoE[,1:5]
B
utilfisher <- function(d, B) {
theta <- rnorm(B)
ui <- matrix(rep(d[, 1] ^ 2, B), ncol = B) * exp(outer(d[, 1], theta))
apply(ui, 2, sum)
}
set.seed(1)
start.d <- matrix(0, nrow = noverall, ncol = d)
start.d <- matrix(0, nrow = 20, ncol = 2)
start.d <- optimumLHS(n = noverall, k = d)
start.d <- optimumLHS(n = 20, k = 2)
pairs(start.d)
ex22a <- ace(utility = utilfisher, start.d = start.d)
ex22a <- ace(utility = utilfisher, start.d = start.d, Q = 5, N1 = 4, N2 = 10)
ex22a$start.d
str(ex22a$phase1.d)
pairs(ex22a$phase1.d)
pairs(ex22a$phase2.d)
ex22a$phase2.d
#################################
### construct a dabble data set
dWP <- 2
nWP <- 6
dSP <- 2
noverall <- 18
nsim <- 500
seed <- 1002
printsteps <- 5
DoE <- SplitPlotFFFDesign2(dWP = dWP, nWP = nWP, dSP = dSP, noverall = noverall, nsim = nsim, seed = seed, printsteps = printsteps)

WP <- DoE[,1]
DoE <- DoE[,-1]
str(DoE)
y <- as.matrix(branin(DoE[,c(1,2)])*branin(DoE[,c(3,4)]))

### implement a likelihood function:
theta <- 1.2
sigma_WP <- 0.4
sigma <- 1
sigma_Spatial <- 0.6
COR <- "gaussk"
n <- nrow(DoE)



Paramintial <- runif(d+1 + d+ 3)
#Parameters <- list(beta, theta, sigma_Spatial, sigma_WP, sigma)
LLHmanual(Paramintial, DoE, WP, y,COR)
### alternativer for llh:
### Optimizing Likelihood:
n_initial <- 100
betahatOLS <- as.numeric(lm(as.numeric(y)~as.matrix(DoE))$coefficients)
ParamintialMatrix <- cbind(matrix(betahatOLS, ncol = d + 1, nrow = n_initial, byrow = TRUE), matrix(runif((d+3)*n_initial), nrow = n_initial) + 0.01)

#ParamintialMatrix <- matrix(runif((d+1+d+3)*n_initial), nrow = n_initial) + 0.01
#LLstartValues <- apply(ParamintialMatrix, 1, LLH2, DoE = DoE, WP = WP, y = y, COR = COR)
str(ParamintialMatrix)
LLstartValues <- apply(ParamintialMatrix, 1, LLHmanual, DoE = DoE, WP = WP, y = y, COR = COR)

imax <- which.max(LLstartValues)

thetastart <- ParamintialMatrix[imax,]
lower <- c(rep(-100000, d + 1), rep(0.0001, d + 3) )
upper <- c(rep(100000, d + 1), rep(14000, d), rep(100000, 3) )
Optim <- optim(thetastart, LLHmanual, DoE = DoE, y = y, WP = WP, COR = COR, method = "L-BFGS-B", control = list(fnscale = -1, trace = 2, REPORT = 5, maxit = 1000), lower = lower, upper = upper)


print(Optim$convergence)
ParOptim <- Optim$par
ParOptim
ParTemp <- list(ParOptim[1:(d+1)], ParOptim[(d+ 2):(2*d+1)], ParOptim[2*d+2], ParOptim[2*d+3], ParOptim[2*d+4])
beta <- ParamTemp[[1]]
theta <- ParamTemp[[2]]
sigma_Spatial <- ParamTemp[[3]]
sigma_WP <- ParamTemp[[4]]
sigma <- ParamTemp[[5]]
beta
theta

### compare to km

KM <- km(~., design = data.frame(DoE), response = y, covtype = "gauss")
plot(KM)
KM
#library(spatstat) # crossdist

### building prediction:
npred <- 500
DoEVal <- data.frame(matrix(runif(npred*(dWP + dSP)), ncol = dWP + dSP))
yval <- as.matrix(branin(DoEVal[,c(1,2)])*branin(DoEVal[,c(3,4)]))

CorSpatial <- matrix(1, ncol = n, nrow = n)
CorSpatial
for(k in 1:d){
  H <- as.matrix(dist(DoE[,k]))
  CorSpatial <- CorSpatial* Rfunc1(H, theta[k], COR)
}
WPF <- as.factor(WP)
MM <- model.matrix(~WPF-1)
CovWP <- sigma_WP * MM%*% t(MM)
R <- CorSpatial*sigma_Spatial + CovWP + sigma*diag(n)
CholR <- chol(R)
str(CholR)
Rinvs <- solve(R)
trend <- model.matrix(~., data = DoE)
betadach <- solve(t(trend) %*% Rinvs %*% trend) %*% t(trend) %*% Rinvs %*% y
Faktor2 <- Rinvs %*% (y - trend %*% betadach)
p <- ncol(trend)
sigmadachhoch2 <- (1 / (n - p)) * t((y - trend %*% betadach)) %*% Rinvs %*% (y - trend %*% betadach)
ydach <- numeric(npred)
sigmadachx <- numeric(npred)
CorSpatialPred <- matrix(1, ncol = npred, nrow = n)
CorSpatialPred
for(k in 1:d){
  H <- as.matrix(dist(c( DoEVal[,k], DoE[,k])))[(npred + 1):(n + npred), 1:npred]
  print(theta[k])
  print(Rfunc1(H, theta[k], COR))
  CorSpatialPred <- CorSpatialPred * Rfunc1(H, theta[k], COR)
}
trendpred <- model.matrix(~., data = data.frame(DoEVal))
predpart1 <- trendpred %*% betadach
CholR_t_inv <- solve(t(CholR))
v <- CholR_t_inv  %*% CorSpatialPred
z <- CholR_t_inv %*% (y - trend %*%betadach)
#%*% 
str(CorSpatialPred)
ypred <- predpart1 + t(v)%*%z

plot(ypred, yval)
abline(0,1)
