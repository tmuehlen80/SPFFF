### assumes, that the functions in "SPFFF_Algo_Git.R" and "SPFFF_utilityFunctions_Git.R" are loaded

setwd("~/08_Learning/Split Plot Fast Flexible Filling Design")

source(file = "SPFFF_Algo_Git.R")
source(file = "SPFFF_utilityFunctions_Git.R")

### compute a SPFFF design:
dWP <- 2
dSP <- 2
n <- 60
nWP <- 10
nsim <- n*10
time1 <- Sys.time()
SPFFF <- SplitPlotFFFDesign(dWP = dWP, nWP = nWP, dSP = dSP, noverall = n, nsim = nsim, seed = 450, printsteps = 10)
time2 <- Sys.time()
time2 - time1
SPFFF

Mins <- apply(SPFFF[,-1], 2, min)
Maxs <- apply(SPFFF[,-1], 2, max)

SPFFFscaled <- data.frame(ClusterassignmentWP = SPFFF[,1],(SPFFF[,-1] - Mins) * 2/ (Maxs - Mins) - 1)

apply(SPFFF, 2, min)
apply(SPFFFscaled, 2, min)
apply(SPFFF, 2, max)
apply(SPFFFscaled, 2, max)

### create LHS :

LHS <- data.frame(maximinLHS(n, dWP + dSP, method = "iterative", maxIter = 10000) * 2 - 1)

### Evaluate:

Minimax(SPFFF[,-1], nsim = 100000)[5]
min(dist(unique(SPFFF[,-1]))) # Maximin
sum(1/dist(SPFFF[,-1])^2) # Phi_2
Ioptimal(SPFFF[,-1])

Minimax(SPFFFscaled[,-1], nsim = 100000)[5]
min(dist(unique(SPFFFscaled[,-1]))) # Maximin
sum(1/dist(SPFFFscaled[,-1])^2) # Phi_2
Ioptimal(SPFFFscaled[,-1])

Minimax(LHS, nsim = 100000)[5]
min(dist(unique(LHS))) # Maximin
sum(1/dist(LHS)^2) # Phi_2
Ioptimal(LHS)

### create prediction for a test function:

DoEVal <- data.frame(matrix(runif((dWP + dSP)*10000) * 2 - 1, ncol = dWP + dSP))
YVal <- Testbed1(DoEVal)

YSPFFF <- Testbed1(SPFFF[,-1])   
YSPFFFscaled <- Testbed1(SPFFFscaled[,-1])   
YLHS <- Testbed1(LHS)   


### fit a model:
DataReg <- data.frame(Y = YSPFFF[,2], SPFFF[,-1])
DataClass <- data.frame(Y=as.factor(YSPFFF[,1]), SPFFF[,-1])

svm.modelReg <- svm( Y ~ ., data = DataReg, cost = 1000, gamma = 0.0001)
svm.modelClass <- svm( Y ~ ., data = DataClass, cost = 100, gamma = 1)

PredValClassSVM <- predict(svm.modelClass, DoEVal)
PredValRegSVM <- predict(svm.modelReg, DoEVal)

plot(PredValRegSVM, YVal[,2])
table(PredValClassSVM, YVal[,1])

