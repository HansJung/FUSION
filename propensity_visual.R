# Step 0. Install necessary packages. 
library(ggplot2) 
library(SuperLearner)
require(stats); require(graphics)
library(splines)
library(cowplot)


# Seperating data frame into Y,X,Z,XZ
## Small noise are added for smoothing over the discrete value 
noiseAdd = function(X,eps=1e-8){
  disc.criterion = 10 
  X = data.frame(X)
  if ((length(unique(X))) > disc.criterion ){
    return(X)
  }
  else{
    N = length(X)
    X = X + eps*rnorm(N)
    X = data.frame(X)
    return(X)  
  }
}

setColname = function(X, myColName){
  colnames(X) = myColName
  return(X)
}

sepData = function(DF){
  Yobs = DF[,1] # First column of 'data' is Y. 
  X = DF[,2] # X variable in 'data'.
  Z = DF[,-c(1,2)] # Z variables in 'data'. 
  XZ = data.frame(cbind(X,Z))
  return(list(Yobs,X,Z,XZ))
}

noiseSepData = function(DF,eps=1e-8){
  Yobs = DF[,1]
  X = DF[,2]
  Z = DF[,-c(1,2)]
  X.noise = noiseAdd(X,eps)
  Z.noise = noiseAdd(Z,eps)
  XZ.noise = cbind(X.noise, Z.noise)
  return(list(Yobs,X.noise,Z.noise,XZ.noise))
}

intvData = function(X,eps,do.val){
  if (do.val == 0){
    X = matrix(rep(0,N),nrow=N) # N length of 0 vector
    X.noise = noiseAdd(X,eps)
  }
  else if (do.val==1){
    X = matrix(rep(1,N),nrow=N) # N length of 0 vector
    X.noise = noiseAdd(X,eps)
  }
  return(X)
}

intvSepData = function(X,Z,eps,do.val){
  X = intvData(X,eps,do.val)
  Z = noiseAdd(Z,eps)
  return(cbind(X,Z))
}

# Training SuperLearner
trainSL = function(Y,X,family, numFold){
  if (family == 'cont'){
    train_model = SuperLearner(Y = Y, X = X, family = gaussian(), cvControl = list(V=numFold),
                               SL.library = c("SL.xgboost","SL.gam","SL.glm"))  
  }
  else if (family == 'disc'){
    train_model = SuperLearner(Y = Y, X = X, family = binomial(), cvControl = list(V=numFold),
                               SL.library = c("SL.xgboost","SL.gam","SL.glm"))  
  }
  return(train_model)
}

propensityScore = function(ps_xz,Z,eps, do.val){
  if (do.val == 0){
    prop_score = (1-predict(ps_xz, Z, onlySL = T)$pred) + eps # P(X=0|Z)  
  }
  else if (do.val==1){
    prop_score = predict(ps_xz, Z, onlySL = T)$pred + eps # P(X=1|Z)  
  }
  return(prop_score)
}

################## MAIN ##################
# Data generation 
# source('datagen/data_generation_cano1_cont.R')
# source('datagen/data_generation_simpson_cont.R') 
# source('datagen/data_generation_simpson_disc.R') # Load data_obs 
source('datagen/data_generation_M.R') # Load data_obs 
eps = 1e-8
resultSep = sepData(data)
Yobs = resultSep[[1]]

ylim_min = min(Yobs)
ylim_max = max(Yobs)

X = resultSep[[2]]
Z = resultSep[[3]]
XZ = resultSep[[4]]

resultNoiseSep = noiseSepData(data,eps)
X.noise = resultNoiseSep[[2]]
Z.noise = resultNoiseSep[[3]]
XZ.noise = resultNoiseSep[[4]]

# X0Z.noise = intvSepData(data$X,data$Z,eps,do.val=0)
# X1Z.noise = intvSepData(data$X,data$Z,eps,do.val=1)

X0Z.noise = intvSepData(X,Z,eps,do.val=0)
X1Z.noise = intvSepData(X,Z,eps,do.val=1)

numFold = 3
ps_xz = trainSL(c(X),X=Z.noise,'cont',numFold) # P(X|z)
propScore_X1 = propensityScore(ps_xz,Z.noise,eps, do.val=1)

qplot(Z.noise$X3,propScore_X1)


