# Step 0. Install necessary packages. 
## If you already installed all, then skip the Step 0. 

# install.packages('car')
# install.packages('MASS')
# install.packages('SuperLearner')
# install.packages(c("caret", "glmnet", "randomForest", "RhpcBLASctl", "xgboost", "gam"))

library(ggplot2) 
library(SuperLearner)
require(stats); require(graphics)
library(splines)

# Step 1. Recall the data. Please refer the comment on data_generation.R 

## This generate 'data' (observations)
### 'data' is observational (Y,X,Z)
sepData = function(DF){
  Yobs = DF[,1] # First column of 'data' is Y. 
  X = DF[,2] # X variable in 'data'.
  Z = DF[,-c(1,2)] # Z variables in 'data'. 
  XZ = data.frame(cbind(X,Z))
  return(list(Yobs,X,Z,XZ))
}

trainSL = function(Y,X,family, numFold){
  if (family == 'cont'){
    train_model = SuperLearner(Y = Y, X = X, family = gaussian(), cvControl = list(V=numFold),
                               SL.library = c("SL.xgboost","SL.lm","SL.randomForest"))  
  }
  else if (family == 'disc'){
    train_model = SuperLearner(Y = Y, X = X, family = binomial(), cvControl = list(V=numFold),
                               SL.library = c("SL.xgboost","SL.lm","SL.randomForest"))  
  }
  return(train_model)
}

sampleEst = function(trainModel, XZ){
  predResult = predict(trainModel, XZ, onlySL=T)$pred
  return(predResult)
}

reducedX = function(X, cutN = 100){
  cutN = min(N/10, cutN)
  sample_X = sample(X,cutN)
  return(sample_X)
}

sampleIntvEst = function(trainModel, sample_X, Z, cutN = 100){
  Yhat_intv = rep(0, cutN)
  for (idx in 1:cutN){
    xi = sample_X[idx]
    Zi = Z[sample(N,cutN),]
    xi_mat = matrix(rep(xi,cutN),nrow=cutN) # N length of 0 vector
    XZi = cbind(xi_mat, Zi)
    colnames(XZi)[1] = 'X'
    Yhat_intv[idx] = mean(sampleEst(g_xz, XZi))
  }
  return(Yhat_intv)
}

sampleObsEst = function(trainModel, sample_X, cutN = 100){
  sample_X_input = data.frame(sample_X)
  sample_X_input = data.frame(sample_X)
  colnames(sample_X_input)[1] = 'X'
  Yhat_obs = sampleEst(g_x,sample_X_input)
  return(Yhat_obs)
  
}


resampler <- function(data) {
  n <- nrow(data)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(data[resample.rows,])
}

spline_estimator <- function(data,n) {
  fit = smooth.spline(x=data[,1],y=data[,2],cv=TRUE)
  eval_grid = seq(from=min(data[,1]),to=max(data[,1]),length.out=n)
  return(predict(fit,x=eval_grid)$y) # We only want the predicted values
}

spline_CI <- function(data,n,B,alpha) {
  spline_main = spline_estimator(data,n)
  spline_boots <- replicate(B,spline_estimator(resampler(data),n))
  cis_lower = 2*spline_main - apply(spline_boots,1,quantile,probs=1-alpha/2)
  cis_upper <- 2*spline_main - apply(spline_boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline_main,lower.ci=cis_lower,upper.ci=cis_upper,
              x=seq(from=min(data[,1]),to=max(data[,1]),length.out=n)))
}

plotResult = function(X,Y, N, B, alpha, title){
  set.seed(1)
  df1 = data.frame(X, Y)
  result_spline = spline_CI(df1,N, B=B,alpha=alpha)
  df2 = data.frame( x=result_spline$x, y=result_spline$main.curve   )
  df3 = data.frame( x=result_spline$x, ymin=result_spline$lower.ci, ymax=result_spline$upper.ci   )
  
  gg1 = ggplot(df1, aes(X, Y)) + geom_point(size=1.5)
  gg1 = gg1 + geom_line(data=df2,aes(x=x, y=y))
  gg1 = gg1 + geom_ribbon(data=df3, aes(x=x, ymin=ymin,ymax=ymax),alpha=0.4)
  
  gg1 = gg1 + ggtitle(title)
  gg1 = gg1 + theme_bw()
  gg1 = gg1 + scale_x_continuous(name = "X=x") + scale_y_continuous(name = "Estimated E[Y|do(x)]",limits=c(min(Y), A=max(Y)))
  gg1 = gg1 + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    text=element_text(size = 16))
  return(gg1)
}

################## MAIN ##################

source('datagen/data_generation_cano1_cont.R') 

result_sep = sepData(data_obs)
Yobs = result_sep[[1]]
X = result_sep[[2]]
Z = result_sep[[3]]
XZ = result_sep[[4]]
N = length(Yobs)

numFold = 3
n = 50

g_xz = trainSL(c(Yobs),XZ,'cont',numFold) # E[Y|x,z]
g_x = trainSL(c(Yobs),data.frame(X),'cont',numFold)

sample_X = reducedX(X,cutN = n)
Yhat_intv = sampleIntvEst(g_xz, sample_X, Z, cutN = n)
Yhat_obs = sampleObsEst(g_x, sample_X, cutN = n)

B = 300
alpha = 0.01

gg_intv = plotResult(sample_X,Yhat_intv,n,B,alpha,"Interventional")
gg_obs = plotResult(sample_X, Yhat_obs,n,B,alpha,"Observational")
library("gridExtra")
grid.arrange(arrangeGrob(gg_intv,gg_obs,nrow=1))

# 
# 
# 
# df1 = data.frame(sample_X, Yhat_obs)
# result_spline = spline_CI(df1, n, B=B,alpha=alpha)
# df2 = data.frame( x=result_spline$x, y=result_spline$main.curve   )
# df3 = data.frame( x=result_spline$x, ymin=result_spline$lower.ci, ymax=result_spline$upper.ci   )
# 
# gg1 = ggplot(df1, aes(sample_X, Yhat_obs)) + geom_point(size=1.5)
# gg1 = gg1 + geom_line(data=df2,aes(x=x, y=y))
# gg1 = gg1 + geom_ribbon(data=df3, aes(x=x, ymin=ymin,ymax=ymax),alpha=0.4)
# 
# gg1 = gg1 + ggtitle(title)
# gg1 = gg1 + theme_bw()
# gg1 = gg1 + scale_x_continuous(name = "X=x") + scale_y_continuous(name = "Estimated E[Y|do(x)]",limits=c(min(Y), A=max(Y)))
# gg1 = gg1 + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                   axis.line.y = element_line(size = 0.5, colour = "black"),
#                   axis.line = element_line(size=1, colour = "black"),
#                   panel.border = element_blank(),
#                   panel.background = element_blank(),
#                   plot.title=element_text(size = 20),
#                   text=element_text(size = 16))
# PLOT 
## Y_intv 




