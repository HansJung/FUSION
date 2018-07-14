# REFERENCE 
## https://stackoverflow.com/questions/23852505/how-to-get-confidence-interval-for-smooth-spline
## http://www.stat.cmu.edu/~cshalizi/402/lectures/11-splines/lecture-11.pdf

# Step 0. Install necessary packages. 
library(ggplot2) 
library(SuperLearner)
require(stats); require(graphics)
library(splines)
library(cowplot)

# Seperating data frame into Y,X,Z,XZ
sepData = function(DF){
  Yobs = DF[,1] # First column of 'data' is Y. 
  X = DF[,2] # X variable in 'data'.
  Z = DF[,-c(1,2)] # Z variables in 'data'. 
  XZ = data.frame(cbind(X,Z))
  return(list(Yobs,X,Z,XZ))
}

# Training SuperLearner
trainSL = function(Y,X,family, numFold){
  if (family == 'cont'){
    train_model = SuperLearner(Y = Y, X = X, family = gaussian(), cvControl = list(V=numFold),
                               SL.library = c("SL.lm","SL.randomForest", "SL.gam","SL.xgboost"))  
  }
  else if (family == 'disc'){
    train_model = SuperLearner(Y = Y, X = X, family = binomial(), cvControl = list(V=numFold),
                               SL.library = c("SL.lm","SL.randomForest","SL.gam","SL.xgboost"))  
  }
  return(train_model)
}

# Predict using the trained SL model given XZ 
sampleEst = function(trainModel, XZ){
  predResult = predict(trainModel, XZ, onlySL=T)$pred
  return(predResult)
}

# Random sample from X (for reducing the size of X)
reducedX = function(X, cutN = 100){
  cutN = min(N/10, cutN)
  sample_X = sample(X,cutN)
  return(sample_X)
}

# Generate Y|do(X), where X are randomly sampled. 
sampleIntvEst = function(trainModel, sample_X, Z, cutN = 100, N){
  Yhat_intv = rep(0, cutN)
  for (idx in 1:cutN){
    xi = sample_X[idx] # For each sample xi 
    xi_mat = matrix(rep(xi,N),nrow=N) # N length of xi vector
    XZi = cbind(xi_mat, Z) # Matrix such that jth row is (xi, zj) for all j
    colnames(XZi)[1] = 'X' # Change the column name to match
    Yhat_intv[idx] = mean(sampleEst(g_xz, XZi)) # Average(E[Y|X=xi,Z]  ), which is causal mean estimate
  }
  return(Yhat_intv)
}

# Generate Y|x
sampleObsEst = function(trainModel, sample_X, cutN = 100){
  sample_X_input = data.frame(sample_X)
  colnames(sample_X_input)[1] = 'X'
  Yhat_obs = sampleEst(g_x,sample_X_input)
  return(Yhat_obs)
}


###### Spline and plotting functinos ########
resampler <- function(data) {
  n <- nrow(data)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(data[resample.rows,])
}

# Estimating spline 
spline_estimator <- function(data,n) {
  fit = smooth.spline(x=data[,1],y=data[,2],cv=TRUE)
  eval_grid = seq(from=min(data[,1]),to=max(data[,1]),length.out=n)
  return(predict(fit,x=eval_grid)$y) # We only want the predicted values
}

# Confidence interval of spline 
spline_CI <- function(data,n,B,alpha) {
  spline_main = spline_estimator(data,n)
  spline_boots = replicate(B,spline_estimator(resampler(data),n)) # Column are repeated // Spline_estimator are repeated for each 'column' not 'row'
  cis_lower = 2*spline_main - apply(spline_boots,1,quantile,probs=1-alpha/2) # 2X - (X*0.)
  cis_upper = 2*spline_main - apply(spline_boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline_main,lower.ci=cis_lower,upper.ci=cis_upper,
              x=seq(from=min(data[,1]),to=max(data[,1]),length.out=n)))
}

indivPlotResult = function(X,Y,quant_X, N, B, alpha, mycolor){
  set.seed(1)
  df1 = data.frame(X, Y)
  result_spline = spline_CI(df1,N, B=B,alpha=alpha)
  df2 = data.frame( x=result_spline$x, y=result_spline$main.curve   )
  df3 = data.frame( x=result_spline$x, ymin=result_spline$lower.ci, ymax=result_spline$upper.ci   )
  
  gg1 = ggplot(df1, aes(X, Y)) + geom_point(size=1.5, color = mycolor)
  gg1 = gg1 + geom_line(data=df2,aes(x=x, y=y), color=mycolor)
  gg1 = gg1 + geom_ribbon(data=df3, aes(x=x, ymin=ymin,ymax=ymax),alpha=0.4, fill=mycolor)
  gg1 = gg1 + geom_vline(data=df1,xintercept=quant_X, linetype='dotted')
  return(gg1)
}

twoPlotResult = function(X,Y_intv,Y_obs,quant_X,n, B, alpha, color1, color2){
  df1_intv = data.frame(x=X, y=Y_intv)
  df1_obs = data.frame(x=X,y=Y_obs)
  
  result_spline_intv = spline_CI(df1_intv, n, B=B,alpha=alpha)
  result_spline_obs = spline_CI(df1_obs, n, B=B,alpha=alpha)
  
  df2_intv = data.frame( x=result_spline_intv$x, y=result_spline_intv$main.curve)
  df2_obs = data.frame( x=result_spline_obs$x, y=result_spline_obs$main.curve)
  
  df3_intv = data.frame( x=result_spline_intv$x, y=result_spline_intv$main.curve, ymin=result_spline_intv$lower.ci, ymax=result_spline_intv$upper.ci   )
  df3_obs = data.frame( x=result_spline_obs$x, y=result_spline_obs$main.curve, ymin=result_spline_obs$lower.ci, ymax=result_spline_obs$upper.ci   )
  
  gg1 = ggplot()
  gg1 = gg1 + geom_point(data=df1_intv, aes(x=x,y=y), size=1.5,color=color1)
  gg1 = gg1 + geom_point(data=df1_obs, aes(x=x,y=y), size=1.5, color=color2)
  gg1 = gg1 + geom_line(data=df2_intv,aes(x=x, y=y),color=color1)
  gg1 = gg1 + geom_line(data=df2_obs,aes(x=x, y=y), color=color2)
  gg1 = gg1 + geom_ribbon(data=df3_intv, aes(x=x, ymin=ymin, ymax=ymax),alpha=0.4,fill=color1)
  gg1 = gg1 + geom_ribbon(data=df3_obs, aes(x=x, ymin=ymin, ymax=ymax),alpha=0.4, fill=color2)

  gg1 = gg1 + theme_bw()
  gg1 = gg1 + scale_x_continuous(name = "X=x") + scale_y_continuous(name = "Estimated Y",limits=c(min(Y), A=max(Y)))
  gg1 = gg1 + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    text=element_text(size = 16))
  return(gg1)
}

applyTheme = function(gg1, title){
  gg1 = gg1 + ggtitle(title)
  gg1 = gg1 + theme_bw()
  gg1 = gg1 + scale_x_continuous(name = "X") + scale_y_continuous(name = "Estimated Y",limits=c(min(Y), A=max(Y)))
  gg1 = gg1 + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    text=element_text(size = 16))
  return(gg1)
}

drawHistogram = function(X,Y){
  df1 = data.frame(X,Y)
  mytheme = theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                  axis.line.y = element_line(size = 0.5, colour = "black"),
                  axis.title.x= element_blank(),
                  axis.line = element_line(size=1, colour = "black"),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.title=element_text(size = 20))
  hist_X = ggplot(df1, aes(x=X)) + geom_histogram(binwidth = 0.1, fill='#FF9999', colour='black') + scale_y_reverse() + mytheme
  hist_X = hist_X + scale_x_continuous(position='top')
  return(hist_X)
}

mergePlot = function(gg1,gg2,vert_hori){
  grid.newpage()
  if (vert_hori == "vert"){
    gg = plot_grid(gg1, gg2, align = "v", nrow =2, rel_heights = c(3/4, 1/4),axis='tb')
  }
  else if (vert_hori == "hori"){
    gg = plot_grid(gg1, gg2, align = "v", nrow =1, rel_heights = c(1/2, 1/2),axis='tb')
  }
  return(gg)
}



################## MAIN ##################
# Data generation 
source('datagen/data_generation_cano1_cont.R') 

result_sep = sepData(data_obs)
Yobs = result_sep[[1]] 
X = result_sep[[2]]
Z = result_sep[[3]]
XZ = result_sep[[4]] # [X,Z]
N = length(Yobs)

numFold = 5
n = 50 # number of sampling from X

g_xz = trainSL(c(Yobs),XZ,'cont',numFold) # E[Y|x,z]
g_x = trainSL(c(Yobs),data.frame(X),'cont',numFold)

sample_X = reducedX(X,cutN = n)
sample_X_intv = reducedX(X_intv,cutN=n)
quant_X = quantile(X,probs=c(0,0.25,0.5,0.75,1))
Yhat_intv = sampleIntvEst(g_xz, sample_X, Z, n, N)
Yhat_obs = sampleObsEst(g_x, sample_X, cutN = n)

B = 100 # number of bootstraping iteration for generating confidence interval. 
alpha = 0.01 # Confidence interval (100-alpha)%
color_intv = "red"
color_obs = "blue"

data = data.frame(sample_X,Yhat_intv)
spline_main = spline_estimator(data,n)


# Individual 
gg_intv = indivPlotResult(sample_X,Yhat_intv,quant_X,n,B,alpha, color_intv)
gg_intv = applyTheme(gg_intv,"Interventional")
gg_obs = plotResult(sample_X, Yhat_obs,quant_X, n,B,alpha,color_obs)
gg_obs= applyTheme(gg_obs, "Observational")
histX = drawHistogram(X,Y)

gg_merged1 = mergePlot(gg_intv,histX,'vert')
gg_merged2 = mergePlot(gg_obs,histX,'vert')
gg_merged_indiv = mergePlot(gg_merged1,gg_merged2,'hori')

# Merged into one plot 
gg_intv = indivPlotResult(sample_X,Yhat_intv,quant_X,n,B,alpha, color_intv)
gg_obs = plotResult(sample_X, Yhat_obs,quant_X, n,B,alpha,color_obs)
gg_merged = twoPlotResult(sample_X,Yhat_intv,Yhat_obs,quant_X,n,B,alpha,color_intv, color_obs)
gg_merged = applyTheme(gg_merged,'Merged')
gg_merged = mergePlot(gg_merged,histX,'vert')



