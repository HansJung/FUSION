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
  # Input
    # DF: data.frame, in the order of (Y,X,Z)
  # Output 
    # list(Y,X,Z,XZ)
  Yobs = DF[,1] # First column of 'data' is Y. 
  X = DF[,2] # X variable in 'data'.
  Z = DF[,-c(1,2)] # Z variables in 'data'. 
  XZ = data.frame(cbind(X,Z))
  return(list(Yobs,X,Z,XZ))
}

# Training using SuperLearner 
trainSL = function(Y,X,family, numFold){
  # Input 
    # Y: output of interest 
    # X: input of interest 
    # family: if Y is continuous, then 'cont'; otherwise 'disc' 
    # numFold: number of cross validation 
  # Output 
    # trained model 
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
  # Input 
    # trainModel: output of function 'trainSL'
    # XZ: data.frame (X,Z)
  # Output 
    # prediced y 
  predResult = predict(trainModel, XZ, onlySL=T)$pred
  return(predResult)
}

# Random sample from X (for reducing the size of X)
reducedX = function(X, cutN = 100){
  # Input 
    # X 
    # cutN: number of samples to be drawn 
  # Output 
    # uniformly drawn samples from X 
  cutN = min(N/10, cutN)
  sample_X = sample(X,cutN)
  return(sample_X)
}

# Generate Y|do(X), where X are randomly sampled. 
sampleIntvEst = function(trainModel, sample_X, Z, cutN = 100, N){
  # Role
    # Using randomly sampled X, estimate E[Y|do(x)] 
  # Input
    # trainModel: output of trainSL (E[Y|x,z]) (g_xz)
    # sample_X: randomly sampled X 
    # Z 
    # cutN: number of samples 
    # N: number of observation 
  # Output 
    # yi=E[Y|do(xi)] for all i=1,2,...,cutN
  
  Yhat_intv = rep(0, cutN)
  for (idx in 1:cutN){
    xi = sample_X[idx] # For each sample xi 
    xi_mat = matrix(rep(xi,N),nrow=N) # N length of xi vector
    XZi = data.frame(cbind(xi_mat, Z)) # Matrix such that jth row is (xi, zj) for all j
    colnames(XZi)[1] = 'X' # Change the column name to match
    Yhat_intv[idx] = mean(sampleEst(g_xz, XZi)) # Average(E[Y|X=xi,Z]  ), which is causal mean estimate
  }
  return(Yhat_intv)
}

# Generate Y|x, where X are randomly sampled. 
sampleObsEst = function(trainModel, sample_X, cutN = 100){
  # Role
    # Using randomly sampled X, estimate E[Y|x] 
  # Input
    # trainModel: output of trainSL (E[Y|x]) (g_x)
    # sample_X: randomly sampled X 
    # Z 
    # cutN: number of samples 
    # N: number of observation 
  # Output 
    # yi=E[Y|xi] for all i=1,2,...,cutN
  sample_X_input = data.frame(sample_X)
  colnames(sample_X_input)[1] = 'X'
  Yhat_obs = sampleEst(trainModel,sample_X_input)
  return(Yhat_obs)
}


###### Spline and plotting functinos ########
# Confidence interval of spline 
spline_CI <- function(data,n,B,alpha) {
  # Role 
    # given a (x,y) pairs, 
    # draw a spline line 
    # and shade an uncertainty area 
  # Input 
    # data: data frame, (x,y) pairs 
    # n: number of points for estimating spline value 
    # B: number of experiment for computing uncertainty. 
      # recommed B >= 100, b/c it doesn't take much computation 
    # alpha: p-value 
  # output
    # main.curve: spline line 
    # lower.ci: lower confidence interval
    # upper.ci: upper confidence interval 
    # x: grid of x (estimating location on X)
  
  spline_main = spline_estimator(data,n) # Spilne estimator from given (x,y) pairs 
  spline_boots = replicate(B,spline_estimator(resampler(data),n)) 
    # spline boots are repeating the following works B times: 
      # resample (x,y) pairs with replacement 
      # given resampled (x,y), run spline 
    # each splined result are stored column by column (not row by row)
    # that is, each column of spline_boots are result of spline_estimator using resampled data. 
  
  cis_lower = 2*spline_main - apply(spline_boots,1,quantile,probs=1-alpha/2) # 2X - (X*0.95) = 
  cis_upper = 2*spline_main - apply(spline_boots,1,quantile,probs=alpha/2) # 2X - (X*)
    # suppose alpha = 0.05, 
    # and we are interested in having upper and lower confidence of vector W
    # where W = spline_main. 
      # Note that W is approximately around at median (0.5). 
      # therefore, 2W is approximated 100% (max).   
      # cis_lower approximates 2W - (0.97.5W) approximates (100% - 95%), 2.5% 
      # cis_upper approxiamtes 2W - (0.025W) approximates (100% - 5%), 97.5%. 
      # resulting approximates 95% confidence interval band 
  
  return(list(main.curve=spline_main,lower.ci=cis_lower,upper.ci=cis_upper,
              x=seq(from=min(data[,1]),to=max(data[,1]),length.out=n)))
}

resampler <- function(data) {
  # Role 
    # given data = (x,y), resample (x,y) with replacement 
  n <- nrow(data)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(data[resample.rows,])
}

# Estimating spline 
spline_estimator <- function(data,n) {
  # given data (x,y), and n points to be estimating the spline, 
  # estimate spline 
  fit = smooth.spline(x=data[,1],y=data[,2],cv=TRUE)
  eval_grid = seq(from=min(data[,1]),to=max(data[,1]),length.out=n)
  return(predict(fit,x=eval_grid)$y) # We only want the predicted values
}

indivPlotResult = function(X,Y,quant_X, n, B, alpha_confidence, mycolor){
  # Role 
    # draw individual plot 
  # Input 
    # X 
    # Y 
    # quant_X: quantile of X 
    # n: number of estimating points of X or Y for spline 
    # B: number of experiment for computing uncertainty.
    # alpha_confidence: confidence interval for spline_CI
    # mycolor: plot color 
  # Output 
    # plot (ggplot output)
  
  set.seed(1)
  df1 = data.frame(X, Y) # data frame of (X,Y)
  result_spline = spline_CI(df1, n, B=B,alpha=alpha_confidence) # resulting spline 
  df2 = data.frame( x=result_spline$x, y=result_spline$main.curve ) 
    # data frame (spline_X, spilne_Y)
  df3 = data.frame( x=result_spline$x, ymin=result_spline$lower.ci, ymax=result_spline$upper.ci   )
    # data frame (spline_X, spline_ymin, spline_ymax)
  
  gg1 = ggplot(df1, aes(X, Y)) + geom_point(size=1.5, color = mycolor)
  gg1 = gg1 + geom_line(data=df2,aes(x=x, y=y), color=mycolor)
  gg1 = gg1 + geom_ribbon(data=df3, aes(x=x, ymin=ymin,ymax=ymax),alpha=0.4, fill=mycolor)
  gg1 = gg1 + geom_vline(data=df1,xintercept=quant_X, linetype='dotted')
  return(gg1)
}

twoPlotResult = function(X,Y_intv,Y_obs,quant_X, n, B, alpha_confidence, color1, color2){
  # Role 
    # draw two line on one plot 
  # Input 
    # X 
    # Y_intv (estimated E[Y|do(x)])
    # Y_obs (estimated E[Y|x])
    # quant_X: quantile of X 
    # n: number of points of X or Y for spline estimation 
    # B: number of experiment for computing uncertainty.
    # alpha_confidence: confidence interval for spline_CI
  # Output 
    # plot (ggplot output)
  
  df1_intv = data.frame(x=X, y=Y_intv) # (X,Y_intv)
  df1_obs = data.frame(x=X,y=Y_obs) # (X,Y_obs)
  
  result_spline_intv = spline_CI(df1_intv, n, B=B,alpha=alpha_confidence)
  result_spline_obs = spline_CI(df1_obs, n, B=B,alpha=alpha_confidence)
  # spline result using df1_intv and df1_obs
  
  df2_intv = data.frame( x=result_spline_intv$x, y=result_spline_intv$main.curve)
  df2_obs = data.frame( x=result_spline_obs$x, y=result_spline_obs$main.curve)
  # spline pairs (x_grid, y_grid), where those pair are estimation location for spline curve 
  
  df3_intv = data.frame( x=result_spline_intv$x, y=result_spline_intv$main.curve, ymin=result_spline_intv$lower.ci, ymax=result_spline_intv$upper.ci   )
  df3_obs = data.frame( x=result_spline_obs$x, y=result_spline_obs$main.curve, ymin=result_spline_obs$lower.ci, ymax=result_spline_obs$upper.ci   )
  # spline pairs (x_grid, y_min, y_max), where (y_max-y_min) shows the confidence interval 
  
  gg1 = ggplot()
  gg1 = gg1 + geom_point(data=df1_intv, aes(x=x,y=y), size=1.5,color=color1)
  gg1 = gg1 + geom_point(data=df1_obs, aes(x=x,y=y), size=1.5, color=color2)
  gg1 = gg1 + geom_line(data=df2_intv,aes(x=x, y=y),color=color1)
  gg1 = gg1 + geom_line(data=df2_obs,aes(x=x, y=y), color=color2)
  gg1 = gg1 + geom_ribbon(data=df3_intv, aes(x=x, ymin=ymin, ymax=ymax),alpha=0.4,fill=color1)
  gg1 = gg1 + geom_ribbon(data=df3_obs, aes(x=x, ymin=ymin, ymax=ymax),alpha=0.4, fill=color2)

  gg1 = gg1 + theme_bw()
  gg1 = gg1 + scale_x_continuous(name = "X=x",limits=c(min(X),A=max(X)),breaks =seq(0,1,length=10) ) + scale_y_continuous(name = "Estimated Y",limits=c(min(Y), A=max(Y)))
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
  # Apply theme
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

drawHistogram = function(X,Y,mybinwidth){
  # Role: draw histogram for X-axis
  # Histogram using (X,Y) pair. 
  
  df1 = data.frame(X,Y)
  mytheme = theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                  axis.line.y = element_line(size = 0.5, colour = "black"),
                  axis.title.x= element_blank(),
                  axis.line = element_line(size=1, colour = "black"),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.title=element_text(size = 20))
  hist_X = ggplot(df1, aes(x=X)) + geom_histogram(binwidth = mybinwidth, fill='#FF9999', colour='black') + scale_y_reverse() + mytheme
  hist_X = hist_X + scale_x_continuous(position='top',limits=c(min(X),A=max(X)))
  return(hist_X)
}

mergePlot = function(gg1,gg2,vert_hori){
  # Merge two plot 
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
# source('datagen/data_generation_cano1_cont.R')
# source('datagen/data_generation_simpson_cont.R') 
source('datagen/data_generation_simpson_cont.R') # Load data_obs 
data_obs = data 
result_sep = sepData(data_obs)
Yobs = result_sep[[1]] 
X = result_sep[[2]]
Z = result_sep[[3]]
XZ = result_sep[[4]] # (X,Z)
N = length(Yobs) # Number of samples of data_obs

numFold = 3
n = 50 # number of sampling from X. This number of points are shown in the resulting plot. 

g_xz = trainSL(c(Yobs),XZ,'cont',numFold) # estimated E[Y|x,z]
g_x = trainSL(c(Yobs),data.frame(X),'cont',numFold) # estimated E[Y|x]

sample_X = reducedX(X,cutN = n) # n points are randomly sampled from X 
sample_X_intv = reducedX(X_intv,cutN=n) 
# n points are randomly sampled from X_intv 
# Note that X_intv are uniformly sampled from the support of X. 

quant_X = quantile(X,probs=c(0,0.25,0.5,0.75,1)) # percetile points of X (0%, 25%, 50%, 75%, 100%)
Yhat_intv = sampleIntvEst(g_xz, sample_X, Z, n, N) # estimated E[Y|do(xi)] for all i = 1,2,...,n
Yhat_obs = sampleObsEst(g_x, sample_X, cutN = n) # estimated E[Y|xi] for all i=1,2,...,n

B = 100 # number of experiment for computing uncertainty (called bootstrap iteration)
alpha = 0.01 # Confidence interval (100-alpha)%
color_intv = "red"
color_obs = "blue"

# Plot drawoing 
## Two individiaul plots side by side 
gg_intv = indivPlotResult(sample_X,Yhat_intv,quant_X,n,B,alpha, color_intv) # individual plot: E[Y|do(x)]
gg_intv = applyTheme(gg_intv,"Interventional") # theme applied for prettier visualization  
gg_obs = indivPlotResult(sample_X, Yhat_obs,quant_X, n,B,alpha,color_obs) # individual plot: E[Y|x]
gg_obs= applyTheme(gg_obs, "Observational") # theme applied for prettier visualization  

mybinwidth = max(1/length(unique(X)),0.2) # binwidth for histogram. If X is binary (0,1), then width = 0.5
histX = drawHistogram(X,Y,mybinwidth) # histogram for X_obs
histX.intv = drawHistogram(X_intv,Y,mybinwidth) # histogram for X_intv

gg_merged1 = mergePlot(gg_intv,histX.intv,'vert') # gg_intv and histX.intv are vertically merged
gg_merged2 = mergePlot(gg_obs,histX,'vert') # gg_obs and histX are vertically merged
gg_merged_indiv = mergePlot(gg_merged1,gg_merged2,'hori') # gg_merged1 and gg_merged2 are horizontally merged (side-by-side)
gg_merged_indiv

## two spline curves (E[Y|do(x)], E[Y|x]) on one plot side by side 
gg_merged = twoPlotResult(sample_X,Yhat_intv,Yhat_obs,quant_X,n,B,alpha,color_intv, color_obs)
gg_merged = applyTheme(gg_merged,'Merged')
gg_merged = mergePlot(gg_merged,histX,'vert')

gg_merged


