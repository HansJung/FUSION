# Step 0. Install necessary packages. 
library(ggplot2) 
library(SuperLearner)
require(stats); require(graphics)
library(splines)

################## Introduction ##################

# This code assumes the following,
## 0. X is binary (X=0,1).
## 1. Z is 1 dimensional (1D) // Multivariate version will be implemented 
## 2. Z is either continuous or discrete 
## 3. If unique(Z) > 2, then consider it as continuous. 


################## Functions  ##################
# Extract (X,Z) from data frame with (Y,X,Z)
extractXZ = function(df){
  # Input: data frame (Y,X,Z)
  # Output: (X,Z) // X is a vector, but Z is a data frame. 
  # From data frame, extract X and Z. 
  X = df[,2]
  Z = df[,-c(1,2)]
  if (is.null(dim(Z))){
    Z = data.frame(Z)
  }
  return(list(X,Z))
}

dimZ = function(Z){
  # dimension of Z 
  return(ncol(Z))
}

addSmallNoise = function(Z, eps){
  # Input Z  in 1D 
  N = dim(Z)[1] # Number of samples 
  Z = Z+eps*rnorm(N)
  return(Z)
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


# Draw a spline plot if Z is continuous. 
drawSplinePlot = function(Z,px.z.prediction, n, B, alpha_confidence, mycolor){
  # Input: Z, P(X=1|Z=zi), for i=1,2,...,N
  # - Z: Z 
  # - px.z.prediction: estimated propensity score 
  # - n: number of estimating points of X or Y for spline 
  # - B: number of experiment for computing uncertainty.
  # - alpha_confidence: confidence interval for spline CI
  # - mycolor: plot of color 
  # Output: plot
  df1 = data.frame(Z,px.z.prediction)
  colnames(df1)[1] = 'x'
  colnames(df1)[2] = 'y'
  # 
  # n = 50 
  # B = 200
  # alpha = 0.01
  # mycolor = 'red'
  result_spline = spline_CI(df1,n,B,alpha_confidence)
  df2 = data.frame( x=result_spline$x, y=result_spline$main.curve ) 
  df3 = data.frame( x=result_spline$x, ymin=result_spline$lower.ci, ymax=result_spline$upper.ci, y=result_spline$main.curve   )
  
  gg1 = ggplot(df1, aes(x=x, y=y)) + geom_point(size=1.5, color = mycolor)
  gg1 = gg1 + geom_line(data=df2,aes(x=x, y=y), color=mycolor)
  gg1 = gg1 + geom_ribbon(data=df3, aes(x=x, ymin=ymin,ymax=ymax),alpha=0.4, fill=mycolor)
  return(gg1)
}


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









# applyTheme = function(gg1, title,TF_discZ){
#   # Apply theme
#   gg1 = gg1 + ggtitle(title)
#   gg1 = gg1 + coord_cartesian(ylim=c(0,1))
#   gg1 = gg1 + theme_bw()
#   if (TF_discZ == FALSE){ # continuous Z 
#     gg1 = gg1 + scale_x_continuous(name = "Z=z") + scale_y_continuous(name = "P(X=1|Z=z)")  
#   }
#   else{
#     gg1 = gg1 + scale_x_discrete(name = "Z=z") + scale_y_continuous(name = "P(X=1|Z=z)")  
#   }
#   
#   gg1 = gg1 + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                     axis.line.y = element_line(size = 0.5, colour = "black"),
#                     axis.line = element_line(size=1, colour = "black"),
#                     panel.border = element_blank(),
#                     panel.background = element_blank(),
#                     plot.title=element_text(size = 20),
#                     text=element_text(size = 16))
#   return(gg1)
# }


# min.mean.sd.max <- function(x) {
#   # Extract min, mean-sd, mean, mean+sd, max 
#   r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
#   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#   return(r)
# }

# discCheckZ = function(Z){
#   # Check whether Z is continuous or discrete. 
#   ## If unique(Z) > 2, then Z is continuous 
#   discreteSwitch = TRUE 
#   
#   # In this function, Z is either 1D or higher. 
#   # Suppose Z is ncol(Z) dimensional. 
#   # Let Zi be i'th column of Z. 
#   for (colidx in 1:ncol(Z)){ # For each column of Z (Zi)
#     zi = Z[,colidx] # i'th dimensional 
#     if (length(unique(zi)) > 2){
#       discreteSwitch = FALSE 
#       return(FALSE)
#     }
#     else{
#       next
#     }
#   }
#   return(discreteSwitch)
# }


# # Draw a box plot (using mean and vairance) if Z is discrete. 
# drawBoxPlot = function(dfXZ, unique_Z){
#   dfXZ.X.Z1 = dfXZ[dfXZ[2]==unique_Z[1],1]
#   dfXZ.X.Z2 = dfXZ[dfXZ[2]==unique_Z[2],1]
#   
#   groupZ = rep(1:2,c(length(dfXZ.X.Z1),length(dfXZ.X.Z2)))
#   groupZ = factor(groupZ,labels=c('Z1','Z2'))
#   
#   mydata = data.frame(c(dfXZ.X.Z1,dfXZ.X.Z2), groupZ)
#   names(mydata) = c("value","group")
#   
#   fill <- "#56B4E9"
#   line <- "#1F3552"
#   
#   gg1 = ggplot(aes(y = value, x = factor(group)), data = mydata)
#   gg1 = gg1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
#                                  outlier.colour = "#1F3552", outlier.shape = 20)
#   return(gg1)
# }




################## MAIN: Case 1, if Z is 1D ##################

# Step 0. Configuration 
eps = 1e-4 # constant for adding small noise 
n = 50 # number of sampling from X. This number of points are shown in the resulting plot. 
B = 100 # number of experiment for computing uncertainty (called bootstrap iteration)
alpha_confidence = 0.01 # Confidence interval (100-alpha)%
plot_color = "red"

# Step 1. Load data 
source('datagen/data_generation_simpson_disc.R') # Z 1D, discrete
# source('datagen/data_generation_contZ_1D.R') # Z 1D, continuous

# Step 2. Extract XZ from data frame (Y,X,Z)
XZ = extractXZ(data)
X = XZ[[1]]
Z = XZ[[2]]
## Add small noises to Z for making Z continuous. 
Z = addSmallNoise(Z,eps)
dimZ(Z) # Check that dimension of Z (in this example, dim(Z) must be 1)
## Wrap up (X,Z)
dfXZ = data.frame(X,Z) 

# Computing propensity score
px.z = trainSL(X,Z,'disc',5)
# obtain P(X=1|Z=zi) for each all i=1,2,...,N
px.z.prediction = predict(px.z,Z,onlySL = T)$pred

# Draw a plot 
g_cont = drawSplinePlot(Z,px.z.prediction, n, B, alpha_confidence, plot_color) # Only a middle spline line without uncertainty shaded area. 
g_cont = applyTheme(g_cont,'Propensity score')
g_cont


# if (discCheckZ(Z) == TRUE){ # If Z is discrete 
#   unique_Z = unique(Z)[,1] # Unique value of Z
#   g_disc = drawBoxPlot(dfXZ, unique_Z) # Draw a box plot 
#   g_disc = applyTheme(g_disc,"Propensity score",discCheckZ(Z))
#   g_disc   
# }else{ # Z is continuous,
#   # Computing propensity score
#   px.z = trainSL(X,Z,'disc',5)
#   # obtain P(X=1|Z=zi) for each all i=1,2,...,N
#   px.z.prediction = predict(px.z,Z,onlySL = T)$pred
#   
#   # Draw a plot 
#   g_cont = drawSpilnePlot(Z,px.z.prediction) # Only a middle spline line without uncertainty shaded area. 
#   g_cont = applyTheme(g_cont,'Propensity score',discCheckZ(Z))
#   g_cont
# }



