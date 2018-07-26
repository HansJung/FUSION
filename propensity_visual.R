# Step 0. Install necessary packages. 
library(ggplot2) 
library(SuperLearner)
require(stats); require(graphics)
library(splines)
library(cowplot)

################## Introduction ##################

# This code assumes the following,
# 1. Z is 1 dimensional (1D)
# 2. Z is either continuous or discrete 
# 3. If unique(Z) > 2, then consider it as continuous. 


################## Functions  ##################
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

discCheckZ = function(Z){
  # Check whether Z is continuous or discrete. 
  ## If unique(Z) > 2, then Z is continuous 
  discreteSwitch = TRUE 
  
  # In this function, Z is either 1D or higher. 
  # Suppose Z is ncol(Z) dimensional. 
  # Let Zi be i'th column of Z. 
  for (colidx in 1:ncol(Z)){ # For each column of Z (Zi)
    zi = Z[,colidx] # i'th dimensional 
    if (length(unique(zi)) > 2){
      discreteSwitch = FALSE 
      return(FALSE)
    }
    else{
      next
    }
  }
  return(discreteSwitch)
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

min.mean.sd.max <- function(x) {
  # Extract min, mean-sd, mean, mean+sd, max 
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

drawBoxPlot = function(dfXZ, unique_Z){
  # Currently assume Z is binary 
  dfXZ.X.Z1 = dfXZ[dfXZ[2]==unique_Z[1],1]
  dfXZ.X.Z2 = dfXZ[dfXZ[2]==unique_Z[2],1]
  
  groupZ = rep(1:2,c(length(dfXZ.X.Z1),length(dfXZ.X.Z2)))
  groupZ = factor(groupZ,labels=c('Z1','Z2'))
  
  mydata = data.frame(c(dfXZ.X.Z1,dfXZ.X.Z2), groupZ)
  names(mydata) = c("value","group")
  
  fill <- "#56B4E9"
  line <- "#1F3552"
  
  gg1 = ggplot(aes(y = value, x = factor(group)), data = mydata)
  gg1 = gg1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
                                 outlier.colour = "#1F3552", outlier.shape = 20)
  return(gg1)
}

applyTheme = function(gg1, title,TF_discZ){
  # Apply theme
  gg1 = gg1 + ggtitle(title)
  gg1 = gg1 + coord_cartesian(ylim=c(0,1))
  gg1 = gg1 + theme_bw()
  if (TF_discZ == FALSE){ # continuous Z 
    gg1 = gg1 + scale_x_continuous(name = "Z=z") + scale_y_continuous(name = "P(X=1|Z=z)")  
  }
  else{
    gg1 = gg1 + scale_x_discrete(name = "Z=z") + scale_y_continuous(name = "P(X=1|Z=z)")  
  }
  
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

# Load data 
source('datagen/data_generation_simpson_disc.R')
# source('datagen/data_generation_M.R')

XZ = extractXZ(data)
X = XZ[[1]]
Z = XZ[[2]]
dfXZ = data.frame(X,Z)

dimZ(Z) # Check that dimension of Z is 1 
discCheckZ(Z) # Check whether Z is discrete or Not. Fo 

# Z is in 1D 
## Z is discrete 
unique_Z = unique(Z)[,1]
# prop_score_zval = rep(0,length(unique_Z))
# idx = 1
# for (z_idx in 1:length(unique_Z)){
#   z_val = unique_Z[z_idx]
#   p.x1.zval = mean(dfXZ[dfXZ[2]==z_val,1])
#   prop_score_zval[idx] = p.x1.zval
#   idx = idx + 1 
# }

g_disc = drawBoxPlot(dfXZ, unique_Z)
g_disc = applyTheme(g_disc,"Propensity score",discCheckZ(Z))
g_disc 

# dfXZ.X.Z1 = dfXZ[dfXZ[2]==unique_Z[1],1]
# dfXZ.X.Z2 = dfXZ[dfXZ[2]==unique_Z[2],1]
# 
# groupZ = rep(1:2,c(length(dfXZ.X.Z1),length(dfXZ.X.Z2)))
# groupZ = factor(groupZ,labels=c('Z1','Z2'))
# mydata = data.frame(c(dfXZ.X.Z1,dfXZ.X.Z2), groupZ)
# names(mydata) = c("value","group")
# 
# fill <- "#56B4E9"
# line <- "#1F3552"
# 
# g_disc = ggplot(aes(y = value, x = factor(group)), data = mydata)
# g_disc = g_disc + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
#                                outlier.colour = "#1F3552", outlier.shape = 20)




# Z is continuous
# source('datagen/data_generation_contZ_1D.R')
# 
# XZ = extractXZ(data)
# X = XZ[[1]]
# Z = XZ[[2]]
# dfXZ = data.frame(X,Z)
# 
# dimZ(Z)
# discCheckZ(Z)
# 
# # Computing propensity score 
# px.z = trainSL(X,Z,'disc',5)
# # obtain P(X=1|Z=zi) for each all i=1,2,...,N
# px.z.prediction = predict(px.z,Z,onlySL = T)$pred
# 
# df.plot = data.frame(Z,px.z.prediction)
# colnames(df.plot)[1] = 'z'
# colnames(df.plot)[2] = 'ps'
# 
# g_cont = ggplot(data=df.plot,aes(z,ps)) +
#   geom_point(size=1.5,color='red') +
#   geom_line(data=data.frame(spline(df.plot,n=50)),aes(x,y), color='red')
# 
# g_cont = applyTheme(g_cont,'Propensity score',discCheckZ(Z))
# g_cont
