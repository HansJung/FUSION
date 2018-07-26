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


################## Functions  ##################
extractXZ = function(df){
  X = df[,2]
  Z = df[,-c(1,2)]
  if (is.null(dim(Z))){
    Z = data.frame(Z)
  }
  return(list(X,Z))
}

dimZ = function(Z){
  return(ncol(Z))
}

discCheckZ = function(Z){
  discreteSwitch = TRUE 
  for (colidx in 1:ncol(Z)){
    zi = Z[,colidx]
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
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

applyTheme = function(gg1, title){
  # Apply theme
  gg1 = gg1 + ggtitle(title)
  gg1 = gg1 + theme_bw()
  gg1 = gg1 + scale_x_continuous(name = "Z=z") + scale_y_continuous(name = "P(X=1|Z=z)",limits=c(min(Y), A=max(Y)))
  gg1 = gg1 + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    text=element_text(size = 16))
  return(gg1)
}

# boxplot_discZ = function(Z0, Z1, px.z0, px.z1, ylim_min, ylim_max){
#   fill <- "#56B4E9"
#   line <- "#1F3552"
#   
#   merged_data = data.frame(cbind(rbind(Z0,Z1),rbind(px.z0, px.z1)))
#   colnames(merged_data) =c('Z','prop_score')
#   merged_data$Z = factor(merged_data$Z,labels=c("Z0","Z1"))
#   
#   merged_data_005 = quantile(merged_data$prop_score, probs = 0.05)
#   merged_data_095 = quantile(merged_data$prop_score, probs = 0.95)
#   
#   bp = ggplot(merged_data, aes(x = Z, y = prop_score)) + 
#     stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
#                  outlier.colour = "#1F3552", outlier.shape = 20)
#   
#   bp = bp + scale_x_discrete(name = "Z=z") + scale_y_continuous(name = "P(X=1|z)")
#   bp = bp + coord_cartesian(ylim=c(ylim_min,ylim_max))
#   #### Customize axis-tick
#   #### Adding title
#   bp = bp + ggtitle("Propensity scores")
#   
#   #### Using white theme
#   bp = bp + theme_bw()
#   bp = bp +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                     axis.line.y = element_line(size = 0.5, colour = "black"),
#                     axis.line = element_line(size=1, colour = "black"),
#                     panel.border = element_blank(),
#                     panel.background = element_blank(),
#                     plot.title=element_text(size = 20),
#                     text=element_text(size = 16)
#   )
#   return(bp) 
# }

################## MAIN ##################

# Load data 
source('datagen/data_generation_simpson_disc.R')
# source('datagen/data_generation_M.R')

XZ = extractXZ(data)
X = XZ[[1]]
Z = XZ[[2]]
dfXZ = data.frame(X,Z)

dimZ(Z)
discCheckZ(Z)

# Z is in 1D 
## Z is discrete 
unique_Z = unique(Z)[,1]
prop_score_zval = rep(0,length(unique_Z))
idx = 1
for (z_idx in 1:length(unique_Z)){
  z_val = unique_Z[z_idx]
  p.x1.zval = mean(dfXZ[dfXZ[2]==z_val,1])
  prop_score_zval[idx] = p.x1.zval
  idx = idx + 1 
}
dfXZ.X.Z1 = dfXZ[dfXZ[2]==unique_Z[1],1]
dfXZ.X.Z2 = dfXZ[dfXZ[2]==unique_Z[2],1]

groupZ = rep(1:2,c(length(dfXZ.X.Z1),length(dfXZ.X.Z2)))
groupZ = factor(groupZ,labels=c('Z1','Z2'))
mydata = data.frame(c(dfXZ.X.Z1,dfXZ.X.Z2), groupZ)
names(mydata) = c("value","group")

fill <- "#56B4E9"
line <- "#1F3552"

g_disc = ggplot(aes(y = value, x = factor(group)), data = mydata)
g_disc = g_disc + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
                               outlier.colour = "#1F3552", outlier.shape = 20)
g_disc = g_disc + coord_cartesian(ylim=c(0,1))
g_disc = g_disc + scale_x_discrete(name = "Z=z") + scale_y_continuous(name = "P(X=1|z)")
g_disc = g_disc + ggtitle("Propensity score when Z disc")
g_disc = g_disc + theme_bw()
g_disc = g_disc + theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                        axis.line.y = element_line(size = 0.5, colour = "black"),
                        axis.line = element_line(size=1, colour = "black"),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        plot.title=element_text(size = 20),
                        text=element_text(size = 16)
)
                        


# Z is continuous 
source('datagen/data_generation_contZ_1D.R')

XZ = extractXZ(data)
X = XZ[[1]]
Z = XZ[[2]]
dfXZ = data.frame(X,Z)

dimZ(Z)
discCheckZ(Z)

px.z = trainSL(X,Z,'disc',5)
px.z.prediction = predict(px.z,Z,onlySL = T)$pred

df.plot = data.frame(Z,px.z.prediction)
colnames(df.plot)[1] = 'z'
colnames(df.plot)[2] = 'ps'

g_cont = ggplot(data=df.plot,aes(z,ps)) + 
  geom_point(size=1.5,color='red') + 
  geom_line(data=data.frame(spline(df.plot,n=50)),aes(x,y), color='red')

g_cont = applyTheme(g_cont,'Propensity score')

