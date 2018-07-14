# Step 0. Install necessary packages. 
library(ggplot2) 
library(SuperLearner)
require(stats); require(graphics)
library(splines)
library(cowplot)
library(grid)

# Seperating data frame into Y,X,Z,XZ
noiseAdd = function(X,eps=1e-8,myColName){
  N = length(X)
  X = X + eps*rnorm(N)
  X = data.frame(X)
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
  X.noise = noiseAdd(X,eps,'X')
  Z.noise = noiseAdd(Z,eps,'Z')
  XZ.noise = cbind(X.noise, Z.noise)
  return(list(Yobs,X.noise,Z.noise,XZ.noise))
}

intvData = function(X,eps,do.val){
  if (do.val == 0){
    X = matrix(rep(0,N),nrow=N) # N length of 0 vector
    X.noise = noiseAdd(X,eps,'X')
  }
  else if (do.val==1){
    X = matrix(rep(1,N),nrow=N) # N length of 0 vector
    X.noise = noiseAdd(X,eps,'X')
  }
  return(X)
}

intvSepData = function(X,Z,eps,do.val){
  X = intvData(X,eps,do.val)
  Z = noiseAdd(Z,eps,'Z')
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

reweightY = function(Y, X, propScore, do.val){
  if (do.val == 0){
    reweighted_Y = (Y/propScore) * (1-X) # Y/P(X=1|Z) * I(X=0)  
  }
  else if (do.val==1){
    reweighted_Y = (Y/propScore) * (X) # Y/P(X=1|Z) * I(X=0)  
  }
  return(reweighted_Y)
}

regressMethod = function(g_xz,X0Z,X1Z,do.val){
  if (do.val==0){
    hat_Yx = predict(g_xz, X0Z, onlySL=T)$pred # g(0,zi) for all i
  }
  else if (do.val==1){
    hat_Yx = predict(g_xz, X1Z, onlySL=T)$pred # g(0,zi) for all i
  }
  return(hat_Yx)
}

IPWMethod = function(Y, X, Z, ps_xz, do.val){
  eps = 1e-8
  propScore = propensityScore(ps_xz,Z,eps, do.val)
  reweighted_Y = reweightY(Y,X,propScore, do.val)
  return(sum(reweighted_Y)/colSums(reweighted_Y != 0))
}

AIPWMethod = function(Y,X,prop_scores_X0,prop_scores_X1,hat_Yx0,hat_Yx1,do.val){
  if (do.val == 0){
    Y_AIPW = (Y/prop_scores_X0) * (1-X) - ((1-X) - prop_scores_X0)/(prop_scores_X0)*mean(hat_Yx0)  
  }
  else if (do.val==1){
    Y_AIPW = (Yobs/prop_scores_X1) * (X) - ((X) - prop_scores_X1)/(prop_scores_X1)*mean(hat_Yx1)
  }
  return(Y_AIPW)
}

TMLEMethod = function(X,XZ,X0Z, X1Z,Y,g_xz,prop_scores_X0, prop_scores_X1, do.val){
  hat_Y_xz = predict(g_xz, XZ, onlySL = T)$pred
  H = X/prop_scores_X1 - (1-X)/prop_scores_X0
  TMLE_update = lm(Y ~ -1 +offset(hat_Y_xz) + H)
  TMLE_coef = TMLE_update$coef
  
  hat_Yx = regressMethod(g_xz,X0Z,X1Z,do.val)
  if (do.val==0){
    Hx = -(1-X)/prop_scores_X0 #H(0,z)
    Y_TMLE =  hat_Yx+ TMLE_coef*Hx
  }
  else if (do.val==1){
    Hx = X/prop_scores_X1 #H(1,z)
    Y_TMLE =  hat_Yx+ TMLE_coef*Hx 
  }
  return(Y_TMLE)
}

tableConstruction = function(Yx0,Yx1,Yobs_X0,Yobs_X1,regress_Yx0,
                             regress_Yx1, ipw_Yx0, ipw_Yx1,aipw_Yx0,aipw_Yx1, 
                             tmle_Yx0, tmle_Yx1){
  result_table_mean = data.frame(true = c(mean(Yx0),mean(Yx1)),
                             obs = c(mean(Yobs_X0),mean(Yobs_X1)),
                             regressed = c(mean(regress_Yx0),mean(regress_Yx1)),
                             IPW  = c(mean(ipw_Yx0),mean(ipw_Yx1)),
                             AIPW = c(mean(aipw_Yx0),mean(aipw_Yx1)),
                             TMLE = c(mean(tmle_Yx0),mean(tmle_Yx1)))
  
  result_table_sd = data.frame(true = c(sd(Yx0),sd(Yx1)),
                                 obs = c(sd(Yobs_X0),sd(Yobs_X1)),
                                 regressed = c(sd(regress_Yx0),sd(regress_Yx1)),
                                 IPW  = c(sd(ipw_Yx0),sd(ipw_Yx1)),
                                 AIPW = c(sd(aipw_Yx0),sd(aipw_Yx1)),
                                 TMLE = c(sd(tmle_Yx0),sd(tmle_Yx1)))
  
  return(list(result_table_mean,result_table_sd))
}

plotCausal = function(X0, X1,tmle_Yx0, tmle_Yx1){
  fill <- "#56B4E9"
  line <- "#1F3552"
  
  Y_TMLE_plot = data.frame(cbind(rbind(X0,X1),rbind(tmle_Yx0, tmle_Yx1)))
  colnames(Y_TMLE_plot) =c('X','Yx')
  Y_TMLE_plot$X = factor(Y_TMLE_plot$X,labels=c("X0","X1"))
  
  Y_TMLE_plot_005 = quantile(Y_TMLE_plot$Yx, probs = 0.05)
  Y_TMLE_plot_095 = quantile(Y_TMLE_plot$Yx, probs = 0.95)
  
  bp = ggplot(Y_TMLE_plot, aes(x = X, y = Yx)) + 
    stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
                 outlier.colour = "#1F3552", outlier.shape = 20)
  
  bp = bp + scale_x_discrete(name = "X=x") + scale_y_continuous(name = "E[Y|do(x)]")
  bp = bp + coord_cartesian(ylim=c(0,1))
  #### Customize axis-tick
  #### Adding title
  bp = bp + ggtitle("Causal Effect")
  
  #### Using white theme
  bp = bp + theme_bw()
  bp = bp +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    text=element_text(size = 16)
  )
  return(bp) 
}

plotObs = function(data){
  fill <- "#56B4E9"
  line <- "#1F3552"
  
  Yobs_plot = data.frame(cbind(data$Y,data$X))
  colnames(Yobs_plot) =c('X','Yobs')
  Yobs_plot$X = factor(Yobs_plot$X,labels=c("X0","X1"))
  bpobs = ggplot(Yobs_plot, aes(x = X, y = Yobs)) + 
    stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", fill=fill, colour=line,
                 outlier.colour = "#1F3552", outlier.shape = 20)
  bpobs = bpobs + coord_cartesian(ylim=c(0,1))
  # bpobs = ggplot(Yobs_plot, aes(x = X, y = Yobs)) + geom_boxplot(fill = fill, colour = line, alpha = 0.7,
  #                                                             outlier.colour = "#1F3552", outlier.shape = 20)
  bpobs = bpobs + scale_x_discrete(name = "X=x") + scale_y_continuous(name = "E[Y|x]")
  bpobs = bpobs + ggtitle("Observational Effect E[Y|x]")
  bpobs = bpobs + theme_bw()
  bpobs = bpobs +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                          axis.line.y = element_line(size = 0.5, colour = "black"),
                          axis.line = element_line(size=1, colour = "black"),
                          panel.border = element_blank(),
                          panel.background = element_blank(),
                          plot.title=element_text(size = 20),
                          text=element_text(size = 16)
  )
  return(bpobs)
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

min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}





####################### MAIN ############################
# Data generation 
source('datagen/data_generation_simpson.R') 
eps = 1e-8
resultSep = sepData(data)
Yobs = resultSep[[1]]
X = resultSep[[2]]
Z = resultSep[[3]]
XZ = resultSep[[4]]

resultNoiseSep = noiseSepData(data,eps)
X.noise = resultNoiseSep[[2]]
Z.noise = resultNoiseSep[[3]]
XZ.noise = resultNoiseSep[[4]]

X0Z.noise = intvSepData(data$X,data$Z,eps,do.val=0)
X1Z.noise = intvSepData(data$X,data$Z,eps,do.val=1)

# For validation purpose 
Yx0 = dataX0[,1]
Yx1 = dataX1[,1]

numFold = 3

g_xz = trainSL(c(Yobs),XZ.noise,'disc',numFold) # E[Y|x,z]
ps_xz = trainSL(c(X),X=Z.noise,'disc',numFold) # P(X|z)

# True 
Yx0 = dataX0[,1]
Yx1 = dataX1[,1]

# Obs 
Yobs_X0 = data[data$X==0, 1]
Yobs_X1 = data[data$X==1, 1]

# Regression based method 
regress_Yx0 = regressMethod(g_xz,X0Z.noise,X1Z.noise,do.val=0)
regress_Yx1 = regressMethod(g_xz,X0Z.noise,X1Z.noise,do.val=1)

# IPW 
ipw_Yx0 = IPWMethod(Y,X,Z.noise,ps_xz, do.val=0)
ipw_Yx1 = IPWMethod(Y,X,Z.noise,ps_xz, do.val=1)

propScore_X0 = propensityScore(ps_xz,Z.noise,eps, do.val=0)
propScore_X1 = propensityScore(ps_xz,Z.noise,eps, do.val=1)

# AIPW 
aipw_Yx0 = AIPWMethod(Y,X,propScore_X0,propScore_X1,regress_Yx0,regress_Yx1,do.val=0)
aipw_Yx1 = AIPWMethod(Y,X,propScore_X0,propScore_X1,regress_Yx0,regress_Yx1,do.val=1)

# TMLE X,XZ,X0Z, X1Z,Y,g_xz,prob_scores_X0, prob_scores_X1, do.val
tmle_Yx0= TMLEMethod(X,XZ.noise,X0Z.noise,X1Z.noise,Y,g_xz,propScore_X0,propScore_X1,do.val=0)
tmle_Yx1= TMLEMethod(X,XZ.noise,X0Z.noise,X1Z.noise,Y,g_xz,propScore_X0,propScore_X1,do.val=1)

resultTable = tableConstruction(Yx0,Yx1,Yobs_X0,Yobs_X1,regress_Yx0,
                                regress_Yx1, ipw_Yx0, ipw_Yx1,aipw_Yx0,aipw_Yx1, 
                                tmle_Yx0, tmle_Yx1)
resultTable.mean = resultTable[[1]]
resultTable.sd = resultTable[[2]]

bp = plotCausal(X0, X1,tmle_Yx0, tmle_Yx1)
bpobs = plotObs(data)

gg = mergePlot(bp,bpobs,'hori')
gg
resultTable.mean
resultTable.sd


