# Step 0. Install necessary packages. 
## If you already installed all, then skip the Step 0. 

# install.packages('car')
# install.packages('MASS')
# install.packages('SuperLearner')
# install.packages(c("caret", "glmnet", "randomForest", "RhpcBLASctl", "xgboost", "gam"))

prob_X1 = function(data){
  prob_x1.1 = mean(data$X1)
  prob_x1.0 = 1- prob_x1.1
  return(c(prob_x1.0,prob_x1.1))
}

prob_X2_X1 = function(data,x1val ){
  data_filtered = data[data$X1==x1val,]
  prob_x2.1 = mean(data_filtered$X2)
  prob_x2.0 = 1- prob_x2.1
  return(c(prob_x2.0,prob_x2.1))
}

prob_X3_X1X2 = function(data,x1val,x2val){
  data_filtered = data_filtered = data[(data$X1==x1val) & (data$X2 == x2val),]
  prob_x3.1 = mean(data_filtered$X3)
  prob_x3.0 = 1- prob_x3.1
  return(c(prob_x3.0,prob_x3.1))
}

# Step 1. Recall the data. Please refer the comment on data_generation.R 
source('data_generation_seq.R') 
## This generate 'data' (observations)
### 'data' is observational (Y,X,Z)

Yobs = data[,1] # First column of 'data' is Y. 
X1 = data[,2] # X variable in 'data'.
X2 = data[,3] # X variable in 'data'.
X3 = data[,4] # X variable in 'data'.
Z1 = data[,5] # X variable in 'data'.
Z2 = data[,6] # X variable in 'data'.
Z3 = data[,7] # X variable in 'data'.
XZ = data.frame(cbind(X1,X2,X3,Z1,Z2,Z3))
N = dim(data)[1] # Number of samples. 

# X0 = matrix(rep(0,N),nrow=N) # N length of 0 vector.
# X1 = matrix(rep(1,N),nrow=N) # N length of 1 vector. 
# X000_Z = cbind(X0,Z) # (X0, Z)
# X1Z = cbind(X1,Z) # (X1, Z).

# colnames(X0Z)[1] = 'X' # Matching the column name as XZ.
# colnames(X1Z)[1] = 'X' # Matching the column name as XZ.

# Step 2. Loading SuperLearner and Estimate causal effect E[Y|do(x)]
library(SuperLearner)
# List of supported regression method 
listWrappers()

# Generate model 
## X1 | Z1 
p_x1_z1 = SuperLearner(Y = X1, X = data.frame(Z1), family = binomial(), cvControl = list(V=3),
                       SL.library = c("SL.xgboost","SL.glm","SL.randomForest"))
## X1 
### Nothing to model 
## X2 | X1, Z1, Z2 
p_x2_x1z1z2 = SuperLearner(Y = X1, X = data.frame(cbind(X1,Z1,Z2)), family = binomial(), cvControl = list(V=3),
                       SL.library = c("SL.xgboost","SL.glm","SL.randomForest"))
## X2 | X1 
### Nothing to model 
## X3 | X1, X2, Z1, Z2, Z3
p_x3_x1x2z1z2z3 = SuperLearner(Y = X1, X = data.frame(cbind(X1, X2, Z1, Z2, Z3)), family = binomial(), cvControl = list(V=3),
                           SL.library = c("SL.xgboost","SL.glm","SL.randomForest"))
## X3 | X1, X2 
### Nothing to model 

# Causal query (x1=0,x2=0,x3=0)
## condition 
x1=0
x2=0
x3=0

data_filtered = data[(data$X1 == x1)&(data$X2==x2)&(data$X3==x3),]
Z1_filtered = data.frame(data_filtered$Z1)
colnames(Z1_filtered)[1] = 'Z1'
hat_p_x1_z1 = predict(p_x1_z1, Z1_filtered, onlySL=T)$pred


# 
# 
# ####### Regression based method #######
# # Step 2-1. Construct SuperLearner Model for estimating E[Y|x,z], which denoted as g(x,z). 
# ## Note that E[Y|x,z] = g(x,z) is a function of x and z. 
# ## If Y is discrete, that set family = binary(); 
# ## otherwise Y is continous, set family=gaussian()
# g_xz = SuperLearner(Y = c(Yobs), X = XZ, family = gaussian(), cvControl = list(V=5),
#                     SL.library = c("SL.xgboost","SL.lm","SL.randomForest"))
# 
# g_xz
# ## Note that g_xz is a model for g(x,z). 
# ## We can obtain predicted E[Y|x,z], an output of g(x,z), by plugging in (X,Z) onto the model. 
# 
# # Step 2-2. Averaging g(x',z) over z, with fixed x'. 
# ## Estimate E[Y|X=0,z] and $E[Y|X=1,Z]. 
# hat_Yx0 = predict(g_xz, X0Z, onlySL=T)$pred # g(0,zi) for all i
# hat_Yx1 = predict(g_xz, X1Z, onlySL=T)$pred # g(1,zi) for all i 
# 
# # Step 2-3. Validation: averaging over sample, and compare with the true estimate 
# c(mean(Yx0),mean(Yx1)) # True (hidden)
# c(mean(hat_Yx0),mean(hat_Yx1)) # Regression-based estimate 
# 
# 
# ####### IPW - Propensity score based estimate #######
# # Note that IPW can be used only if X is discrete. 
# 
# # Step 2-1. Construct the SuperLearner Model for training P(x|Z) ('Propensity score')
# ## Training is conducted by regressing Z onto X. 
# ## Note the resulting output of SuperLearner is P(X=1 | Z).
# ## As X is discrete, set 'family = binomial()'. 
# Xobs = X
# ps_xz = SuperLearner(Y = Xobs, X = Z, family = binomial(), cvControl = list(V=5),
#                      SL.library = c("SL.gam", "SL.glm", "SL.xgboost"))
# 
# c = 1e-8 # This is a small quantity added to P(X=1 | Z) for avoding 'dividing to 0' error. 
# ## P(X=0|Z)
# prop_scores_X0 = (1-predict(ps_xz, Z, onlySL = T)$pred) + c # P(X=0|Z)
# ## P(X=1|Z)
# prop_scores_X1 = predict(ps_xz, Z, onlySL = T)$pred + c # P(X=1|Z)
# 
# # Step 2-2. Re-weight Y by inverse probabiltiy (inverse of propensity score). 
# ## For binary X, the indicator that X=0, denoted as I(X=0), is equivalent to 1-X. 
# ## In similar, the indicator that X=1, denoted as I(X=1), is equivalent to X. 
# Y_reweighted_X0 = (Yobs/prop_scores_X0) * (1-Xobs) # Y/P(X=1|Z) * I(X=0)
# Y_reweighted_X1 = (Yobs/prop_scores_X1) * Xobs # Y/P(X=1 | Z) * I(X=1)
# 
# # Step 2-3. Let the user to download the re-weighted data. 
# ## If user's are interested in E[Y|do(X=0)], 
# ## the user might want (Y_reweighted_X0, X, Z)
# downloadable_data_X0 = cbind(Y_reweighted_X0, X, Z)
# ## If user's are interested in E[Y|do(X=1)], 
# ## the user might want (Y_reweighted_X1, X, Z)
# downloadable_data_X1 = cbind(Y_reweighted_X1, X, Z)
# 
# # Step 2-4. Validation: averaging over sample, and compare with the true estimate 
# ## Note that the average is conducted only for nonzero Y. 
# ## Instead of dividing over N, we are dividing over the number of nonzero row. 
# c(mean(Yx0),mean(Yx1)) # True estimate
# c(sum(Y_reweighted_X0)/colSums(Y_reweighted_X0 != 0),sum(Y_reweighted_X1)/colSums(Y_reweighted_X1 != 0)) # IPW estimate
# 
# 
# ####### Doubly Robust method: AIPW   #######
# # Training for Propensity score (P(X=x|Z)) and E[Y|x,z] using regression-based method. 
# # Refer the above. 
# Y_AIPW_X0 = (Yobs/prop_scores_X0) * (1-Xobs) - ((1-Xobs) - prop_scores_X0)/(prop_scores_X0)*mean(hat_Yx0)
# Y_AIPW_X1 = (Yobs/prop_scores_X1) * (Xobs) - ((Xobs) - prop_scores_X1)/(prop_scores_X1)*mean(hat_Yx1)
# 
# c(mean(Yx0),mean(Yx1)) # True 
# c(mean(Y_AIPW_X0),mean(Y_AIPW_X1)) # Regression-based estimate
# 
# 
# 
# ####### Doubly Robust method: TMLE  #######
# # Step 2-1. Estimate E[Y|x,z].
# ## Training model 
# g_xz_pred = predict(g_xz, XZ, onlySL = T)
# ## Estimated E[Y|x,z]
# hat_Y_xz = g_xz_pred$pred 
# 
# # Step 2-2. Estimate H(x,z) = I(X=1)/P(X=1|Z) - I(X=0)/P(X=0|Z)
# H = Xobs/prop_scores_X1 - (1-Xobs)/prop_scores_X0
# H1 = Xobs/prop_scores_X1 #H(1,z)
# H0 = -(1-Xobs)/prop_scores_X0 #H(0,z)
# 
# # Step 2-3. Update the estimated E[Y|x,z] with H. 
# ## Let g1(x,z) be the estimated E[Y|x,z]. 
# ## then g2(x,z) = g1(x,z)+eps*H(x,z). 
# ### Compute the estimated eps by residual of the estimated E[Y|x,z] (i.e., Y-E[Y|x,z]) onto H. 
# TMLE_update = lm(Yobs ~ -1 +offset(hat_Y_xz) + H)
# eps = TMLE_update$coef
# ### TMLE updated (g2(x,z)) is following: 
# hat_y_xz_update = hat_Y_xz + eps * H
# 
# ### Average g2(x,z) over z, after fixing x to do(x')
# Y_TMLE_X0 =  hat_Yx0+ eps*H0
# Y_TMLE_X1 =  hat_Yx1+ eps*H1
# 
# c(mean(Yx0),mean(Yx1)) # True 
# c(mean(Y_TMLE_X0),mean(Y_TMLE_X1)) # Regression-based estimate
# 
# # Step 2.5. Compare the esimating result of Regression based, IPW, AIPW, and TMLE. 
# ## Note that AIPW is most accurate with minimum variance. 
# ## Hence, I would recommend to use AIPW as doubly robust method. 
# result_table <- data.frame(true = c(mean(Yx0),mean(Yx1)),
#                            regressed = c(mean(hat_Yx0),mean(hat_Yx1)),
#                            IPW  = c(mean(Y_reweighted_X0),mean(Y_reweighted_X1)), 
#                            AIPW = c(mean(Y_AIPW_X0),mean(Y_AIPW_X1)), 
#                            TMLE = c(c(mean(Y_TMLE_X0),mean(Y_TMLE_X1)))
# )
# 
# result_table_var <- data.frame(
#   regressed = c(var(hat_Yx0),var(hat_Yx1)),
#   IPW  = c(var(Y_reweighted_X0),var(Y_reweighted_X1)), 
#   AIPW = c(var(Y_AIPW_X0),var(Y_AIPW_X1)), 
#   TMLE = c(c(var(Y_TMLE_X0),var(Y_TMLE_X1)))
# )
# 
# result_table
# result_table_var
# 
# # Step 3. Plot! 
# ## In this example, the AIPW is used. 
# 
# ## Step 3-1. Load the necessary library
# library(ggplot2)
# 
# ## Step 3-2. Encode data with (X,Yx). 
# Y_AIPW_plot = data.frame(cbind(rbind(X0,X1),rbind(Y_AIPW_X0, Y_AIPW_X1)))
# colnames(Y_AIPW_plot) =c('X','Yx')
# Y_AIPW_plot$X = factor(Y_AIPW_plot$X,labels=c("X0","X1"))
# 
# ## Step 3-3. Draw a boxplot 
# ### Refer: http://t-redactyl.io/blog/2016/04/creating-plots-in-r-using-ggplot2-part-10-boxplots.html
# #### Basic Box plot 
# #### colored box plot 
# fill <- "#56B4E9"
# line <- "#1F3552"
# bp = ggplot(Y_AIPW_plot, aes(x = X, y = Yx)) +geom_boxplot()
# bp = ggplot(Y_AIPW_plot, aes(x = X, y = Yx)) + geom_boxplot(fill = fill, colour = line, alpha = 0.7,
#                                                             outlier.colour = "#1F3552", outlier.shape = 20) 
# #### Customize X-axis and Y-axis
# Y_AIPW_plot_005 = quantile(Y_AIPW_plot$Yx, probs = 0.05)
# Y_AIPW_plot_095 = quantile(Y_AIPW_plot$Yx, probs = 0.95)
# 
# 
# bp = bp + scale_x_discrete(name = "X=x") + scale_y_continuous(name = "E[Y|do(x)]", limits=c(min(Yobs), A=max(Yobs)))
# #### Customize axis-tick 
# #### Adding title
# bp = bp + ggtitle("Causal Effect E[Y|do(x)]")
# 
# #### Using white theme
# bp = bp + theme_bw()
# bp = bp +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                   axis.line.y = element_line(size = 0.5, colour = "black"),
#                   axis.line = element_line(size=1, colour = "black"), 
#                   panel.border = element_blank(),
#                   panel.background = element_blank(),
#                   plot.title=element_text(size = 20),
#                   text=element_text(size = 16)
# )
# 
# 
# 
# # Step 4. Comparison with E[Y|X=x]
# ## Y | X=x
# ## X2 is an index representing for X. 
# ## In data matrix, the first column is Y 
# Yobs_X0 = data[data$X2==0, 1]
# Yobs_X1 = data[data$X2==1, 1]
# 
# Yobs_plot = data.frame(cbind(data$X2,data$X1))
# colnames(Yobs_plot) =c('X','Yobs')
# Yobs_plot$X = factor(Yobs_plot$X,labels=c("X0","X1"))
# bpobs = ggplot(Yobs_plot, aes(x = X, y = Yobs)) + geom_boxplot(fill = fill, colour = line, alpha = 0.7,
#                                                                outlier.colour = "#1F3552", outlier.shape = 20) 
# bpobs = bpobs + scale_x_discrete(name = "X=x") + scale_y_continuous(name = "E[Y|x]",limits=c(min(Yobs), A=max(Yobs)))
# bpobs = bpobs + ggtitle("Observational Effect E[Y|x]")
# bpobs = bpobs + theme_bw()
# bpobs = bpobs +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
#                         axis.line.y = element_line(size = 0.5, colour = "black"),
#                         axis.line = element_line(size=1, colour = "black"), 
#                         panel.border = element_blank(),
#                         panel.background = element_blank(),
#                         plot.title=element_text(size = 20),
#                         text=element_text(size = 16)
# )
# # source('multiplot.R')
# # multiplot(bp, bpobs, cols=2)
# # install.packages("gridExtra")
# library("gridExtra")
# grid.arrange(arrangeGrob(bp,bpobs,nrow=1))

