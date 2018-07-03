# Reference 
## https://cran.r-project.org/web/packages/SuperLearner/vignettes/SuperLearnerPresent.pdf


# Preparation and installation 
# install.packages("SuperLearner")
install.packages(c("caret", "glmnet", "randomForest", "RhpcBLASctl", "xgboost", "gam"))
 

# Loading data 
install.packages('car')
source('data_generation.R')
set.seed(1) # Set a seed for reproducibility in this random sampling.
colSums(is.na(data)) # Check missingness for each colume 

Yobs = data[,1] # First columne is Y 
Yx0 = dataX0[,1] 
Yx1 = dataX1[,1]

XZ = data[,-1] # X, Z
Xobs = XZ[,1] # X from observational dataset 
Z = XZ[,-1] 
X0Z = dataX0[,-1] # X0
X1Z = dataX1[,-1] # X1 

# Loading SuperLearner packages 
library(SuperLearner)
# List of supported regression method 
listWrappers()

####### Regression based method #######
# Construct SuperLearner Model for estimating E[Y|x,z]
## If Y is discrete, that set family = binary(); 
## otherwise Y is continous, set family=gaussian()
g_xz = SuperLearner(Y = c(Yobs), X = XZ, family = gaussian(), cvControl = list(V=5),
                     SL.library = c("SL.xgboost","SL.lm","SL.randomForest"))

g_xz

# Check how fit the estimate to the truth by compariing the histogram 
g_xz_pred = predict(g_xz, XZ, onlySL = T)
par(mfrow=c(2,1))
plot(density(Yobs))
plot(density(g_xz_pred$pred))

# Estimate E[Y|X=0,z] and $E[Y|X=1,Z]. 
hat_Yx0 = predict(g_xz, X0Z, onlySL=T)$pred # g(0,zi)
hat_Yx1 = predict(g_xz, X1Z, onlySL=T)$pred # g(1,zi)

# Averaging over sample, and compare with the true estimate 
c(mean(Yx0),mean(Yx1)) # True 
c(mean(hat_Yx0),mean(hat_Yx1)) # Regression-based estimate 




####### IPW - Propensity score based estimate #######
# Construct SuperLearner Model for estimating 
## If Y is discrete, that set family = binary(); 
# P(X=x | Z). Regress X on Z 
# P(X=1| Z) trained....
## Note that X is discrete! 
ps_xz = SuperLearner(Y = Xobs, X = Z, family = binomial(), cvControl = list(V=5),
                    SL.library = c("SL.gam"))

ps_xz # trained model 

# Compute propensity scores 
## P(X=1|Z)
c = 1e-8
prop_scores_X1 = predict(ps_xz, Z, onlySL = T)$pred + c 
prop_scores_X0 = (1-predict(ps_xz, Z, onlySL = T)$pred) + c 

# Estimate 

Y_reweighted_X0 = (Yobs/prop_scores_X0) * (1-Xobs)
Y_reweighted_X1 = (Yobs/prop_scores_X1) * Xobs

## Now is the time that we can let users download the data! 
## "CAUSAL" dataset 

# Averaging over sample, and compare with the true estimate 
c(mean(Yx0),mean(Yx1)) # True estimate
c(mean(Y_reweighted_X0),mean(Y_reweighted_X1)) # Regression-based estimate


####### Doubly Robust method: AIPW   #######
Y_AIPW_X0 = (Yobs/prop_scores_X0) * (1-Xobs) - ((1-Xobs) - prop_scores_X0)/(prop_scores_X0)*mean(hat_Yx0)
Y_AIPW_X1 = (Yobs/prop_scores_X1) * (Xobs) - ((Xobs) - prop_scores_X1)/(prop_scores_X1)*mean(hat_Yx1)

c(mean(Yx0),mean(Yx1)) # True 
c(mean(Y_AIPW_X0),mean(Y_AIPW_X1)) # Regression-based estimate


####### Doubly Robust method: TMLE  #######
hat_Y_xz = g_xz_pred$pred

H1 = Xobs/prop_scores_X1
H0 = (1-Xobs)/prop_scores_X0
H = H1 + H0 

# Y - ghat(x,z) = eps * H; 
# eps will be estimated using linear regression. 

# g*(x,z) = ghat(x,z) + eps*H

TMLE_update = lm(Yobs ~ -1 +offset(hat_Y_xz) + H)
eps = TMLE_update$coef

Y_TMLE_X0 =  hat_Yx0+ eps*H0
Y_TMLE_X1 =  hat_Yx1+ eps*H1

c(mean(Yx0),mean(Yx1)) # True 
c(mean(Y_TMLE_X0),mean(Y_TMLE_X1)) # Regression-based estimate

result_table <- data.frame(true = c(mean(Yx0),mean(Yx1)),
                             regressed = c(mean(hat_Yx0),mean(hat_Yx1)),
                             IPW  = c(mean(Y_reweighted_X0),mean(Y_reweighted_X1)), 
                             AIPW = c(mean(Y_AIPW_X0),mean(Y_AIPW_X1)), 
                             TMLE = c(c(mean(Y_TMLE_X0),mean(Y_TMLE_X1)))
                             )

result_table_var <- data.frame(
                           regressed = c(var(hat_Yx0),var(hat_Yx1)),
                           IPW  = c(var(Y_reweighted_X0),var(Y_reweighted_X1)), 
                           AIPW = c(var(Y_AIPW_X0),var(Y_AIPW_X1)), 
                           TMLE = c(c(var(Y_TMLE_X0),var(Y_TMLE_X1)))
)

result_table
result_table_var
