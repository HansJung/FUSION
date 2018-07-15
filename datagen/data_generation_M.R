# Note that we are considering the case G = {(X <- Z -> Y),(X -> Y)}, 
# where Z is a confounding variable of X and Y. 

# Loading necessary document 
library(MASS)
library(car)

# Step 0. Setting
## Set random seed for reproducibility. 
set.seed(1)

## Configuration for data generation 
N = 2000 # Number of data points 

# Step 1. Generate Z1, Z4
mu1 = 1*rnorm(1) # Set mu random. 
sigma1 = abs(1*rnorm(1)) + runif(1)
Z1 = mvrnorm(N, mu1, sigma1) 

mu4 = 1*rnorm(1) # Set mu random. 
sigma4 = abs(1*rnorm(1)) + runif(1)
Z4 = mvrnorm(N, mu4, sigma4) 

# Step 2. Generate Z2, Z3, Z5
## Z1 -> Z2 
mu2 = 1*rnorm(1) # Set mu random. 
sigma2 = abs(1*rnorm(1)) + runif(1)
Z2 = log(abs(Z1)) + Z1 + mvrnorm(N, mu4, sigma4) 

## Z4 -> Z5 
mu5 = 1*rnorm(1) # Set mu random. 
sigma5 = abs(1*rnorm(1)) + runif(1)
Z5 = exp(Z4) - Z4 - mvrnorm(N, mu5, sigma5) - max(Z4)

## (Z1,Z4) -> Z3 
mu3 = 1*rnorm(1) # Set mu random. 
sigma3 = abs(1*rnorm(1)) + runif(1)
Z3 = 2*Z1 - Z4^2 + abs(Z1 + Z4) + mvrnorm(N, mu3, sigma3)

# Step 3. Generate X
## (Z2, Z3) -> X
muX = 1*rnorm(1) # Set mu random. 
sigmaX = abs(1*rnorm(1)) + runif(1)
X = sqrt(exp(Z1)) - log(abs(Z4)) - abs(Z1 + Z4) + mvrnorm(N, muX, sigmaX)
X = exp(X)/(1+exp(X))
X = round(X)
# X = round(X)
# 
# mean_X = mean(X)
# Q1 = quantile(X)[2]
# Q3 = quantile(X)[4]
# 
# X_box = c()
# for (idx in (1:length(X))){
#   xi = X[idx]
#   if (xi > mean_X){
#     if (abs(mean_X - xi) < abs(Q3 - xi) ){
#       X_box = c(X_box, mean_X)
#     }
#     else{
#       X_box = c(X_box, Q3)
#     }
#   }
#   else{
#     if (abs(mean_X - xi) < abs(Q1 - xi) ){
#       X_box = c(X_box, mean_X)
#     }
#     else{
#       X_box = c(X_box, Q1)
#     }
#   }
# }
# 
# X = X_box

# Step 4. Generate Y
## (Z3, Z5, X) -> Y
muY = 1*rnorm(1) # Set mu random. 
sigmaY = abs(1*rnorm(1)) + runif(1)
Y_noise = mvrnorm(N, muY, sigmaY)
Y = Z3^2 + sqrt(exp(Z5)) + mean(Z3+Z5)*X + Y_noise
Y = exp(Y)/(exp(Y)+1)

# Step 5. Counterfactual Y 
### Note that X is replaced to X0 and X1. 
X0 = matrix(rep(0,N),nrow=N) # N length of 0 vector
X1 = matrix(rep(1,N),nrow=N) # N length of 1 vector

Y0 = Z3^2 + sqrt(exp(Z5)) + mean(Z3+Z5)*X0 + Y_noise
Y0 = exp(Y0)/(exp(Y0)+1)

Y1 = Z3^2 + sqrt(exp(Z5)) + mean(Z3+Z5)*X1 + Y_noise
Y1 = exp(Y1)/(exp(Y1)+1)


# Step 6. Combine all dataset 
### Observational data (Y,X,Z)
data = data.frame(cbind(Y,X,Z1,Z2,Z3,Z4,Z5)) # Observational 

### Intervened data (Y0,X0,Z) and (Y1,X1,Z)
# dataX0 = data.frame(cbind(Y0,X0,Z1,Z2,Z3,Z4,Z5)) # Experimental when X=0
# dataX1 = data.frame(cbind(Y1,X1,Z1,Z2,Z3,Z4,Z5)) # Experimental when X=1
