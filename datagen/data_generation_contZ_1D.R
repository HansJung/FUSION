# Note that we are considering the case G = {(X <- Z -> Y),(X -> Y)}, 
# where Z is a confounding variable of X and Y. 

# Loading necessary document 
library(MASS)
library(car)

# Step 0. Setting
## Set random seed for reproducibility. 
set.seed(1)

## Configuration for data generation 
N = 1000 # Number of data points 
D = 1 # Dimension of the confounding variable Z 

# Step 1. Generate Z 
## In here, Z is D-dimensional multivariate normal random variable. 
fZ = function(D,N){
  mu = 1*rnorm(D) # Set mu random. 
  sigma = replicate(D, abs(rnorm(D))) + 2*diag(D) # Set sigma random. 
  Z = mvrnorm(N, mu, sigma) # Generate N numbers of D-dimensional multivariate normal random variable.
  Z = (Z-mu)/max(sigma)
  return(Z)
}

# Step 2. Generate X 
## In here, X such that Z -> X. 
## Specifically, X = round( exp(Bx * Z) / exp(Bx * Z) + 1 )
fX = function(D,N,Z){
  Bx = -1*(abs(rnorm(D) + rgamma(D,1,1))) # Coefficient of Z. 
  X = +0.5*matrix(Z,ncol=D) %*% matrix(Bx,ncol=1) - 0.1
  X = exp(X)/(exp(X)+1) # Note X is continuous ranged from [0,1]
  X = round(X)
  return(X)
}


# Step 3. Generate Y
## Observational Y such that (X,Z) -> Y
fY = function(D,N, Z,X){
  Byz = -1*(rnorm(D))
  Byx = 1*(rgamma(1,3,3))
  part1 = (matrix(Z,ncol=D) %*% matrix(Byz,ncol=1))
  part2 = 2*(X>0.5)-1
  Y = part1 + 3*part2*X + 2*X
  return(Y)
}

Z = fZ(D,N)
X = fX(D,N,Z)
Y = fY(D,N,Z,X)
data = data.frame(Y,X,Z) # Observational 

X_intv = runif(N)
Y_intv = fY(D,N,Z,X_intv)
data_intv = data.frame(cbind(Y_intv,X_intv,Z))

# plot(data_intv[,2],data_intv[,1])
# plot(data_obs[,2],data_obs[,1])
# hist(X)
