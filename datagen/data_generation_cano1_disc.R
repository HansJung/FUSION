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
D = 3 # Dimension of the confounding variable Z 

# Step 1. Generate Z 
## In here, Z is D-dimensional multivariate normal random variable. 
fZ = function(D){
  mu = 1*rnorm(D) # Set mu random. 
  sigma = replicate(D, abs(rnorm(D))) + 10*diag(D) # Set sigma random. 
  Z = mvrnorm(N, mu, sigma) # Generate N numbers of D-dimensional multivariate normal random variable.
  return(Z)
}
# Step 2. Generate X 
## In here, X such that Z -> X. 

### Observational X 
fX = function(D,Z){
  Bx = -2*(rnorm(D) + rgamma(D,1,1)) - 3 # Coefficient of Z. 
  X = -(matrix(Z,ncol=D) %*% matrix(Bx,ncol=1) + 2*rnorm(N))
  X = round(exp(X)/(exp(X)+1)) # Note X is binary.  
  return(X)
}

### Intervened X 
### Note that the intervened X (X0 and X1) are NOT affected by Z. 
X0 = matrix(rep(0,N),nrow=N) # N length of 0 vector
X1 = matrix(rep(1,N),nrow=N) # N length of 1 vector

# Step 3. Generate Y
## Observational Y such that (X,Z) -> Y
fY = function(D,Z,X){
  Byz = -3*(rnorm(D))
  Byx = 2*(rgamma(1,3,3))
  # Y = 2*(matrix(Z,ncol=D) %*% matrix(Byz,ncol=1)) - 8*Byx*X + 1*rnorm(N)
  # Y = 2*abs((matrix(Z,ncol=D) %*% matrix(Byz,ncol=1)) +
  #       - 5*Byx*X) + 1*rnorm(N)
  
  Y = abs((matrix(Z,ncol=D) %*% matrix(Byz,ncol=1)) -
            5*(abs((matrix(Z,ncol=D) %*% matrix(Byz,ncol=1)))) +
            - 2*Byx*X) + 1*abs(2*X) + 1*rnorm(N)
  # Y = 10*(exp(Y)/(1+exp(Y)))
  return(Y)
}


# Step 4. Combine all dataset 
## Combine X,Y,Z 
### Observational data (Y,X,Z)

Z = fZ(D)
X = fX(D,Z)
Y = fY(D,Z,X)
data = data.frame(cbind(Y,X,Z)) # Observational 

## Intervened Y (= counterfactual Y) 
### Note that X is replaced to X0 and X1. 
#### (X0, Z) -> Y0; and (X1, Z) -> Y1 
Y0 = fY(D,Z,X0)
Y1 = fY(D,Z,X1)
### Intervened data (Y0,X0,Z) and (Y1,X1,Z)
dataX0 = data.frame(cbind(Y0,X0,Z)) # Experimental when X=0
dataX1 = data.frame(cbind(Y1,X1,Z)) # Experimental when X=1
