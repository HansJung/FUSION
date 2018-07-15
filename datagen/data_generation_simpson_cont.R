# Loading necessary document 
library(MASS)
library(car)

# Step 0. Setting
## Set random seed for reproducibility. 
set.seed(1)

## Configuration for data generation 
N = 500 # Number of data points 

# Step 1. Generate Z 
## In here, Z is D-dimensional multivariate normal random variable. 
fZ = function(N){
  N.young = N/2
  N.old = N - N.young 
  Z = c(rep(1,N.young),rep(0,N.old))
  Z = sample(Z,N)
  return(Z)
}

fX = function(Z,N){
  X = c()
  for (idx in 1:N){
    zi = Z[idx]
    if (zi == 1){
      xi = rnorm(1,1,0.5)
    }
    else if (zi == 0){
      xi = rnorm(1,2,0.5)
    }
    X = c(X,xi)
  }
  return(X)
}

fY = function(Z,X,N){
  Y = c()
  for (idx in 1:N){
    zi = Z[idx]
    xi = X[idx]
    if (zi==1){
      yi = rnorm(1,5-(zi+1)*(xi),2) 
    }
    else if (zi == 0){
      yi = rnorm(1,10-(zi+1)*(xi),2) 
    }
    Y = c(Y,yi)
  }
  return(Y)
}

Z = fZ(N)
X = fX(Z,N)
Y = fY(Z,X,N)
data = data.frame(cbind(Y,X,Z)) # Observational 

X.min = min(X)
X.max = max(X)

qplot(data[,2],data[,1])
qplot(data[data$Z==0,2],data[data$Z==0,1])
qplot(data[data$Z==1,2],data[data$Z==1,1])

### Intervened X 
### Note that the intervened X (X0 and X1) are NOT affected by Z. 
X_intv = runif(N,X.min,X.max)
Y_intv = fY(Z,X_intv,N)

data_obs = data.frame(cbind(Y,X,Z)) # Observational 
data.intv = data.frame(cbind(Y_intv,X_intv,Z)) # Observational 
qplot(data.intv[,2],data.intv[,1])
