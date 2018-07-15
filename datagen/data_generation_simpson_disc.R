# Loading necessary document 
library(MASS)
library(car)

# Step 0. Setting
## Set random seed for reproducibility. 
set.seed(1)

## Configuration for data generation 
N = 500 # Number of data points 
treatProb.male = 0.25
treatProb.female = 0.75

survProb.treatMale = 0.93
survProb.treatFemale = 0.5
survProb.nonTreatMale = 0.8
survProb.nonTreatFemale = 0.4

# Step 1. Generate Z 
## In here, Z is D-dimensional multivariate normal random variable. 
fZ = function(N){
  N.male = N/2
  N.female = N - N.male 
  Z = c(rep(1,N.male),rep(0,N.female))
  Z = sample(Z,N)
  return(Z)
}

fX = function(Z,N,treatProb.male,treatProb.female){
  X = c()
  for (idx in 1:N){
    zi = Z[idx]
    if (zi == 1){
      xi = rbinom(1,1,treatProb.male)  
    }
    else if (zi == 0){
      xi = rbinom(1,1,treatProb.female)
    }
    X = c(X,xi)
  }
  return(X)
}

fY = function(Z,X,N,survProb.treatMale,survProb.nonTreatMale, survProb.treatFemale, survProb.nonTreatFemale){
  Y = c()
  for (idx in 1:N){
    zi = Z[idx]
    xi = X[idx]
    if ((zi==1)&&(xi==1)){
      yi = rbinom(1,1,survProb.treatMale)  
    }
    else if ((zi==1)&&(xi==0)){
      yi = rbinom(1,1,survProb.nonTreatMale)  
    }
    else if ((zi==0)&&(xi==1)){
      yi = rbinom(1,1,survProb.treatFemale)  
    }
    else if ((zi==0)&&(xi==0)){
      yi = rbinom(1,1,survProb.nonTreatFemale)  
    }
    Y = c(Y,yi)
  }
  return(Y)
}

Z = fZ(N)
X = fX(Z,N,treatProb.male, treatProb.female)
Y = fY(Z,X,N,survProb.treatMale,survProb.nonTreatMale, survProb.treatFemale, survProb.nonTreatFemale)

### Intervened X 
### Note that the intervened X (X0 and X1) are NOT affected by Z. 
X0 = matrix(rep(0,N),nrow=N) # N length of 0 vector
X1 = matrix(rep(1,N),nrow=N) # N length of 1 vector

## Intervened Y (= counterfactual Y) 
### Note that X is replaced to X0 and X1. 
#### (X0, Z) -> Y0; and (X1, Z) -> Y1 
Y0 = fY(Z,X0,N,survProb.treatMale,survProb.nonTreatMale, survProb.treatFemale, survProb.nonTreatFemale)
Y1 = fY(Z,X1,N,survProb.treatMale,survProb.nonTreatMale, survProb.treatFemale, survProb.nonTreatFemale)

data = data.frame(cbind(Y,X,Z)) # Observational 
dataX0 = data.frame(cbind(Y0,X0,Z)) # Experimental when X=0
dataX1 = data.frame(cbind(Y1,X1,Z)) # Experimental when X=1


c(mean(dataX0[,1]),mean(dataX1[,1]))
c(mean(data[data$X==0,1]), mean(data[data$X==1,1]))

