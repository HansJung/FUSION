# Loading necessary document 
library(MASS)
library(car)
# Step 0. Setting
## Set random seed for reproducibility. 
set.seed(1)

## Configuration for data generation 
N = 2000 # Number of data points 
D = 5 # Dimension of the confounding variable Z 

## SCM
gen_Z1 = function(U4,U5){
  return(U4 - U5)
}
gen_X1 = function(Z1, U2, U6){
  X1 = Norm01(U2 + 5*U6 - 3*Z1)
  X1 = round(X1)
  return(X1)
}
gen_Z2 = function(U5, Z1, X1, U1){
  Z2 = Norm01(exp(U5) + log(abs(Z1)+0.1) - 3*X1 + U1)
  return(Z2)
}
gen_X2 = function(Z2,U2,U3,Z1,X1){
  X2 = Norm01( -exp(Z2)-abs(U2) - U3 - 1.5*Z1 + 2*X1   )
  X2 = round(X2)
  return(X2)
}
gen_Z3 = function(U1, U4, Z2, Z1, X1, X2){
  Norm01( -exp(U1) - exp(U4) + 2*Z2 + Z1 - 2*X1 + log(abs(X2)+0.5) )
}
gen_X3 = function(U3, U6, Z1,Z2,Z3,X1,X2){
  X3 = round(Norm01(exp(U3) - exp(U6) - 2*Z1 - Z2 + 2*Z3 + 2*log(abs(X2+X1)+0.8)  ))
  return(X3)
}

gen_Y = function(Z1,Z2,Z3,X1,X2,X3){
  Norm01(exp(Z1)-exp(Z2)+2*exp(Z3) + log(abs(X1+X2+X3+0.5)) - 1*X1 + 2*X2    )
}

gen_data = function(U1,U2,U3,U4,U5,U6){
  Z1 = gen_Z1(U4,U5)
  X1 = gen_X1(Z1, U2, U6)
  Z2 = gen_Z2(U5, Z1, X1, U1)
  X2 = gen_X2(Z2,U2,U3,Z1,X1)
  Z3 = gen_Z3(U1, U4, Z2, Z1, X1, X2)
  X3 = gen_X3(U3, U6, Z1,Z2,Z3,X1,X2)
  Y = gen_Y(Z1,Z2,Z3,X1,X2,X3)
  
  data = data.frame(cbind(Y,X1,X2,X3,Z1,Z2,Z3)) # Observational 
  return(data)
}

gen_data_cf = function(U1,U2,U3,U4,U5,U6,X1,X2,X3){
  Z1 = gen_Z1(U4,U5)
  Z2 = gen_Z2(U5, Z1, X1, U1)
  Z3 = gen_Z3(U1, U4, Z2, Z1, X1, X2)
  Y = gen_Y(Z1,Z2,Z3,X1,X2,X3)
  
  data_cf = data.frame(cbind(Y,X1,X2,X3,Z1,Z2,Z3)) # Observational 
  return(data_cf)
}


# Step 1. Generate latent variable 

U_gen = function(N){
  mu = 1*rnorm(1) # Set mu random. 
  sigma = abs(1/rgamma(1,1)) + 1 
  U = rnorm(N,mu,sigma)
  U = Norm01(U)
  return(U)
}

Norm01 = function(X){
  return( exp(X)/(1+exp(X)) )
}

U1 = U_gen(N)
U2 = U_gen(N)
U3 = U_gen(N)
U4 = U_gen(N)
U5 = U_gen(N)
U6 = U_gen(N)

# Step 2. Generate Z and X in sequence 
Z1 = gen_Z1(U4,U5)
X1 = gen_X1(Z1, U2, U6)
Z2 = gen_Z2(U5, Z1, X1, U1)
X2 = gen_X2(Z2,U2,U3,Z1,X1)
Z3 = gen_Z3(U1, U4, Z2, Z1, X1, X2)
X3 = gen_X3(U3, U6, Z1,Z2,Z3,X1,X2)
Y = gen_Y(Z1,Z2,Z3,X1,X2,X3)

X1_0 = matrix(rep(0,N),nrow=N) # N length of 0 vector
X2_0 = matrix(rep(0,N),nrow=N) # N length of 0 vector
X3_0 = matrix(rep(0,N),nrow=N) # N length of 0 vector

X1_1 = matrix(rep(1,N),nrow=N) # N length of 1 vector
X2_1 = matrix(rep(1,N),nrow=N) # N length of 1 vector
X3_1 = matrix(rep(1,N),nrow=N) # N length of 1 vector

# Step 3. Combine all dataset 
### Observational data (Y,X1,X2,X3,Z1,Z2,Z3)
data = data.frame(cbind(Y,X1,X2,X3,Z1,Z2,Z3)) # Observational 
data000 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_0,X2_0,X3_0)
data001 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_0,X2_0,X3_1)
data010 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_0,X2_1,X3_0)
data011 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_0,X2_1,X3_1)
data100 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_1,X2_0,X3_0)
data101 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_1,X2_0,X3_1)
data110 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_1,X2_1,X3_0)
data111 = gen_data_cf(U1,U2,U3,U4,U5,U6,X1_1,X2_1,X3_1)


