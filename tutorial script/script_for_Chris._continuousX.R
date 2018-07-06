# Step 0. Install necessary packages. 
## If you already installed all, then skip the Step 0. 

# install.packages('car')
# install.packages('MASS')
# install.packages('SuperLearner')
# install.packages(c("caret", "glmnet", "randomForest", "RhpcBLASctl", "xgboost", "gam"))

# Step 1. Recall the data. Please refer the comment on data_generation.R 
source('data_generation_continuous_X.R') 
## This generate 'data' (observations)
### 'data' is observational (Y,X,Z)

Yobs = data[,1] # First column of 'data' is Y. 
X = data[,2] # X variable in 'data'.
Z = data[,-c(1,2)] # Z variables in 'data'. 
XZ = data.frame(cbind(X,Z))
N = dim(data)[1] # Number of samples. 

if (length(unique(X)) > 20){
  cat('X is continuous!')
}

# Step 1.5 Prepartaion of X domain  
## This is a list of value of x' from E[Y|do(x')]
## X domain (in this example, [0,1]) are divided into 20 equal intervals. 
### i.e., [0, 0.05, 0.1, 0.15, .., 0.9, 0.95, 1]
quant_X = quantile(X,c(1:10)/10)
# quant_X = seq(min(X),max(X),length.out=20)

## Construct X' matrix (= doX matrix)
### Suppose quant_X = [x1,x2,...,x19,x20].
### Repeat quant_X N times to generate N by 20 matrix. 
### Therefore, doX matrix would be N repeated rows of [x1,x2,...,x19,x20] (N by 20 matrix)
doX = matrix(rep(quant_X,each=N),nrow=N)


# Step 2. Loading SuperLearner and Estimate causal effect E[Y|do(x)]
library(SuperLearner)
# List of supported regression method 
listWrappers()

####### Regression based method #######
# Step 2-1. Construct SuperLearner Model for estimating E[Y|x,z], which denoted as g(x,z). 
## Note that E[Y|x,z] = g(x,z) is a function of x and z. 
## If Y is discrete, that set family = binary(); 
## otherwise Y is continous, set family=gaussian()
g_xz = SuperLearner(Y = c(Yobs), X = XZ, family = gaussian(), cvControl = list(V=5),
                    SL.library = c("SL.xgboost","SL.lm","SL.randomForest"))

g_xz
## Note that g_xz is a model for g(x,z). 
## We can obtain predicted E[Y|x,z], an output of g(x,z), by plugging in (X,Z) onto the model. 

# Step 2-2. Averaging g(x',z) over z, with fixed x' in quant_X
hat_Yx = c() # Initialize empty vector. It will be matrix such that each column is predicted g(xi, z) for xi in quant_X. 

for (col_idx in 1:length(quant_X)){ # For each column in doX
  XiZ = cbind(doX[,col_idx],Z) # Construct (Xi, Z) vector to estimate g(x', z) for each x' in quant_X
  colnames(XiZ)[1] = 'X' # For matching the colume names to XZ (from observations)
  hat_Yxi = predict(g_xz, XiZ, onlySL=T)$pred # vector consisting of g(x',zi) for all i=1,2,...,N
  hat_Yx = cbind(hat_Yx,hat_Yxi)
}

hat_Yx_mean = apply(hat_Yx,2,mean) # Column mean of hat_Yx matrix. This is a vector of estimate of E[Y|do(x'),z] for x' in quant_X
hat_Yx_sd = apply(hat_Yx,2,sd) # Column standard deviation of hat_Yx matrix. This is a vector of estimate of std[Y|do(x'),z] for x' in quant_X
hat_Yx_upper95 = hat_Yx_mean+1.96*hat_Yx_sd/sqrt(N) # Column 95% upper confidence interval
hat_Yx_lower95 = hat_Yx_mean-1.96*hat_Yx_sd/sqrt(N) # Column 95% lower confidence interval

# Step 3. Plot! 
## Loading necessary library 
library(ggplot2) 
### Preparation of graph object 
hat_Yx_graphobj = data.frame(cbind(round(quant_X,4), hat_Yx_mean, hat_Yx_upper95, hat_Yx_lower95))

causal_plot = ggplot(hat_Yx_graphobj, aes(x=quant_X, y=hat_Yx_mean)) +                    # basic graphical object
  geom_smooth(aes(y=hat_Yx_mean), colour="black", span=0.5) + 
  geom_point(size=2)
  # geom_smooth(aes(y=hat_Yx_upper95), colour="gray") + 
  # geom_smooth(aes(y=hat_Yx_lower95), colour="gray")

### Set a X,Y-axis tile 
causal_plot = causal_plot + scale_x_continuous(name = "X=x") + scale_y_continuous(name = "E[Y|do(x)]",limits=c(min(Yobs), A=max(Yobs)))

### Set a graph title 
causal_plot = causal_plot + ggtitle("Causal Effect E[Y|do(x)]")

### Set a background as white 
causal_plot = causal_plot + theme_bw()

### Configurate font-size, etc...
causal_plot = causal_plot +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                  axis.line.y = element_line(size = 0.5, colour = "black"),
                  axis.line = element_line(size=1, colour = "black"), 
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.title=element_text(size = 20),
                  text=element_text(size = 16)
)

### Shading 95% interval 
causal_plot = causal_plot + geom_ribbon(aes(ymin=hat_Yx_lower95, ymax=hat_Yx_upper95), alpha=0.1)

# Step 4. Comparison with E[Y|x]
y_x = SuperLearner(Y = Yobs, X = data.frame(X), family = gaussian(), cvControl = list(V=5),
                    SL.library = c("SL.xgboost","SL.lm","SL.randomForest"))

Yobs_x_mat = c() # Initialize empty vector. It will be matrix such that each column is predicted g(xi, z) for xi in quant_X. 
Xobs_rep = matrix(rep(quant_X,each=N),nrow=N)
for (col_idx in 1:length(quant_X)){ # For each column in doX
  Xi = data.frame(Xobs_rep[,col_idx]) # Construct (Xi, Z) vector to estimate g(x', z) for each x' in quant_X
  colnames(Xi)[1] = 'X' # For matching thTUe colume names to XZ (from observations)
  hat_Yobs_xi = predict(y_x, Xi, onlySL=T)$pred # vector consisting of g(x',zi) for all i=1,2,...,N
  Yobs_x_mat = cbind(Yobs_x_mat,hat_Yobs_xi)
}

Yobs_x_mean = apply(Yobs_x_mat,2,mean) # Column mean of hat_Yx matrix. This is a vector of estimate of E[Y|do(x'),z] for x' in quant_X
Yobs_x_sd = apply(Yobs_x_mat,2,sd) # Column standard deviation of hat_Yx matrix. This is a vector of estimate of std[Y|do(x'),z] for x' in quant_X
Yobs_x_upper95 = Yobs_x_mean+1.96*Yobs_x_sd/sqrt(N) # Column 95% upper confidence interval
Yobs_x_lower95 = Yobs_x_mean-1.96*Yobs_x_sd/sqrt(N) # Column 95% lower confidence interval

# Step 4-1. Comparison with E[Y|x] by plotting 
## Note that NO confidence interval, because confidence interval is 0 
### Because we are not taking over expectation over zi. 

Yobs_x_graphobj = data.frame(cbind(quant_X, Yobs_x_mean, Yobs_x_upper95, Yobs_x_lower95))

obs_plot = ggplot(Yobs_x_graphobj, aes(x=quant_X, y=Yobs_x_mean)) +                    # basic graphical object
  geom_smooth(aes(y=Yobs_x_mean), colour="black", span=0.5) + geom_point(size=2)

### Set a X,Y-axis tile 
obs_plot = obs_plot + scale_x_continuous(name = "X=x") + scale_y_continuous(name = "E[Y|x]",limits=c(min(Yobs), A=max(Yobs)))

### Set a graph title 
obs_plot = obs_plot + ggtitle("Observational Effect E[Y|x]")

### Set a background as white 
obs_plot = obs_plot + theme_bw()

### Configurate font-size, etc...
obs_plot = obs_plot +   theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                                    axis.line.y = element_line(size = 0.5, colour = "black"),
                                    axis.line = element_line(size=1, colour = "black"), 
                                    panel.border = element_blank(),
                                    panel.background = element_blank(),
                                    plot.title=element_text(size = 20),
                                    text=element_text(size = 16)
)

# source('multiplot.R')
# multiplot(bp, bpobs, cols=2)
# install.packages("gridExtra")
library("gridExtra")
grid.arrange(arrangeGrob(causal_plot,obs_plot,nrow=1))
