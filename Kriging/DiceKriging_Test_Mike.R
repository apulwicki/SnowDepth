library(R.matlab)
library(DiceKriging)
library(foreach)

# Comparing prediction with ordinary and simple kriging with 
#   single model fit

#
# For clarity, the model here is being fit with ordinary kriging
#   (because the intercept is being estimated from data instead 
#   of just being specified by the user). 
#
# The prediction with simple or ordinary kriging just effects
#   the prediction intervals, not the predicted values as
#   shown below
#

residuals = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/residuals.mat')
res       = residuals$res
utm       = data.frame(residuals$utm)
sizexy    = residuals$sizexy

m = km(~1, design       = utm, 
           response     = res, 
		   covtype      = "matern5_2",
		   iso          = TRUE, 
		   multistart   = 5,
		   nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m  = predict(m,grid,"UK", se.compute = TRUE)
pred.m2 = predict(m,grid,"SK", se.compute = TRUE)

summary(pred.m$mean - pred.m2$mean)

#####################################

# Comparing prediction with ordinary and simple kriging  
#   while fitting twice
#
# In both cases the fitting is ordinary kriging (because
#   the mean is being estimated from data, rather than being 
#   simply stated as being known)
#
# This case should look the same as the above if the estimation 
#   of the model is stable (that you typically get the same 
#   or similar optimal values for the parameters).
#
m = km(~1, design       = utm, 
           response     = res, 
		   covtype      = "matern5_2",
		   iso          = TRUE, 
		   multistart   = 5,
		   nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m  = predict(m,grid,"UK", se.compute = TRUE)

m = km(~1, design       = utm, 
           response     = res, 
		   covtype      = "matern5_2",
		   iso          = TRUE, 
		   multistart   = 5,
		   nugget.estim = TRUE)

maxLL2     = -m@logLik
intercept2 = m@trend.coef
nugget2    = m@covariance@nugget
theta2     = m@covariance@range.val
model2     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m2 = predict(m,grid,"SK", se.compute = TRUE)

summary( pred.m$mean - pred.m2$mean ) 

#
# This seems to indicate that the issue is with the model
#   fitting (the stability of the solution). 
#


#####################################

# Comparing prediction with ordinary and simple kriging  
#   while fitting twice with larger number for multistart

m = km(~1, design       = utm, 
           response     = res, 
		   covtype      = "matern5_2",
		   iso          = TRUE, 
		   multistart   = 50,
		   nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m  = predict(m,grid,"UK", se.compute = TRUE)

m = km(~1, design       = utm, 
           response     = res, 
		   covtype      = "matern5_2",
		   iso          = TRUE, 
		   multistart   = 50,
		   nugget.estim = TRUE)

maxLL2     = -m@logLik
intercept2 = m@trend.coef
nugget2    = m@covariance@nugget
theta2     = m@covariance@range.val
model2     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m2 = predict(m,grid,"SK", se.compute = TRUE)

summary( pred.m$mean - pred.m2$mean ) 

#
# The difference in predictions is now much smaller, indicating
#   that 5 was too few restarts and a larger number is preferred
#

#####################################

# Comparing doing estimation and prediction with universal 
#   kriging with a linear term in both X and Y (named X1 and X2 
#   in the utm variable)
#
# This will likely need even more multistarts because of the 
#   extra difficulty fitting the coefficients for X1 and X2

m = km(~., design       = utm, 
           response     = res, 
		   covtype      = "matern5_2",
		   iso          = TRUE, 
		   multistart   = 200,
		   nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m  = predict(m,grid,"UK", se.compute = TRUE)

## Plotting the data ##
#library(rgl)
pred = matrix(pred.m$mean, sizexy[1,1], sizexy[1,2], byrow = TRUE)
trend = matrix(pred.m$trend, sizexy[1,1], sizexy[1,2], byrow = TRUE)

plot3d(utm[,1],utm[,2],res, xlim=c(0,5000),ylim=c(0,5000),zlim=0:1)
surface3d(x,y, pred,col="light blue", alpha=0.5)
surface3d(x,y, trend, xlim=c(0,5000),ylim=c(0,5000), col="dark blue", alpha=0.25)
#
# This is equivalent (up to the same sort of variability due
#   to the difficult optimization problem) to doing the 
#   universal kriging by:
#

m = km(~ X1 + X2, design       = utm, 
                  response     = res, 
		          covtype      = "matern5_2",
       		      iso          = TRUE, 
		          multistart   = 200,
		          nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
pred.m  = predict(m,grid,"UK", se.compute = TRUE)
