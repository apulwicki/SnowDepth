## Load my data ##
library(R.matlab)

#residuals = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/residuals.mat')
#residuals = readMat('/home/glaciology1/Documents/Data/SnowDepth/Kriging/residuals.mat')
residuals = readMat('D:/Kriging/residuals.mat')
res = residuals$res
utm = residuals$utm


############################################################################
#Dice Kriging
############################################################################

#install.packages('DiceKriging')
#install.packages('DiceOptim')
library(DiceKriging)
library(DiceOptim)

# 1D Example
inputs = c(-1,-0.5,0,0.5,1)
output = c(9,-5,1,-9,11)
theta = 0.4
sigma = 5
trend = c(0,11,2)
model = km(formula = ~., design = data.frame(x=inputs), 
           response = output, covtype = "matern5_2")
           #coef.trend = trend, coef.cov = theta, coef.var = sigma^2)
model

t = seq(from = -2, to = 2, length = 200)
p = predict(model, newdata = data.frame(x=t), type = "SK")

plot(t, p$mean, type = "l", xlim = c(-2, 2), ylim = c(-30, 30),
     xlab = "x", ylab = "y")
lines(t, p$lower95, col = "black", lty = 2)
lines(t, p$upper95, col = "black", lty = 2)
points(inputs, output, col = "red", pch = 19)

## My data 1D ##

inputs = utm[,1]
output = res
model = km(formula = ~., design = data.frame(x=inputs), nugget.estim = TRUE,
           response = output, covtype = "matern5_2")
t = seq(from = 0, to = max(inputs)+100, length = 200)
p = predict(model, newdata = data.frame(x=t), type = "SK")

plot(t, p$mean, type = "l", xlim = c(0, max(t)), ylim = c(0, max(p$mean)+0.1),
     xlab = "x", ylab = "y")
  lines(t, p$lower95, col = "black", lty = 2)
  lines(t, p$upper95, col = "black", lty = 2)
  points(inputs, output, col = "red", pch = 19)



## 2D Example
X <- expand.grid(x1 = seq(0, 1, length = 4), x2 = seq(0, 1, length = 4))
y <- apply(X, 1, branin)
m = km(~., design = X, response = y, covtype = "gauss",nugget.estim = TRUE)

 #concentrates likelihood
n.grid = 30
x.grid = seq(0.01,2,length = n.grid)
X.grid = expand.grid(x.grid,x.grid)
logLik.grid = apply(X.grid,1,logLikFun,m)

contour(x.grid, x.grid, matrix(logLik.grid, n.grid, n.grid), 40,
            xlab = expression(theta[1]), ylab = expression(theta[2]))
opt <- m@covariance@range.val
points(opt[1], opt[2], pch = 19, col = "red")

 #draw kriging mean and visualize pediction accuracy
n.grid = 50
x.grid = seq(0,1,length = n.grid)
X.grid = expand.grid(x1=x.grid, x2 = x.grid)
y.grid = apply(X.grid,1,branin)
pred.m = predict(m,X.grid,"UK")
par(mfrow = c(1,3))

contour(x.grid, x.grid, matrix(y.grid, n.grid, n.grid), 50,
            main = "Branin")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
contour(x.grid, x.grid, matrix(pred.m$mean, n.grid, n.grid), 50,
               main = "Kriging mean")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
contour(x.grid, x.grid, matrix(pred.m$sd^2, n.grid, n.grid), 15,
               main = "Kriging variance")
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

plot(m)


## My data 2D ##
X <- data.frame(utm)
y <- res
m = km(~.,design = X, response = y, covtype = "matern5_2",nugget.estim = TRUE)
# m = km(~.^2, design = X, response = y, covtype = "gauss", nugget = 1e-8 * var(y))
#plot(m)
m

 #concentrated likelihood
n.grid = 30
x.grid = seq(0,300,length = n.grid)
X.grid = expand.grid(x.grid,x.grid)
logLik.grid = apply(X.grid,1,logLikFun,m)

contour(x.grid, x.grid, matrix(logLik.grid, n.grid, n.grid), 40,
        xlab = expression(theta[1]), ylab = expression(theta[2]))
opt <- m@covariance@range.val
points(opt[1], opt[2], pch = 19, col = "red")

#draw kriging mean and visualize pediction accuracy
n.grid = 51
x.grid = seq(0,3000,length = n.grid)
X.grid = expand.grid(X1=x.grid, X2 = x.grid)
pred.m = predict(m,X.grid,"SK")

#contour(x.grid, x.grid, matrix(y.grid, n.grid, n.grid), 50,
 #       main = "Branin")
#points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")
contour(x.grid, x.grid, matrix(pred.m$mean, n.grid, n.grid), 60,
        main = "Kriging mean")
res_range <- range(res)
points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, 
       col = rgb((res-res_range[1])/diff(res_range),0,1-(res-res_range[1])/diff(res_range),1))
# Color above gradients from blue to red as data varies from min to max.

#contour(x.grid, x.grid, matrix(pred.m$sd^2, n.grid, n.grid), 15,
  #      main = "Kriging variance")
#points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

# Install and load rgl package
install.packages("rgl")
library(rgl)

# Plot surface and observations
plot3d(X[,1],X[,2],res, xlim=c(0,3000),ylim=c(0,3000),zlim=0:1)
surface3d(x.grid,x.grid, matrix(pred.m$mean,n.grid,n.grid),col="light blue", alpha=0.5)

# Plot surface and observations with intervals
plot3d(X[,1],X[,2],res, xlim=c(0,3000),ylim=c(0,3000),zlim=0:1)
surface3d(x.grid,x.grid, matrix(pred.m$mean,n.grid,n.grid),col="light blue", alpha=0.5)
surface3d(x.grid,x.grid, matrix(pred.m$upper95,n.grid,n.grid),col="light blue", alpha=0.25)
surface3d(x.grid,x.grid, matrix(pred.m$lower95,n.grid,n.grid),col="light blue", alpha=0.25)
rgl.snapshot("filename.png")

##############################################################################################
#GP fit
##############################################################################################
library(GPfit)

#1D Example
x = seq(from = 0, to = 2*pi, length = 10)
x = (x-min(x))/(max(x)-min(x))
y = sin(x)+runif(10)
y = (y-min(y))/(max(y)-min(y))
plot(x,y)
GPmodel = GP_fit(x,y)
print(GPmodel)
plot(GPmodel)
xnew = seq(from = pi, to = 2.5*pi, length = 10)
GPpredict = predict(GPmodel,xnew)
plot(xnew, GPpredict$Y_hat)
points(x, y, col = "red", pch = 9)

## My data 1D ##
x = utm[,1]
x = (x-min(x))/(max(x)-min(x))
y = res
y = (y-min(y))/(max(y)-min(y))
plot(x,y)
GPmodel = GP_fit(x,y)
print(GPmodel)
plot(GPmodel)
xnew = seq(from = 0, to = max(x), length = 200)
GPpredict = predict(GPmodel,xnew)
plot(xnew, GPpredict$Y_hat)
points(x, y, col = "red", pch = 10)


#2D Example
X = matrix(c(1:10,1:10), ncol = 2)
X[,1] = (X[,1]-min(X[,1]))/(max(X[,1])-min(X[,1]))
X[,2] = (X[,2]-min(X[,2]))/(max(X[,2])-min(X[,2]))
y = sin(1:10)+runif(10)
y = (y-min(y))/(max(y)-min(y))
GPmodel = GP_fit(X,y)
print(GPmodel)
plot(GPmodel)

Xnew = matrix(c(1:8, 1:8), ncol = 2)
Xnew[,1] = (Xnew[,1]-min(Xnew[,1]))/(max(Xnew[,1])-min(Xnew[,1]))
Xnew[,2] = (Xnew[,2]-min(Xnew[,2]))/(max(Xnew[,2])-min(Xnew[,2]))
GPpredict = predict(GPmodel,Xnew)



## My data 2D ##
X = utm
X[,1] = (X[,1]-min(X[,1]))/(max(X[,1])-min(X[,1]))
X[,2] = (X[,2]-min(X[,2]))/(max(X[,2])-min(X[,2]))
y = res
y = (y-min(y))/(max(y)-min(y))
GPmodel = GP_fit(X,y)
print(GPmodel)
plot(GPmodel)

xnew = seq(from = 0, to = max(x), length = 100)
x.grid <- expand.grid(xnew,xnew)
GPpredict = predict(GPmodel,x.grid)
plot(GPmodel, response = FALSE, contour = TRUE)


plot(xnew, GPpredict$Y_hat)
points(x, y, col = "red", pch = 9)


contour(xnew, xnew, matrix(GPpredict$Y_hat, 100, 100), 60,
        main = "Kriging mean")
res_range <- range(res)
points(X, pch = 19, cex = 1.5, 
       col = rgb((res-res_range[1])/diff(res_range),0,1-(res-res_range[1])/diff(res_range),1))
# Color above gradients from blue to red as data varies from min to max.

#contour(x.grid, x.grid, matrix(pred.m$sd^2, n.grid, n.grid), 15,
#      main = "Kriging variance")
#points(X[ , 1], X[ , 2], pch = 19, cex = 1.5, col = "red")

# Plot surface and observations
plot3d(X[,1],X[,2],res, xlim=0:1,ylim=0:1,zlim=0:1)
surface3d(xnew,xnew, matrix(GPpredict$Y_hat,100,100),col="light blue", alpha=0.5)



