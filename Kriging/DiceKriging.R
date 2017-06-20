############################################################################
#Dice Kriging
############################################################################

## Libraries
library(R.matlab)
library(DiceKriging)
library(DiceOptim)
library(foreach)


## Load my data ##
#residuals = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/residuals.mat')
residuals = readMat('/home/glaciology1/Documents/Data/SnowDepth/Kriging/residuals.mat')
res = residuals$res
utm = data.frame(residuals$utm)
sizexy = residuals$sizexy

## Model ##
m = km(~1,design = utm, response = res, covtype = "matern5_2", nugget.estim = TRUE, multistart = 5, iso = TRUE)
 
 #plot(m)
 #m
 
 #Return model paramaters
 maxLL = -m@logLik
 intercept = m@trend.coef
 nugget = m@covariance@nugget
 theta = m@covariance@range.val
    model = data.frame(intercept, nugget, maxLL, theta)

 #Cross validation (leave one out)
 LOO = leaveOneOut.km(m, "SK",trend.reestim = TRUE)
 
 
## Kriging prediction surface ##
x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)
grid = expand.grid(X1=x, X2 = y)
pred.m = predict(m,grid,"SK")


pred = matrix(pred.m$mean, sizexy[1,1], sizexy[1,2], byrow = TRUE)
lower95 = matrix(pred.m$lower95, sizexy[1,1], sizexy[1,2], byrow = TRUE)
upper95 = matrix(pred.m$upper95, sizexy[1,1], sizexy[1,2], byrow = TRUE)


writeMat('/home/glaciology1/Documents/Data/SnowDepth/Kriging/kriging.mat',
         pred=pred, lower95=lower95, upper95=upper95, model = model, LOO = LOO,
         fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
#  writeMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/kriging.mat',
#           pred=pred, lower95=lower95, upper95=upper95, model = model, LOO = LOO,
#           fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)

## Install and load rgl package
#library(rgl)

## Plot surface and observations
#plot3d(utm[,1],utm[,2],res, xlim=c(0,3000),ylim=c(0,3000),zlim=0:1)
#surface3d(x,y, matrix(pred.m$mean,length(x),length(y)),col="light blue", alpha=0.5)

## Plot surface and observations with intervals
#  rglwidget()
#plot3d(utm[,1],utm[,2],res, xlim=c(0,3000),ylim=c(0,3000),zlim=0:1)
#surface3d(x,y, matrix(pred.m$mean, sizexy[1,1], sizexy[1,2]),col="light blue", alpha=0.5)
#surface3d(x.grid,x.grid, matrix(pred.m$upper95,n.grid,n.grid),col="light blue", alpha=0.25)
#surface3d(x.grid,x.grid, matrix(pred.m$lower95,n.grid,n.grid),col="light blue", alpha=0.25)
#rgl.snapshot("filename.png")
