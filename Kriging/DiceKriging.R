############################################################################
#Dice Kriging
############################################################################

## Libraries
library(R.matlab)
library(DiceKriging)
library(DiceOptim)


## Load my data ##
#residuals = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/residuals.mat')
residuals = readMat('/home/glaciology1/Documents/Data/SnowDepth/Kriging/residuals.mat')
res = residuals$res
utm = data.frame(residuals$utm)
sizexy = residuals$sizexy


## Model ##
m = km(~.,design = utm, response = res, covtype = "matern5_2", nugget.estim = TRUE)
#plot(m)
#m

## Kriging prediction surface ##
x = 0:40:(sizexy[1,2]*40)
y = 1:sizexy[1,1]
grid = expand.grid(X1=x, X2 = y)
pred.m = predict(m,grid,"SK")


pred = matrix(pred.m$mean,n,n)
lower95 = matrix(pred.m$lower95,n,n)
upper95 = matrix(pred.m$upper95,n,n)

writeMat('/home/glaciology1/Documents/Data/SnowDepth/Kriging/kriging.mat',pred=pred, lower95=lower95, upper95=upper95,
         fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
#  writeMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/BMS/R2mat.mat',Gcoeffs=Gcoeffs,
#           fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)

## Install and load rgl package
#library(rgl)

## Plot surface and observations
#plot3d(X[,1],X[,2],res, xlim=c(0,3000),ylim=c(0,3000),zlim=0:1)
#surface3d(x.grid,x.grid, matrix(pred.m$mean,n.grid,n.grid),col="light blue", alpha=0.5)

## Plot surface and observations with intervals
plot3d(utm[,1],utm[,2],res, xlim=c(0,3000),ylim=c(0,3000),zlim=0:1)
surface3d(x,y, matrix(pred.m$mean, sizexy[1,1], sizexy[1,2]),col="light blue", alpha=0.5)
#surface3d(x.grid,x.grid, matrix(pred.m$upper95,n.grid,n.grid),col="light blue", alpha=0.25)
#surface3d(x.grid,x.grid, matrix(pred.m$lower95,n.grid,n.grid),col="light blue", alpha=0.25)
#rgl.snapshot("filename.png")
