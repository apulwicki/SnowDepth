############################################################################
#Dice Kriging - Universal kriging
############################################################################

## Libraries
# library(R.matlab)
# library(DiceKriging)
# library(DiceOptim)
# library(foreach)

# Testing of UK with various number sof multistart
# Taking 1000 runs as the reference

residuals = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/residuals.mat')
res       = residuals$res
utm       = data.frame(residuals$utm)
sizexy    = residuals$sizexy

m = km(~., design   = utm, 
       response     = res, 
       covtype      = "matern5_2",
       iso          = TRUE, 
       multistart   = 100,
       nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
G4pred100    = predict(m,grid,"UK", se.compute = TRUE)

##

m = km(~., design   = utm, 
       response     = res, 
       covtype      = "matern5_2",
       iso          = TRUE, 
       multistart   = 1000,
       nugget.estim = TRUE)

maxLL     = -m@logLik
intercept = m@trend.coef
nugget    = m@covariance@nugget
theta     = m@covariance@range.val
model     = data.frame(intercept, nugget, maxLL, theta)

x = seq(from = 0, to = (sizexy[1,2]-1)*40, by = 40)
y = seq(from = 0, to = (sizexy[1,1]-1)*40, by = 40)

grid    = expand.grid(X1=x, X2 = y)
G4pred1000    = predict(m,grid,"UK", se.compute = TRUE)

sd( G4pred1000$mean - G4pred100$mean )

## Plotting the data ##
#library(rgl)
# pred = matrix(pred.m$mean, sizexy[1,1], sizexy[1,2], byrow = TRUE)
# trend = matrix(pred.m$trend, sizexy[1,1], sizexy[1,2], byrow = TRUE)
# 
# plot3d(utm[,1],utm[,2],res, xlim=c(0,5000),ylim=c(0,5000),zlim=0:1)
# surface3d(x,y, pred,col="light blue", alpha=0.5)
# surface3d(x,y, trend, xlim=c(0,5000),ylim=c(0,5000), col="dark blue", alpha=0.25)

#summary( pred.m$mean - pred.m2$mean ) 


# writeMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/Kriging/UKmultistarttest.mat',
#          pred1000=pred1000,  pred700=pred700, pred500=pred500, pred400=pred400, pred300=pred300,
#          pred200=pred200, pred150=pred150, pred100=pred100, pred70=pred70,
#          fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
