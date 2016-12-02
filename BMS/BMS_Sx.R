library(BMS)
library(R.matlab)
library(caret)

Sx = readMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/SxBMS.mat')

SxG4 = Sx$Sx[1]
SxG4 = data.frame(SxG4)
  findLinearCombos(SxG4)

SxG2 = Sx$Sx[2]
SxG2 = data.frame(SxG2)

SxG13 = Sx$Sx[3]
SxG13 = data.frame(SxG13)


fls1 = bms(SxG4, burn = 50000, iter = 1e+05, g = "BRIC", mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)


data(datafls)
fls1 = bms(datafls, burn = 50000, iter = 1e+05, g = "BRIC", mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)

summary(fls1)
plotConv(fls1[1:100])
