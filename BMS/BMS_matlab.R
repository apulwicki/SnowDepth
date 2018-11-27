
################# Import Data ###############################
#install.packages(...)

library(BMS)
library(R.matlab)

#importG = readMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/mat2R.mat')
importG = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/BMS/mat2R.mat')

#Glacier data
    params = attr(importG$topoG, "dimnames")[[1]]
    params = c("swe",params)
    elevation = importG$topoG[1]
    #centreD = importG$topoG[2]
    #aspect = importG$topoG[3]
    slope = importG$topoG[2]
    #northness = importG$topoG[5]
    curvature = importG$topoG[3]
    Sx = importG$topoG[4]
#G = data.frame(importG$sweG,  elevation, centreD, aspect, slope, northness, curvature, Sx)
G = data.frame(importG$sweG,  elevation, slope, curvature, Sx)
colnames(G) = params

#rm(aspect, elevation, northness, curvature, slope, Sx, centreD, params, importG)
rm(elevation, curvature, slope, Sx, params, importG)
  

####### Glacier BMS ##########
  
    ##uniform prior model
      #mprior is uniform model prior, used UPI
      attG = bms(G, mprior = "uniform", user.int = F)
          #get coefficients
        GC_uni = coef(attG,order.by.pip = F,include.constant = T)
    
    ##binomial prior model
#    att_fixedG = bms(G, mprior = "fixed", mprior.size = 2, user.int = T)
#        GC_fix = coef(att_fixedG, order.by.pip = F, include.constant = T)
    
    ##variable prior model (test how important 'complaints' is)
#    att_pipG = bms(G, mprior = "pip", mprior.size = c(0.01, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), user.int = F)
#        GC_vari = coef(att_pipG,order.by.pip = F, include.constant = T)
    
    ##beta-binomial model priors = default BMS
#    att_randomG = bms(G, mprior = "random", mprior.size = 3, user.int = F) #equivalent att_random = bms(G)
#        GC_rand = coef(att_randomG, order.by.pip = F, include.constant = T)
    
    ##MCMC sampling
#    att_mcmcG = bms(G, burn = 50000, iter = 1e+05, g = "BRIC", mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)
#        GC_mcmc = coef(att_mcmcG,order.by.pip = F,include.constant = T)

    ###return coeffs as structure
#    Gcoeffs = data.frame(GC_uni, GC_fix, GC_rand, GC_vari, GC_mcmc)
#    Gcoeffs = data.frame(GC_uni, GC_fix, GC_rand, GC_vari)
     Gcoeffs = data.frame(GC_uni)

###### Saving to matlab file  
#  writeMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/R2mat.mat',Gcoeffs=Gcoeffs,
#           fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
  writeMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/BMS/R2mat.mat',Gcoeffs=Gcoeffs,
           fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    