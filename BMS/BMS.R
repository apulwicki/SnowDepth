## Work on ##
#image(attG2) and for G13 not working...


################# Import Data ###############################
#install.packages(...)

library(BMS)
library(R.matlab)

T = readMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/mat2R.mat')

#G4 data
    params = attr(T$topoG4, "dimnames")[[1]]
    params = c("swe",params)
    aspect = T$topoG4[1]
    elevation = T$topoG4[2]
    northness = T$topoG4[3]
    profileCurve = T$topoG4[4]
    slope = T$topoG4[5]
    tangentCurve = T$topoG4[6]
    Sx = T$topoG4[7]
    centreD = T$topoG4[8]
  G4 = data.frame(T$sweG4, aspect, elevation, northness, profileCurve, slope, tangentCurve, Sx, centreD)
  colnames(G4) = params
#G2 data
    params = attr(T$topoG2, "dimnames")[[1]]
    params = c("swe",params)
    aspect = T$topoG2[1]
    elevation = T$topoG2[2]
    northness = T$topoG2[3]
    profileCurve = T$topoG2[4]
    slope = T$topoG2[5]
    tangentCurve = T$topoG2[6]
    Sx = T$topoG2[7]
    centreD = T$topoG2[8]
  G2 = data.frame(T$sweG2, aspect, elevation, northness, profileCurve, slope, tangentCurve, Sx, centreD)
  colnames(G2) = params
#G13 data
    params = attr(T$topoG13, "dimnames")[[1]]
    params = c("swe",params)
    aspect = T$topoG13[1]
    elevation = T$topoG13[2]
    northness = T$topoG13[3]
    profileCurve = T$topoG13[4]
    slope = T$topoG13[5]
    tangentCurve = T$topoG13[6]
    Sx = T$topoG13[7]
    centreD = T$topoG13[8]
  G13 = data.frame(T$sweG13, aspect, elevation, northness, profileCurve, slope, tangentCurve, Sx, centreD)
  colnames(G13) = params

rm(aspect, elevation, northness, profileCurve, slope, tangentCurve, Sx, centreD,params,T)


####### G4 BMS ##########

##uniform prior model
  #mprior is uniform model prior, used UPI
  attG4 = bms(G4, mprior = "uniform", user.int = F)
    #get coefficients
    G4C_uni = coef(attG4,order.by.pip = F,include.constant = T)
    #standardized coeffs
    coef(attG4, std.coefs = T, order.by.pip = F, include.constant = T)
    #other results -> mean # of regressors (can also be found using sum(coef(att)[,1]))
    summary(attG4)
    #top models
    topmodels.bma(attG4)[,1:3]
    #get coeffs from these top models
    beta.draws.bma(attG4)[,1:3]
    #overview ->blue = positive, red negative, white noninclusion
    image(attG4)
    image(attG4, yprop2pip = T)
    #examining difference between prior and posterior distibutions
    plotModelsize(attG4) #-> perhas chose different prior?

##binomial prior model
  att_fixedG4 = bms(G4, mprior = "fixed", mprior.size = 2, user.int = T)
    plotModelsize(att_fixedG4)
    G4C_fix = coef(att_fixedG4, order.by.pip = F, include.constant = T)

##variable prior model (test how important 'complaints' is)
  att_pipG4 = bms(G4, mprior = "pip", mprior.size = c(0.01, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), user.int = F)
    G4C_vari = coef(att_pipG4,order.by.pip = F, include.constant = T)
    plotModelsize(att_pipG4)

##beta-binomial model priors = default BMS
  att_randomG4 = bms(G4, mprior = "random", mprior.size = 3, user.int = F) #equivalent att_random = bms(G4)
    plotModelsize(att_randomG4)
    G4C_rand = coef(att_randomG4, order.by.pip = F, include.constant = T)

##MCMC sampling
  att_mcmcG4 = bms(G4, burn = 50000, iter = 1e+05, g = "BRIC", mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)
    plotConv(att_mcmcG4)
    summary(att_mcmcG4)
    G4C_mcmc = coef(att_mcmcG4,order.by.pip = F,include.constant = T)
    plotModelsize(att_mcmcG4)

###comparing models
  #PIP
    plotComp(Uniform = attG4, Fixed= att_fixedG4, PIP = att_pipG4, Random = att_randomG4, MCMC = att_mcmc)
  #Coefs
    plotComp(Uniform = G4C_uni, Fixed= G4C_fix)


###return coeffs as structure
G4coeffs = data.frame(G4C_uni, G4C_fix, G4C_rand, G4C_vari, G4C_mcmc)



####### G2 BMS ##########

##uniform prior model
  #mprior is uniform model prior, used UPI
  attG2 = bms(G2, mprior = "uniform", user.int = F)
    #get coefficients
    G2C_uni = coef(attG2,order.by.pip = F,include.constant = T)
    #standardized coeffs
    coef(attG2, std.coefs = T, order.by.pip = T, include.constant = T)
    #other results -> mean # of regressors (can also be found using sum(coef(att)[,1]))
    summary(attG2)
    #top models
    topmodels.bma(attG2)[,1:3]
    #get coeffs from these top models
    beta.draws.bma(attG2)[,1:3]
    #overview ->blue = positive, red negative, white noninclusion
    image(attG2)
    image(attG2, yprop2pip = T)
    #examining difference between prior and posterior distibutions
    plotModelsize(attG2) #-> perhas chose different prior?

##binomial prior model
  att_fixedG2 = bms(G2, mprior = "fixed", mprior.size = 4, user.int = T)
    plotModelsize(att_fixed)
    G2C_fix = coef(att_fixedG2, order.by.pip = F, include.constant = T)

##variable prior model (test how important 'complaints' is)
  att_pipG2 = bms(G2, mprior = "pip", mprior.size = c(0.01, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), user.int = F)
    G2C_vari = coef(att_pipG2,order.by.pip = F, include.constant = T)
    plotModelsize(att_pipG2)

##beta-binomial model priors = default BMS
  att_randomG2 = bms(G2, mprior = "random", mprior.size = 3, user.int = F) #equivalent att_random = bms(G2)
    plotModelsize(att_randomG2)
    G2C_rand = coef(att_randomG2, order.by.pip = F, include.constant = T)

##MCMC sampling
  att_mcmcG2 = bms(G2, burn = 50000, iter = 1e+05, g = "BRIC", mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)
    plotConv(att_mcmcG2)
    summary(att_mcmcG2)
    G2C_mcmc = coef(att_mcmcG2,order.by.pip = F,include.constant = T)
    plotModelsize(att_mcmcG2)

###comparing models
plotComp(Uniform = attG2, Fixed= att_fixedG2, PIP = att_pipG2, Random = att_randomG2)  

###return coeffs as structure
G2coeffs = data.frame(G2C_uni, G2C_fix, G2C_rand, G2C_vari, G2C_mcmc)  

####### G13 BMS ##########

##uniform prior model
  #mprior is uniform model prior, used UPI
  attG13 = bms(G13, mprior = "uniform", user.int = F)
    #get coefficients
    G13C_uni = coef(attG13,order.by.pip = F,include.constant = T)
    #standardized coeffs
    coef(attG13, std.coefs = T, order.by.pip = T, include.constant = T)
    #other results -> mean # of regressors (can also be found using sum(coef(att)[,1]))
    summary(attG13)
    #top models
    topmodels.bma(attG13)[,1:3]
    #get coeffs from these top models
    beta.draws.bma(attG13)[,1:3]
    #overview ->blue = positive, red negative, white noninclusion
    image(attG13)
    image(attG13, yprop2pip = T)
    #examining difference between prior and posterior distibutions
    plotModelsize(attG13) #-> perhas chose different prior?

##binomial prior model
  att_fixedG13 = bms(G13, mprior = "fixed", mprior.size = 5, user.int = T)
    plotModelsize(att_fixed)
    G13C_fix = coef(att_fixedG13, order.by.pip = F, include.constant = T)

##variable prior model (test how important 'complaints' is)
  att_pipG13 = bms(G13, mprior = "pip", mprior.size = c(0.01, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), user.int = F)
    G13C_vari = coef(att_pipG13,order.by.pip = F, include.constant = T)
    plotModelsize(att_pipG13)

##beta-binomial model priors = default BMS
  att_randomG13 = bms(G13, mprior = "random", mprior.size = 3, user.int = F) #equivalent att_random = bms(G13)
    plotModelsize(att_randomG13)
    G13C_rand = coef(att_randomG13, order.by.pip = F, include.constant = T)

##MCMC sampling
  att_mcmcG13 = bms(G13, burn = 50000, iter = 1e+05, g = "BRIC", mprior = "uniform", nmodel = 2000, mcmc = "bd", user.int = F)
    plotConv(att_mcmcG13)
    summary(att_mcmcG13)
    G13C_mcmc = coef(att_mcmcG13,order.by.pip = F,include.constant = T)
    plotModelsize(att_mcmcG13)

###comparing models
plotComp(Uniform = attG13, Fixed= att_fixedG13, PIP = att_pipG13, Random = att_randomG13)    

###return coeffs as structure
G13coeffs = data.frame(G13C_uni, G13C_fix, G13C_rand, G13C_vari, G13C_mcmc)

  
  
  
########################## Full Topo Prediction #############################################
  
  fullTopo = readMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/fullTopo.mat')
  #G4 full topo data
    params = attr(fullTopo$G4topo, "dimnames")[[1]]
    aspect = fullTopo$G4topo[2]
    elevation = fullTopo$G4topo[3]
    northness = fullTopo$G4topo[4]
    profileCurve = fullTopo$G4topo[5]
    slope = fullTopo$G4topo[6]
    tangentCurve = fullTopo$G4topo[7]
    Sx = fullTopo$G4topo[1]
    centreD = fullTopo$G4topo[8]
  fulltopoG4 = data.frame(aspect, elevation, northness, profileCurve, slope, tangentCurve, Sx, centreD)
  colnames(fulltopoG4) = params
    rm(aspect, elevation, northness, profileCurve, slope, tangentCurve, Sx, centreD,params, fullTopo)

    #prediction for a few values
    pdens = pred.density(attG4, newdata = fulltopoG4[49:52,])
    plot(pdens)
    quantile(pdens, c(0.05, 0.95))
    lps.bma(pdens, fulltopoG4[49:52,1]) #log-predictive score (like rmse?)
    
    
    
###################### cross validation ######################################

#Calibration data (3/4)
    C = G4  
    calInd = C[sample(1:nrow(C),round(nrow(C)*3/4) , replace=FALSE),]

    
    C = matrix(1:20) 
    calInd = matrix(C[sample(1:20,15 , replace=FALSE),])
    
#BMS
    attG4 = bms(calInd, mprior = "uniform", user.int = F)
    #get coefficients
    coefG4 = coef(attG4,order.by.pip = F, include.constant = T)
    coefG4 = coefG4[,2]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    