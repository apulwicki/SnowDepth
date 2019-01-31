
################# Import Data ###############################
#install.packages(...)

library(BMS)
library(R.matlab)

importG = readMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/mat2R.mat')
#importG = readMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/BMS/mat2R.mat')

#Glacier data
    params = attr(importG$topoG, "dimnames")[[1]]
    params = c("swe",params)
    elevation = importG$topoG[1]
    #   centreD = importG$topoG[2]
    aspect = importG$topoG[2]
    slope = importG$topoG[3]
    northness = importG$topoG[4]
    curvature = importG$topoG[5]
    Sx = importG$topoG[6]
#G = data.frame(importG$sweG,  elevation, centreD, aspect, slope, northness, curvature, Sx)
G = data.frame(importG$sweG,  elevation, aspect, slope, northness, curvature, Sx)
#G = data.frame(importG$sweG,  elevation, slope, curvature, Sx)
colnames(G) = params
#colnames(G) = c("swe","elevation","slope","curvature","Sx")

rm(aspect, elevation, northness, curvature, slope, Sx, params, importG)
#rm(elevation, curvature, slope, Sx, params, importG)
  

####### Glacier BMS ##########
  
   
    ##uniform prior model
      #mprior is uniform model prior, used UPI
      attG = tryCatch(bms(G, mprior = "uniform", user.int = F), error=function(e) NA)
          #get coefficients
          bad_form = data.frame(matrix(0, ncol = 5, nrow = 5))
         colnames(bad_form) = c("Post_Mean", "Post_SD", "PIP", "Cond_Pos_Sign")
       GC_uni = tryCatch(coef(attG,order.by.pip = F,include.constant = T), error=function(e) bad_form)
 
   Gcoeffs = data.frame(GC_uni)



###### Saving to matlab file  
  writeMat('/home/glaciology1/Documents/Data/SnowDepth/BMS/R2mat.mat',Gcoeffs=Gcoeffs,
           fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
#  writeMat('/Users/Alexandra/Documents/SFU/Data/SnowDepth/BMS/R2mat.mat',Gcoeffs=Gcoeffs,
 #          fixNames=TRUE, matVersion="5", onWrite=NULL, verbose=FALSE)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    