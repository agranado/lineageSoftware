#compare bit vs trit 
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations. 
#The idea is to compare how well the trit is doing vs the bit. 
#runs in parallel. 

#this comment was added from the RStudioServer 


rm(list=ls())
library(doParallel)
source("simMemoirStrDist2.R")
source("simulation2.R")
registerDoParallel(cores=32)

#SET PARAMETERS
barcodes = c(10,20,30,50,60,70,100,150)
barcodes= c(100)
generations=c(10,12,15)
#generations=c(8)
mus = c(0.3,0.2,0.1,0.05,0.01)
#mus = c(0.7)#,0.6)#,0.5,0.4,0.2,0.1)
mus =c(0.2,0.1)
#barcodes = c(6,7)
#generations = c(7,8)
nRepeats=32

types=c('binary','trit')
types=c('trit')

#RUN CODE
muVariation=list()
for(m in 1:length(mus)){
  
  start_time =Sys.time()
  print(paste("Start time:", as.character(start_time)) )
  
  
  simType=list()
  mu=mus[m]
  for(st in 1:length(types)){  
    simulationType=types[st]
    barcodeData=list()
    for(bc in 1:length(barcodes)){
       barcodeLength=barcodes[bc]
       genData=list()
       for(ng in 1:length(generations)){
          nGen=generations[ng]
          genData[[ng]]= compareDist(simulationType=simulationType,alpha_=1/2,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats)
          print(paste("sim mu=",toString(mu),"; ",simulationType," nG=",toString(nGen)," BC=",toString(barcodeLength),sep=""))
       }
       barcodeData[[bc]]=genData
    }
    simType[[simulationType]]=barcodeData
  }
  muVariation[[m]]=simType
}
# plotting ----------------------------------------------------------------
# 
# a=array()
# for(ng in 1:length(generations)){
#   a[ng]=simMemoirRandomGuess(generations[ng],mu,alpha,barcodes[1],nRepeats)
# }
# 
# 
# manual.idx=12
# upgam.idx=11
# alpha=2/3
# par(mfrow=c(2,3))
# distBit=array(0,dim=c(length(barcodes),length(generations)))
# distTrit=array(0,dim=c(length(barcodes),length(generations)))
# 
# 
# 
# for(ng in 1:length(generations)){
# 
#   for(bc in 1:length(barcodes)){
#     distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
#     distTrit[bc,ng]=apply(simType[['trit']][[bc]][[ng]],2,mean)[12]
#   }
#   
#   plot(barcodes,distBit[,ng],main=paste("gen = ",toString(generations[ng]),sep=""),type="o",ylim=c(0,a[ng]+0.05*a[ng]),col="blue",ylab="RF dist",xlab="N units")
#   lines(barcodes,distTrit[,ng],type="o")
#   
#   
#   abline(h=a[ng],col="red")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
