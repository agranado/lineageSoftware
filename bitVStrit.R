#compare bit vs trit
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations.
#The idea is to compare how well the trit is doing vs the bit.
#runs in parallel.

rm(list=ls())
source("simMemoirStrDist2.R")
source("simulation2.R")
source("MLfunctions.R")
library(doParallel)
registerDoParallel(cores=8)

#SET PARAMETERS
barcodes = c(2,4,6,8,10,15,20,25,50,100)
barcodes= c(6,10)
generations=c(3,4,5,6,7,8)
generations=c(8)
mus = c(0.9999,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0001)
#mus = c(0.7)#,0.6)#,0.5,0.4,0.2,0.1)
mus =c(0.999999)
#barcodes = c(6,7)
#generations = c(7,8)
nRepeats=20

types=c('binary','trit')
types=c('trit')

muVariation=list()
for(m in 1:length(mus)){
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
          genData[[ng]]= compareDist(simulationType=simulationType,alpha_=2/3,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats)

       }
       barcodeData[[bc]]=genData
    }
    simType[[simulationType]]=barcodeData
  }
  muVariation[[m]]=simType
}
