#compare bit vs trit
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations.
#The idea is to compare how well the trit is doing vs the bit.
#runs in parallel.

rm(list=ls())
source("simMemoirStrDist3.R")
source("simulation2.R")
source("MLfunctions.R")
library(doParallel)
os=system("cat ../os.txt",intern = T)
if(os=="mac"){
  registerDoParallel(cores=8)
}else if(os=="linux"){
  registerDoParallel(cores=36)
}

#SET PARAMETERS
barcodes = c(2,4,6,8,10)#,15,20,25,50,100)
barcodes= c(4,6,10,15)
generations=c(3,4,5,6,7,8,9)
generations=c(3,4,5,6,7,8)
mus = c(0.9999,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0001)
#mus = c(0.7)#,0.6)#,0.5,0.4,0.2,0.1)
mus =c(0.9,0.6,0.5,0.4,0.3,0.1)
#barcodes = c(6,7)
#generations = c(7,8)
nRepeats=16

types=c('binary','trit')
#types=c('trit')

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
        
        results= foreach(i=1:nRepeats) %dopar% simMemoirBarcodes(nGen=nGen,alpha=1/2,mu=mu,simulationType=simulationType,barcodeLength = barcodeLength)
        results.mat=do.call(rbind,results)
        genData[[ng]]=results.mat
        
        print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))
      }
      barcodeData[[bc]]=genData
    }
    simType[[simulationType]]=barcodeData
  }
  muVariation[[m]]=simType
}

x11();
#plot average number of unique states for each parameter set: nGen, mu, BC=10
par(mfrow=c(2,3))
for(mIdx in 1:length(mus)){
  #mIdx=5
  bc=3
  ng=1
  
  nRep=dim(muVariation[[4]][['binary']][[3]][[1]])[1]
  
  bit.states=array()
  trit.states=array()
  bit.sdev = array()
  trit.sdev = array()
  for(ng in 1:length(generations)){
    bit.unique=array(); 
    trit.unique=array(); 
    for( i in 1:nRep){
      bit.unique[i]=length(unique(muVariation[[mIdx]][['binary']][[bc]][[ng]][i,]))
      trit.unique[i]=length(unique(muVariation[[mIdx]][['trit']][[bc]][[ng]][i,]))
    }
    bit.states[ng]=mean(bit.unique)
    trit.states[ng]=mean(trit.unique)
    
    trit.sdev[ng]=sd(trit.unique)
    bit.sdev[ng]=sd(bit.unique)
  }
  



  plot(generations,bit.states,ylim=c(1,2^8),type="o",xlab="time (cell divisions)",ylab="Number of unique states",lwd=2,cex.lab=1.5,cex.ax=1.5)
  arrows(generations, bit.states-bit.sdev, generations, bit.states+bit.sdev, length=0.05, angle=90, code=3)
  
  lines(generations,trit.states,type="o",col="blue",lwd=2)
  arrows(generations, trit.states-trit.sdev, generations, trit.states+trit.sdev, length=0.05, angle=90, code=3)
}


