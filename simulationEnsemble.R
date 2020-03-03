#compare bit vs trit
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations.
#The idea is to compare how well the trit is doing vs the bit.
#runs in parallel.

rm(list=ls())
library(gplots)
source("simMemoirStrDist3.R")
source("simulation2.R")
source("MLfunctions.R")
library(doParallel)

#request number of CPUs accordingly
os=system("cat ../os.txt",intern = T)
if(os=="mac"){
  registerDoParallel(cores=7)
  pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"
  #AWS
}else if(os=="linux"){
  registerDoParallel(cores=6)
}

generations=c(3,4,5,6)#,7,8,9,10)
# different range for each generation
params = list()
for(g in generations){ params[[g]] = list()}
params[[3]][['barcodes']] = c(10,15,20,30)
params[[4]][['barcodes']] = c(15,20,20,40)
params[[5]][['barcodes']] = c(20,30,40,50)
params[[6]][['barcodes']] = c(30,40,50,60)

params[[3]][['mu']] = c(0.3,0.1)
params[[4]][['mu']] = c(0.2,0.1)
params[[5]][['mu']] = c(0.2,0.1)
params[[6]][['mu']] = c(0.1,0.05)

nRepeats=7

types=c('binary','trit')
#types=c('trit')
simulationType ='trit'

genData = list()
for(g in generations){
    nGen = g
    genData[[g]]= list()
    barcodes = params[[g]][['barcodes']]
    mus = params[[g]][['mu']]
    for(bc in barcodes){
      barcodeLength =bc
      genData[[g]][[bc]] = list()
      for(mu in mus){
          genData[[g]][[bc]][[as.character(mu)]]= compareDist(simulationType=simulationType,alpha_=1/2,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats)
          print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))

      }
    }

    save(genData, file = paste(pathName,"2020_simResults/",as.character(nGen),"_",as.character(bc),"_",as.character(mu),"_",as.character(Sys.time()),sep=""))
}
