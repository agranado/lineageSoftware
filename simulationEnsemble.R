#compare bit vs trit
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations.
#The idea is to compare how well the trit is doing vs the bit.
#runs in parallel.

rm(list=ls())
library(gplots)
source("simMemoirStrDist2.R")
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
  pathName="/home/ubuntu/alejandrog/Caltech/lineage/"
  registerDoParallel(cores=72)
}

generations=c(3,4,5,6,7,8,9,10,11,12)
generations = c(10,11,12)
# different range for each generation
params = list()
for(g in generations){ params[[g]] = list()}
# params[[3]][['barcodes']] = c(10,15,20,30)
# params[[4]][['barcodes']] = c(15,20,20,40)
# params[[5]][['barcodes']] = c(20,30,35,40)
# params[[6]][['barcodes']] = c(30,40,45,50)
# params[[7]][['barcodes']] = c(40,50,55,60)
# params[[8]][['barcodes']] = c(40,50,55,60)
# params[[9]][['barcodes']] = c(50,60,65,70)
params[[10]][['barcodes']] = c(60,65,70,75)
params[[11]][['barcodes']] = c(70,80,90,100)
params[[12]][['barcodes']] = c(80,90,100,120)

# params[[3]][['mu']] = c(0.4,0.3,0.1)
# params[[4]][['mu']] = c(0.3,0.2,0.1)
# params[[5]][['mu']] = c(0.3,0.2,0.1)
# params[[6]][['mu']] = c(0.2,0.1,0.05)
# params[[7]][['mu']] = c(0.2,0.1,0.05)
# params[[8]][['mu']] = c(0.2,0.1,0.05,0.01)
# params[[9]][['mu']] = c(0.1,0.05,0.01)
params[[10]][['mu']] = c(0.15,0.1,0.05)
params[[11]][['mu']] = c(0.15,0.1,0.05)
params[[12]][['mu']] = c(0.15,0.1,0.05)



nRepeats=72

types=c('binary','trit')
#types=c('trit')
simulationType ='trit'

run_id = "TueNight_test2"
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
        
          print(as.character(Sys.time()))
        
          genData[[g]][[bc]][[as.character(mu)]]= compareDist(simulationType=simulationType,alpha_=1/2,nGen=nGen,barcodeLength=barcodeLength,mu=mu,nRepeats=nRepeats)
          print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))
          print(as.character(Sys.time()))
      }
    }

    save(genData, file = paste(pathName,"2020_simResults/",as.character(nGen),"_",run_id,"_",as.character(Sys.time()),sep=""))
    print(as.character(Sys.time()))
}
