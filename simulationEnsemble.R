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
  registerDoParallel(cores=6)
}

# This is the main reference array
# The simulation will ONLY execute the generation specified here
generations=c(3,4,5)#,7,8,9,10)
# different range for each generation
params = list()
for(g in generations){ params[[g]] = list()}
params[[3]][['barcodes']] = c(10,16,24,38,60)
params[[4]][['barcodes']] = c(10,16,24,38,60)
params[[5]][['barcodes']] = c(10,16,24,38,60)
# params[[6]][['barcodes']] = c(30,40,50,60)

params[[3]][['mu']] = c(0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.0001)
params[[4]][['mu']] = c(0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.0001)
params[[5]][['mu']] = c(0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.01,0.0001)
# params[[6]][['mu']] = c(0.1,0.05)

nRepeats=200

types=c('binary','trit')
#types=c('trit')
simulationType ='trit'

for(simulationType in types){
  genData = list()
  #this is the main parameter to change if we want binary simulations
  alpha = if(simulationType =='trit') 1/2 else 0;

  for(g in generations){
      nGen = g
      genData[[g]]= list()
      barcodes = params[[g]][['barcodes']]
      mus = params[[g]][['mu']]
      for(bc in barcodes){
        barcodeLength =bc
        genData[[g]][[bc]] = list()
        for(mu in mus){
            genData[[g]][[bc]][[as.character(mu)]]= compareDist(simulationType=simulationType,alpha_=alpha,
                                                      nGen=nGen,barcodeLength=barcodeLength,
                                                      mu=mu,nRepeats=nRepeats)
            print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))

        }
      }

      save(genData, file = paste(pathName,"2020_simResults/9March200run_",simulationType,"_",as.character(nGen),"_",as.character(bc),"_",as.character(mu),"_",as.character(Sys.time()),sep=""))
  }
}


#to an existent object genData
# we want to add more runs, for parameters not includes originally
# This function should not:
#   write to hard drive the genData object

# test param set (has to ADD to the previous param set)
addParams<-function(){  # e.g., more edit rates:
  params_ = list()
  # generations here has the same value (since we are NOT adding more generations to the data, we could tho)
  for(g in generations){ params_[[g]] = list()}

  params_[[3]][['barcodes']] = c(10,16,24,38,60)
  params_[[4]][['barcodes']] = c(10,16,24,38,60)
  params_[[5]][['barcodes']] = c(10,16,24,38,60)
  # params[[6]][['barcodes']] = c(30,40,50,60)

  params_[[3]][['mu']] = c(0.99,0.01)
  params_[[4]][['mu']] = c(0.99,0.01)
  params_[[5]][['mu']] = c(0.99,0.01)

  return(params)
}

addSimulations<-function(genData, params,simulationType ='trit', nRepeats = 100){

    #this is the main parameter to change if we want binary simulations
    alpha = if(simulationType =='trit') 1/2 else 0;
    for(g in generations){
        nGen = g
        barcodes = params[[g]][['barcodes']]
        mus = params[[g]][['mu']]
        for(bc in barcodes){
          barcodeLength =bc
          for(mu in mus){
              genData[[g]][[bc]][[as.character(mu)]]= compareDist(simulationType=simulationType,alpha_=alpha,
                                                        nGen=nGen,barcodeLength=barcodeLength,
                                                        mu=mu,nRepeats=nRepeats)
              print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))

          }
        }
    }

    return(genData)


}

# This is important. As long as we save the newick files for each simulation, we don't need to save R objects but just load the trees (takes some time)
library(parallel)
readTrees<-function(params,pathName,simulationType, barcodes,generations, mus){
  genData = list()
  for(g in generations){
      genData[[g]] = list()
      nGen = g
      for(bc in barcodes){
        genData[[g]][[bc]] = list()
        barcodeLength = bc

        for(mu in mus){

          #our 100 simulations
          file_pattern = paste(simulationType,"_nG_",as.character(nGen),"_NBC_",as.character(barcodeLength),"mu_",as.character(mu),"_",sep ="")
          #returns a list of length nRepeats (if all goes well)
          these_trees = list.files(path = paste(pathName,"2020trees/",sep=""),pattern= file_pattern)
          these_trees  = paste(pathName,"2020trees/",these_trees,sep="")
          #depending on how/when did we execute the simulations, we might have different N of repeats
          genData[[nGen]][[barcodeLength]][[as.character(mu)]] = matrix(0,length(these_trees),2)


          if(simulationType=='trit'){
            mu_array = rep(mu, barcodeLength)
            alpha_array = rep(alpha, barcodeLength)
          }else{
            mu_array = mu
            alpha_array = alpha
          }

          # READ files to phylo
          tree_list = lapply(these_trees, read.newick)
          # Reconstruct in batch
          # parallel see http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html
          scores = mclapply(tree_list, reconstructSimulation, mu = mu_array,
                            alpha = alpha_array,simulationType = simulationType,
                            roote = T, mc.cores = 7)
          scores = unlist(scores)

          genData[[nGen]][[barcodeLength]][[as.character(mu)]][,2] = scores

          print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))
        }
      }
    }
    return(genData)
}
