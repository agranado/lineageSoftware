
# Aug 14th 2019
# Integrase paper writing
#
# Simulation is restricted to toplogies that are observed experimentally
# Goal: Establish a baseline for reconstruction accuracy

library(data.table)
library(dplyr)
library(ape)
library(data.tree)
library(stringr)
library(cluster)
library(phangorn)

#source('../GIT/MLfunctions2.R')

# Set up the full path for the data
# Barcode data is in the folder FISH/
# Main folder is ../integrase_folder/10mer_2019/
if(length(grep("linux",read.table("../os.txt")$V1))){
    source("../integrase-data/RF.experiment.R")
    integrase_folder= "integrase-data/" # For Ubuntu
    file.path=paste("/home/agranado/MEGA/Caltech/trees/",integrase_folder,"10mer_2019/",sep="") # Ubuntu

}else{
    source("../GraceData/RF.experiment.R")
    integrase_folder = "GraceData/"     # For Mac
    file.path=paste("/users/alejandrog/MEGA/Caltech/trees/",integrase_folder,"10mer_2019/",sep="") # Mac
}




# Functions for simulation of recording restricted to a given topology.
# we don't need to simulate the lineage but just follow the one that was observed in the experiment.
# use the 10mer (10 recording units) and simulate the edits with experimentally-measured parameteres
# parameters are trit-dependent and therefore each unit has its own mu, alpha

# # old functions (simplified for this particular situation)
# mutationScratchpad_2019 <- function(barcode,mu,alpha,type='trit',recType="integrase",cInts=1,tInts=1,chr.acc=c()){
#
#   originalBarcode=strsplit(barcode,"")[[1]]
#   #now a is an array of char elements representing the scratchpad
#   m_event<-runif(1,0,1)
#   #mutation happens (rate = mutations/generation)
#   #for each element in the array, ask if mutation happens
#   #with constant independent rate (mutation in [1] does not affect [2])
#   if(recType=="integrase"){
#     mutatedBarcode=recIntegrase(originalBarcode,mu=mu,alpha=alpha,type=type,currentInts=cInts,totalInts=tInts)
#   }else if(recType=="epimemoir"){ #epimemoir
#     mutatedBarcode=recEpiMemoir(originalBarcode,mu=mu,alpha=alpha,type=type,totalInts=tInts,chr.acc=chr.acc)
#   }else { #no mutation
#     mutatedBarcode= originalBarcode
#   }
#   #returns the mutated string
#
# }

recIntegrase_2019<-function(barcode,mu,alpha,type = "trit"){

  #barcodeLength=length(a)
  #fullBarcode = a
  #we want to edit only the portion that corresponds to
  #the currently active integrase

  a=strsplit(barcode,"")[[1]]


  for (c in 1:length(a)){
    #only mutate on unchanged elements
    if(a[c]=="1"){
      m_event<-runif(1,0,1) #random number
      if(m_event<mu[c]){
        #mutually exclusive events
        if(type=='trit'){
          trans_pr <-runif(1,0,1) #random number

          if(trans_pr<alpha[c]){
            a[c]="2"
          }else{
            a[c]="0"
          }
        }else if(type=='binary'){
          a[c]="0"
        }
      }
    }
  }

  fullBarcode= a

  return( paste(fullBarcode,sep="",collapse=""))


}


# Get the data.tree object ready starting from a phylo object (ground truth)
# returns this_node
initSimulation<-function(ground_truth, initial_barcode = "1111111111"){

  # name the internal nodes before converting
  ground_truth<-makeNodeLabel(ground_truth,prefix="")
  # convert to data.tree
  ground_node = as.Node(ground_truth)
  # check that things are correct then assign the barcode state to the first cell
  if(ground_node$isRoot)
    ground_node$barcode = initial_barcode

  return(ground_node)


}

# Main Function to simulate recording over an exisiting lineage
# returns ground_sim
simulateGroundTruth <- function(this_node,mu,alpha){

  if(this_node$isLeaf){
    return()
  }else{
    #set the new barcode by simulating a cell division recording
    this_node$children[[1]]$barcode = recIntegrase_2019(this_node$barcode,mu,alpha)
    this_node$children[[2]]$barcode = recIntegrase_2019(this_node$barcode,mu,alpha)

    # set the barcode as the name for the node
  #  this_node$children[[1]]$name =   this_node$children[[1]]$barcode
  #  this_node$children[[2]]$name =   this_node$children[[2]]$barcode


    simulateGroundTruth(this_node$children[[1]],mu,alpha)
    simulateGroundTruth(this_node$children[[2]],mu,alpha)

  }

  return(this_node)
}


# END MAIN SIMULATION




# GLOBAL parameters:
#Experimentally measured edit rates per site

# GLOBAL PARAMETERS
params_global = estim.params.global(estimG = 4, fil = paste("../",integrase_folder,"10mer_2019/editRate/allBarcodes.txt",sep=""))
mu = params_global[[1]]
alpha = params_global[[2]]


# New functions, specialized for working on an exisiting lineage

# pipeline:
# 1. get ground_truth tree from inspect.tree(i)
# 2. this_node = initSimulation()
# 3. ground_sim = simulateGroundTruth(this_node)
# 4. ground_phylo = convertToPhylo(ground_sim, ground_truth)
# 5. reconstruction = reconstructSimulation(ground_phylo)
# 6. dist[i] = RF.dist(ground_sim, reconstruction)

# Main function to get a distance from a simulated lineage
# Run it many times to get a histogram / average dist for THIS ground_truth
pipelineSim<-function(ground_truth,mu,alpha,return.tree = F,n_barcodes = 10){


  if(n_barcodes>10){
      mu = rep(mu, ceiling(n_barcodes/10))
      alpha = rep(alpha,ceiling(n_barcodes/10))
  }
  initial_barcode = paste(rep("1",n_barcodes),collapse="")

  if(return.tree ==F){
    d = initSimulation(ground_truth,initial_barcode) %>%
          simulateGroundTruth(mu = mu, alpha = alpha) %>%
            convertToPhylo(ground_truth = ground_truth) %>%
              reconstructSimulation(mu = mu, alpha = alpha)
  }else{
    #return the simulated tree
    d = initSimulation(ground_truth,initial_barcode) %>%
          simulateGroundTruth(mu = mu, alpha = alpha) %>%
            convertToPhylo(ground_truth = ground_truth)

  }

  return(d)
}

# parallel version of the simulation for a fixed lineage:
# For testing scanning parameters (edit rate, or barcodes)
parallelPipelineSim<-function(ground_truth =list(),n_barcodes = 10,edit_rate = 0.1, nsim =100){


  ground_reps_score = c();
  for(i in 1:nsim){
      ground_reps_score[i] = pipelineSim(ground_truth,rep(edit_rate,50),rep(0.5,50),n_barcodes = n_barcodes )
  }
  return(ground_reps_score)
}



# UPDATE Sep 16th 2019
# From the reconstruction scripts we create a memoirData structure with
# Using default parameters
# process.all.files(files_to_process) %>% reconstruct.all.lineages() -> memoirData
# memoirData is a dataframe with the ground truth trees and some statistics
# To perform simulation we just need to take the lineages and so we can jus extract that from the df
# memoirData includes only trees that have been filtered by quality / number of cells

# if MEMO mode is T , then we remove duplicated barcodes and reconstruct the lineage of "clones"
runAllMetrics<-function(id,memoirData = c(), n_repeats=50,n_barcodes=10,MEMO = F){
  print(paste( "processing tree", toString(id)))
  # in memoirData, lineages with 3 cells have score NA
  if(!MEMO){
    ground_truth = read.newick(text = as.character(memoirData$ground[id]))
  }else{
    ground_truth = read.newick(text = as.character(memoirData$MEMOtrees[id]))
  }

  # After converting the trees to MEMO mode, some are left with only 2 tips
  # We return NULL if the tree already was classified as NA or if they are <3 leaves after MEMO
  if(!is.na(memoirData$RF[id]) & length(ground_truth$tip.label)>3){


      if(n_barcodes>10){
          mu = rep(mu, ceiling(n_barcodes/10))
          alpha = rep(alpha,ceiling(n_barcodes/10))
      }

      d_membow = c()
      d_memoir = c()
      d_recursive = c()
      for(i in 1:n_repeats){
        this_node  =   initSimulation(ground_truth,initial_barcode = paste(rep("1",n_barcodes),collapse=""))
        ground_sim =   simulateGroundTruth(this_node,mu = mu, alpha = alpha)
        ground_phylo = convertToPhylo(ground_sim, ground_truth)
        manualTree = reconstructLineage(ground_phylo, mu,alpha,return_tree = T,clust.method="diana")

        d_memoir[i] = 1-RF.dist(ground_phylo,manualTree,normalize = T)

        #Try recursive reconstruction with the simulation data (update Sep 16th)
        manualTree_recursive = recursiveReconstruction(manualTree,mu = mu,alpha = alpha)
        d_recursive[i] =1- RF.dist( manualTree_recursive,ground_phylo,normalize = T)

        # New membow distance (no second round of reconstruction)
        d_membow[i] = reconstructMembow(ground_phylo,manualTree_recursive,mu, alpha)


      }
      return(list(d_membow,d_memoir,d_recursive , length(ground_phylo$tip.label)))



  }else{
    return(NULL)
  }


}

runAllMetrics_v0.1<-function(id, n_repeats =10,n_barcodes=10){
  #this simulation works at the tree level
  #it has to be automated to analyse all data

  # We first compute the global parameteres
  params_global = estim.params.global(estimG = 4, file = paste("../", integrase_folder,"10mer_2019/editRate/allBarcodes.txt",sep = ""))
  mu = params_global[[1]]
  alpha = params_global[[2]]


  if(n_barcodes>10){
      mu = rep(mu, ceiling(n_barcodes/10))
      alpha = rep(alpha,ceiling(n_barcodes/10))
  }


  res = inspect.tree(id, return.tree = T,plot.all=F,clust.method="diana")
  ground_truth = res[[8]]

  if(length(res)>0){
      d_membow = c()
      d_memoir = c()
      for(i in 1:n_repeats){
        this_node  =   initSimulation(ground_truth,initial_barcode = paste(rep("1",n_barcodes),collapse=""))
        ground_sim =   simulateGroundTruth(this_node,mu = mu, alpha = alpha)
        ground_phylo = convertToPhylo(ground_sim, ground_truth)
        manualTree = reconstructSimulation(ground_phylo, mu,alpha,return_tree = T)
        d_membow[i] = reconstructMembow(ground_phylo,manualTree,mu, alpha)
        d_memoir[i] = 1-RF.dist(ground_phylo,manualTree,normalize = T)
      }
      return(list(d_membow,d_memoir,length(ground_phylo$tip.label)))
  }else{
      return(NULL)
  }
}


#
# This function will get the barcode data from the simulated lineage
printBarcodes<-function(ground_sim){
  # Get all the leaves from the tree
  a<-ground_sim$leaves
  # For each leave, get the barcode
  sim_barcodes=c()
  for(i in 1:length(a)){
    sim_barcodes[i] = a[[i]]$barcode
  }
  return(sim_barcodes)
}

# Get the simulated lineage ready to be reconstructed
# returns ground_phylo
convertToPhylo<-function(ground_sim,ground_truth){

  # When we convert to phylo, the barcode information is lost
  # So we need to rename the leaves
  ground_phylo = as.phylo(ground_sim)
  # We extract leaves fromt he tree and save the barcodes as a list
  # The order is the same as in the phylo object
  sim_barcodes = printBarcodes(ground_sim)
  # So we can replace the names of the tips directly:
  ground_phylo$tip.label<-sim_barcodes

  cell_names = str_split(ground_truth$tip.label,"_",simplify=T)[,1]

  ground_phylo$tip.label<-paste(cell_names,ground_phylo$tip.label,sep="_")

  return(ground_phylo )

}

#returns RF dist for this simulated lineage

reconstructSimulation<-function(ground_phylo,mu,alpha,return_tree = F){

  #get the barcode and cell id (cell id is not neccessarily continuous numbers)
  barcodes = str_split(ground_phylo$tip.label,"_",simplify=T)[,2]
  cell_ids = str_split(ground_phylo$tip.label,"_",simplify=T)[,1]

  # translate the barcode data to the old notation uxr
  barcodes_urx = str_replace_all(barcodes, c("2" = "r", "1" = "u","0"="x"))
  barcodes_urx = str_replace_all(barcodes_urx, c("3" = "u"))

  #recontruct the tree
  matdist_=manualDistML_2(as.character(barcodes_urx),mu = mu ,alpha = alpha,nGen = 4 )

  row.names(matdist_)<- paste(cell_ids,barcodes,sep="_")
  colnames(matdist_)<- paste(cell_ids,barcodes,sep="_")

  hclust.tree = as.phylo(as.hclust( diana(as.dist(t(matdist_)))))

  #here ground_phylo is "the ground truth" which is a simulation but this is what we want
  d = 1- RF.dist(hclust.tree,ground_phylo,normalize =T)
  if(return_tree){
    return(hclust.tree)
  }else{
    return(d)
  }

}

# UPDATE this version only removes leaves and does not perform a second round of Reconstruction
# Only when the MEMO trees are less than min.memo.tree cells, we are going to reconstruct the barcodes again
# DIANA clustering might be biased when too many identical genotypes
# min.memo.tree = 3 means NO second round of reconstruction DEFAULT
reconstructMembow<-function(ground_phylo, manualTree,estimMu,estimAlpha,return_tree = F,
    min.memo.tree = 3,nGen = 4, rooted = T){

  unique.genotypes = unique(do.call(rbind,str_split(manualTree$tip.label,"_"))[,2])

  # check that tree has more that two unique leafs
  if(length(unique.genotypes)>2){

    new.tree = manualTree
    new.alive.tree = ground_phylo
        for(g in 1:length(unique.genotypes)){
            drop.leaves = ground_phylo$tip.label[grep(unique.genotypes[g],ground_phylo$tip.label)]
            #delete all except one

            if(length(drop.leaves)>1){
                for(i in 1:(length(drop.leaves)-1)){
                    #reconstruction
                    new.tree = drop.tip(new.tree,drop.leaves[i])
                    #ground truth:
                    new.alive.tree = drop.tip(new.alive.tree,drop.leaves[i])
                }
            }

        }

        # If there are only 3 unique genotypes there is no point in reconstructing the tree
        # since the RF distance is not well defined for n=3
        if(length(new.tree$tip.label)>3){
          # if after removing leaves we are left with only a few number of unique genotypes
          # we might want to re-reconstruct the barcodes using ward.D2 as DIANA does not work as well with too few cells
          if(length(new.alive.tree$tip.label)< min.memo.tree)
            new.tree = reconstructLineage(new.alive.tree,mu,alpha,return_tree=T,clust.method="ward.D2",nGen = nGen)

          this.score = 1-RF.dist(new.tree,new.alive.tree,normalize = T,rooted = rooted)
        }else{this.score = -1}
        #now we have trees with unique leaves
        # barcodes_ = new.tree$tip.label
        # barcodes.raw  =  substr(barcodes_,4,14)
        # barcodes.raw = str_replace_all(barcodes.raw, c("2" = "r", "1" = "u","0"="x"))
        #
        #
        # # using global parameters
        # matdist_2 = manualDistML_2(barcodes.raw,estimMu,estimAlpha,3)
        # colnames(matdist_2)<- substr(barcodes_,4,14); row.names(matdist_2)<- substr(barcodes_,4,14)
        # new.tree = as.phylo(as.hclust(diana(as.dist(t(matdist_2)))))
        #
        # new.alive.tree_ = new.alive.tree;
        # new.alive.tree_$tip.label = substr(new.alive.tree$tip.label,4,14)
        #
        # this.score = 1-RF.dist(new.alive.tree_,new.tree,normalize = T) #RF distance without number ID may20
        #
        # if(is.na(this.score)){
        #   #add a root to the tree, based on the initial state of the barcodes: 1111...1
        #   #then calcualte the distance using RF(x, root = T)
        #   this.score = 1-RF.dist(root(add.tips(new.tree, "1111111111", 4),"1111111111"),
        #       root(add.tips(new.alive.tree_, "1111111111", 4),"1111111111"),normalize = T,rooted = T)
        # }


  }else{this.score = -1; new.alive.tree = ground_phylo; new.tree = manualTree} #main IF end

  if(!return_tree){
    return(this.score)
  }else{
    return(list(new.alive.tree)) #,new.tree,this.score))
  }
}


randomControl <- function(id,memoirData, nrepeats = 10,MEMO = F){
  # This function will take a ground truth lineage given by the id
  # We are going to do random guess reconstruction which means just shuffling the leaves of the three
  # this, however assumes that we know the topology and we only want to put the barcodes in the right order.
  # Given that some topologies are easier to reconstruct than other, this assumption takes care of trivial reconstruction
  # since random guess will be probably as good as actual reconstruction

  # in memoirData, lineages with 3 cells have score NA
  if(!MEMO){
    ground_truth = read.newick(text = as.character(memoirData$ground[id]))
  }else{
    ground_truth = read.newick(text = as.character(memoirData$MEMOtrees[id]))
  }

  if(!is.na(memoirData$RF[id]) & length(ground_truth$tip.label)>3){



  # res = inspect.tree(id, return.tree = T,plot.all=F)
  # ground_truth = res[[8]]


    rand_dist=c()
    for(i in 1:nrepeats){
        ground_rand = ground_truth
        ground_rand$tip.label = sample(ground_rand$tip.label)
        rand_dist[i]=1-RF.dist(ground_truth,ground_rand,normalize=T)

    }
    return(rand_dist)

  }else{
    return(NULL)
  }

}

# Use this function to reproduce Fig 3-XX
makePlotAccuracy<-function(plot.only = T,nrepeats=10){

  #Run the whole analysis
  # Takes a few minutes
  # If already ran then set plot.only to TRUE (default)
  if(!plot.only){
    # Run simulations with 10mer
    res.sim<-lapply(1:100,runAllMetrics)
    sim_memoir<-do.call(rbind,lapply(res.sim,unlist))[,(nrepeats+1):(2*nrepeats)]
    sim_memoir<-ecdf(sim_memoir)
    #Run simulations with 20mer
    res.sim20mer<-lapply(1:100,runAllMetrics,n_barcodes =20)
    sim_memoir_20<-do.call(rbind,lapply(res.sim20mer,unlist))[,(nrepeats+1):(2*nrepeats)]
    sim_memoir_20<-ecdf(sim_memoir_20)
    #Random control
    res_rand = lapply(1:100,randomControl)
    rand_memoir = ecdf(unlist(res_rand))

    res.mat = runAllTrees()
    int_memoir = ecdf(as.numeric(levels(res.mat$memoir)[ as.numeric(res.mat$memoir) ]))
  }

  x11()
  plot(plot_range,sim_memoir(plot_range),main="reconstruction accuracy",ylab="Fraction of colonies",xlab="Reconstruction score (normalized RF)",lwd=2,type="l",col="grey")
  lines(plot_range,int_memoir(plot_range),lwd=2,col="black")
  lines(plot_range,rand_memoir(plot_range),lwd=2,col="red")
  lines(plot_range,sim_memoir_20(plot_range),lwd=2,col="blue")
  legend(0,0.5,legend = c("intMEMOIR", "sim 10mer","sim 20mer","random guees"), fill = c("black","grey","blue","red"))

}


# calculate entropy in barcodes and see whether that predicts reconstruction accuracy

#simple calculation of shannon entropy per site
# Returns the vector for entropies for each site
barcodeEntropy<-function( ground_truth, matrix_entropy = F ){

  barcodes = str_split(ground_truth$tip.label,"_",simplify=T)[,2]
  barcode_matrix<-do.call(rbind,lapply(lapply(barcodes,str_split,"",simplify=T),as.numeric))
  h = c()

  if(matrix_entropy){
    return(entropy(barcode_matrix))
  }else{

    for(i in 1:dim(barcode_matrix)[2]){
        h[i]=entropy::entropy(table(barcode_matrix[,i])/dim(barcode_matrix)[1])
      }
    return(h)
  }
}

# within a lineage, calculate the muatual information between the barcode states
barcodesMI<-function(ground_truth){
  barcodes = str_split(ground_truth$tip.label,"_",simplify=T)[,2]
  barcode_matrix<-do.call(rbind,lapply(lapply(barcodes,str_split,"",simplify=T),as.numeric))
  n_barcodes = dim(barcode_matrix)[2]

  # sum off diagonal values of the MI matrix, barcodes vs barcodes
  return(  sum(mutinformation(as.data.frame(barcode_matrix),method = 'emp')[upper.tri(matrix(0,n_barcodes,n_barcodes))]) )

}

barcodeEntropyString<-function(ground_truth){
  phylo_ground=read.newick( text = as.character( ground_truth))
  return(barcodeEntropy(phylo_ground) )

}

entropyAllTrees<-function(i){

      res = inspect.tree(i,return.tree = T, clust.method = "diana",plot.all = F)
      if(length(res)>0){
          ground_truth = res[[8]]

          h = sum(barcodeEntropy(ground_truth))
          memoir  = as.numeric(res[[4]] )
          membow = as.numeric( res[[5]])
          ncells = as.numeric( res[[3]])
          id = as.numeric( res[[2]])



      return(list(id,ncells,memoir,membow,h))
    }else{
      return(NULL)
    }

}

# Jan 2020
# # Entropy for simulated lineage:
stringToPhylo <- function(x){ read.newick(text = x) }
editRatePhylo <- function(x){
  tip_list<-lapply(x$tip.label,str_split,'_')
  barcode_list = do.call(rbind,lapply(tip_list,unlist))[,2]
  barcode_matrix = do.call(rbind,lapply(barcode_list,str_split,'',simplify = T))
  return(1-sum(barcode_matrix =='1')/ length(barcode_matrix))
 }

entropySimulation<-function(memoirData){

  phyloList<-lapply(memoirData$MEMOtrees, stringToPhylo)
  sim_phyloList = lapply(phyloList,pipelineSim,mu,alpha,return.tree = T)

  entropy_mat_sim = do.call(rbind,lapply(sim_phyloList,barcodeEntropy))
  entropy_mat_sim2 = do.call(rbind,lapply(sim_phyloList,barcodeEntropy,matrix_entropy=T))

  entropy_sum_sim =  rowSums(entropy_mat_sim)

  #calculate edit rate for simulated trees
  edit_rates = unlist(lapply(sim_phyloList,editRatePhylo))

  raw_entropy  = entropy_sum_sim
  raw_entropy_mat = entropy_mat_sim
  norm_entropy = entropy_sum_sim * edit_rates

  colony_size= unlist(lapply(memoirData$MEMOtrees,n_leaves))

  entropy_df = data.frame(colony_size = colony_size, entropy = norm_entropy,edit = edit_rates, raw_entropy = raw_entropy, raw_entropy2=raw_entropy_mat)

  entropy_df$rec_score = rep(NA,dim(entropy_df)[1])

  large_trees = sim_phyloList[entropy_df$colony_size>3]
  rec_scores_large_trees =unlist(lapply(large_trees,reconstructSimulation,mu,alpha))

  entropy_df$rec_score[entropy_df$colony_size>3] = rec_scores_large_trees

  return(entropy_df)

}

makeEntropyPlots<-function(entropy_res,entro_threshold= 2.3){


  max_memoir_ecdf = ecdf(entropy_res$rec_score[entropy_res$entropy <entro_threshold & entropy_res$rec_score >=0])
  low_memoir_ecdf = ecdf(entropy_res$rec_score[entropy_res$entropy >=entro_threshold & entropy_res$rec_score >=0])
  all_memoir_ecdf = ecdf( entropy_res$rec_score[ entropy_res$rec_score>=0] )


  ax = seq(0,1,0.001)
  plot(ax,all_memoir_ecdf(ax),lwd = 3,type= 'l',col = 'gray',
        ylab = "Fraction of colonies",xlab='RF score',main="Barcode entropy predict reconstruction accuracy",
        cex = 1.2,cex.lab = 1.7, cex.axis = 1.3)
  lines(ax,low_memoir_ecdf(ax),lwd = 3,type= 'l',col = 'red')
  lines(ax,max_memoir_ecdf(ax),lwd = 3,type= 'l',col = 'blue')

  legend(0.01,0.9, legend = c('All','L MEMOIR','H MEMOIR'), fill = c('gray','red','blue' ))
}
