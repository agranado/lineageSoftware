# Aug 14th 2019
# Integrase paper writing
#
# Simulation is restricted to toplogies that are observed experimentally
# Goal: Establish a baseline for reconstruction accuracy


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

# old functions (simplified for this particular situation)
mutationScratchpad_2019 <- function(barcode,mu,alpha,type='trit',recType="integrase",cInts=1,tInts=1,chr.acc=c()){

  originalBarcode=strsplit(barcode,"")[[1]]
  #now a is an array of char elements representing the scratchpad
  m_event<-runif(1,0,1)
  #mutation happens (rate = mutations/generation)
  #for each element in the array, ask if mutation happens
  #with constant independent rate (mutation in [1] does not affect [2])
  if(recType=="integrase"){
    mutatedBarcode=recIntegrase(originalBarcode,mu=mu,alpha=alpha,type=type,currentInts=cInts,totalInts=tInts)
  }else if(recType=="epimemoir"){ #epimemoir
    mutatedBarcode=recEpiMemoir(originalBarcode,mu=mu,alpha=alpha,type=type,totalInts=tInts,chr.acc=chr.acc)
  }else { #no mutation
    mutatedBarcode= originalBarcode
  }
  #returns the mutated string

}

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
pipelineSim<-function(ground_truth,mu,alpha){

  d = initSimulation(ground_truth) %>%
        simulateGroundTruth(mu = mu, alpha = alpha) %>%
          convertToPhylo(ground_truth = ground_truth) %>%
            reconstructSimulation(mu = mu, alpha = alpha)
  return(d)
}


# UPDATE Sep 16th 2019
# From the reconstruction scripts we create a memoirData structure with
# Using default parameters
# process.all.files(files_to_process) %>% reconstruct.all.lineages() -> memoirData
# memoirData is a dataframe with the ground truth trees and some statistics
# To perform simulation we just need to take the lineages and so we can jus extract that from the df
# memoirData includes only trees that have been filtered by quality / number of cells

# if MEMO mode is T , then we remove duplicated barcodes and reconstruct the lineage of "clones"
runAllMetrics<-function(id,memoirData = c(), n_repeats=10,n_barcodes=10,MEMO = F){
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
reconstructMembow<-function(ground_phylo, manualTree,estimMu,estimAlpha,return_tree = F,min.memo.tree = 3,nGen = 4){

  unique.genotypes = unique(substr(manualTree$tip.label,4,14))

  # check that tree has more that two unique leafs
  if(length(unique.genotypes)>2){

    new.tree = manualTree
    new.alive.tree = ground_phylo
        for(g in 1:length(unique.genotypes)){
            drop.leaves = manualTree$tip.label[grep(unique.genotypes[g],manualTree$tip.label)]
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

          this.score = 1-RF.dist(new.tree,new.alive.tree,normalize = T)
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
    return(list(new.alive.tree,new.tree,this.score))
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
barcodeEntropy<-function( ground_truth ){

  barcodes = str_split(ground_truth$tip.label,"_",simplify=T)[,2]
  barcode_matrix<-do.call(rbind,lapply(lapply(barcodes,str_split,"",simplify=T),as.numeric))

  for(i in 1:dim(barcode_matrix)[2]){
      h[i]=entropy(table(barcode_matrix[,i])/dim(barcode_matrix)[1])
  }

  return(h)
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
