
#commands that execute this function in parallel

#200 repeats, using the indicated parameters,
# the list of methods contains the distance metrics that we want to test in the stringdist function
# >  results= foreach(i=1:200) %dopar% simMemoirStrdist(nGen,mu,alpha,barcodeLength,methods)

#convert list to matrix
# >  results.matrix=do.call(rbind,results)
# >  apply(results.matrix,2,mean)

#Arguments of the function:
#barcodeLength<-16
#nGen =5
#mu= (1/barcodeLength)*3;
#mu=1/5
#alpha= 2/3;

library(phangorn)
library(stringdist)
library(doParallel)
source("simulation2.R")

compareDist <- function(simulationType='trit',nGen=3,mu=0.4,alpha_=2/3,barcodeLength=6,nRepeats=20){


  results= foreach(i=1:nRepeats) %dopar% simMemoirStrdist_2020(nGen=nGen,mu=mu,alpha=alpha_,barcodeLength=barcodeLength,simulationType=simulationType,
                                                              write.newick = T)
  results.matrix=do.call(rbind,results)
  #Optional when only interested in the mean
  #apply(results.matrix,2,mean)
  return(results.matrix)
}

#functions to use the proportion of perfect trees as the measure.
#May 8
eq.zero<-function(r,x){sum(r[,x]==0)}





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

#Mar 2nd 2020
# Let's simplify this function
# the actual simulation step is not the most time consuming
# conversion toPhylo takes way longer than the simulation itself
# So let's get rid of useless stuff

simMemoirStrdist_2020<-function(nGen,mu,alpha,barcodeLength,simulationType,write.newick = F){
  #load necessary libraries and functions
  #detection of OS

  #update (trying to create a single branch that works in AWS and in my laptop)
  os=system("cat ../os.txt",intern = T)
  if(os=="mac"){ #local Alejandro's laptop
    pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"
    pathName2="/Users/alejandrog/MEGA/Caltech/trees/simulation"
  }else if(os=="linux"){ #AWS server or any other linux machine (path is for AWS)
    pathName="/home/ubuntu/alejandrog/Caltech/lineage/"
    pathName2="/home/ubuntu/alejandrog/Caltech/lineage"

  }
  #clear the variable (since it behaves as global)
  if(exists("firstCell")){
    rm(firstCell)
  }
  #create intialize the tree using 1 as the ID
  #Node is a function from data.tree
  firstCell <- Node$new("1")

  #number of letters (this will change for simulation)
  #A goal is to understand how the lenght of the pattern affects the ability to reconstruc lineages

  #intialitze the barcode using "u->unchanged
  firstCell$barcode <-paste(rep("u",barcodeLength),collapse="")

  #all variables of the data.tree structure are global
  #any pointer to a node, points to the real tree.

  for (g in 1:nGen){
    #this function simulates one generation of the tree
    divideCellRecursive2(firstCell,mu,alpha,type=simulationType)
  }

  # takes a long time ( for G = 12, 75 min)
  ground_phylo = as.phylo(firstCell)

  barcodes = printBarcodes(firstCell)
  ground_phylo$tip.label<-paste(as.character(1:length(barcodes)),"_", barcodes, sep ="")
  # rename the phylo OBJECT
  #manualdist from MLfunctinos.R
  mu_array = rep(mu,barcodeLength)
  alpha_array = rep(alpha,barcodeLength)

  aa = reconstructSimulation(ground_phylo,mu_array,alpha_array)

  allDistances = c()
  allDistances[1] = 0
  allDistances[2] = aa
  if(write.newick){
    write.tree(ground_phylo,file = paste(pathName,"2020trees/nG_",as.character(nGen),"_NBC_",as.character(barcodeLength),"mu_",as.character(mu),"_",as.character(runif(1)),".nwk",sep =""))
  }

    return(allDistances)
}



reconstructSimulation<-function(ground_phylo,mu,alpha,return_tree = F){

  #get the barcode and cell id (cell id is not neccessarily continuous numbers)
  barcodes = str_split(ground_phylo$tip.label,"_",simplify=T)[,2]
  cell_ids = str_split(ground_phylo$tip.label,"_",simplify=T)[,1]

  # translate the barcode data to the old notation uxr
  #barcodes_urx = str_replace_all(barcodes, c("2" = "r", "1" = "u","0"="x"))
  barcodes_urx = barcodes

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




























#April 8th
#Test stringdistance measures using the stringdist R library
#use the same format as before but testing different methods included in the stringdist function
simMemoirStrdist<-function(nGen,mu,alpha,barcodeLength,methods,simulationType){
  #load necessary libraries and functions
  #detection of OS

  #update (trying to create a single branch that works in AWS and in my laptop)
  os=system("cat ../os.txt",intern = T)
  if(os=="mac"){ #local Alejandro's laptop
    pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"
    pathName2="/Users/alejandrog/MEGA/Caltech/trees/simulation"
  }else if(os=="linux"){ #AWS server or any other linux machine (path is for AWS)
    pathName="/home/ubuntu/alejandrog/Caltech/lineage/"
    pathName2="/home/ubuntu/alejandrog/Caltech/lineage"

  }
  #clear the variable (since it behaves as global)
  if(exists("firstCell")){
    rm(firstCell)
  }
  #create intialize the tree using 1 as the ID
  #Node is a function from data.tree
  firstCell <- Node$new("1")

  #number of letters (this will change for simulation)
  #A goal is to understand how the lenght of the pattern affects the ability to reconstruc lineages

  #intialitze the barcode using "u->unchanged
  firstCell$barcode <-paste(rep("u",barcodeLength),collapse="")

  #all variables of the data.tree structure are global
  #any pointer to a node, points to the real tree.

  for (g in 1:nGen){
    #this function simulates one generation of the tree
    divideCellRecursive2(firstCell,mu,alpha,type=simulationType)
  }

  #prints only the barcodes for all leaves
  #  print(firstCell,"barcode")
  #  print("Tree simulation completed")
  #save to file as newick tree
  #save the length of branches plus the ID (which so far is a number)
  newickTree<-ToNewick(firstCell)
  #Generate unique ID for writing file to disc (we'll erase it later)
  fileID = toString(runif(1))

  #firstCellFile = paste(pathName,"trees/firstCell",fileID,".nwk",sep="")

  firstCellFile =tempfile("trees/firstCell",tmpdir = pathName2)
  firstCellFile =paste(firstCellFile,fileID,".nwk",sep="")


  write(newickTree,file=firstCellFile)
  #load the tree from the file as a tree structure.
  trueTree<-read.tree(file=firstCellFile)

  #file is now deleted
  #  print("True tree read")
  # plot(trueTree,main=paste("True tree ",sep=""))

  #get the sequences from the simulated tree + names
  barcodes<-firstCell$Get("barcode")

  #now we have the patters for ALL cells, but we need only the leaves, since that is what
  #we are going to use for reconstruction.
  #The way the tree was built, only the last 2^g cells are leaves; where g is the number of generations
  #take the number ID for the leaves.
  leavesID<-(length(barcodes)-2^nGen+1):length(barcodes)
  #grab those cells from the tree

  barcodeLeaves = array()
  namesLeaves=array()
  for (l in 1:length(leavesID)){
    barcodeLeaves[l] <-barcodes[names(barcodes)==leavesID[l]]
    namesLeaves[l] <-names(barcodes)[names(barcodes)==leavesID[l]]
  }
  names(barcodeLeaves)<-namesLeaves
  #now barcodeLeaves has all the leaves of the tree, with their original ID from the data.tree structure.
  #create Fasta file using only the leaves of the tree (n= 2^g)
  fastaBarcodes<-convertSimToFasta(barcodeLeaves)
  #convert name of variable to string

  #3  varName<-deparse(substitute(firstCell))

  fasID =toString(runif(1))
  #fasIN <-paste(pathName,"fasta/firstCell_",fasID,".fas",sep="")
  fasIN =tempfile("fasta/firstCell",tmpdir = pathName2)
  fasIN = paste(fasIN,fasID,".fas",sep="")

  write(fastaBarcodes,file=fasIN)
  #  print("writting fasta file, simulated tree")
  #this is the format:
  #>1_uuuuuu
  #uuuuuu
  #>2_uxuxuu
  #uxuxuu
  #>4_uxuxxc
  #uxuxxc
  #>8_uxuxxc
  #uxuxxc

  #we can also convert the data.tree format to igraph format
  #igraph has more algorithms and plotting function, so it's worth.
  #LIBRARY
  #library(igraph)
  #as.igraph.Node converts data.tree to igraph object
  #  firstCell_igraph<-as.igraph.Node(firstCell)
  #rename all vertices in the network using the barcodes.
  # V(firstCell_igraph)$name<-barcodes
  #the object has now all info.
  #REMEBER: this is the REAL tree, still we need to perform alignment and reconstruction


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  #### FROM here comes the alignment and reconstruction
  #previously saved fasta file
  #LIBRARY\


  #this a new type of object, from the package phangorn "phyDat"

  #we can apply lineage reconstruction here so later we can compare with the real tree.
  #before loading the fasta file we need to convert the characters to DNA, so that it is compatible
  #sed.return ==1 if replacement went through (Apr 25)
  sed.return<-convertMemoirToDNA(fasIN)
  #now we can use the phyDat
  if(sed.return){
    memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")
    #for distance based trees
    #from the phangorn tutorial pdf:
  }
  #Apr 8th: this is where the distance comes into place

  #Apr 8th:
  #we calculate string distances only for the leaves ( which is the data we actually get)

  allDistances = array()
  for(m in 1:length(methods)){
    stringTree= upgma(stringdistmatrix(barcodeLeaves,method=methods[m]))
    stringTree$tip.label<-trueTree$tip.label
    allDistances[m]= RF.dist(stringTree,trueTree)
  }

  #control against default dist.ml function + UPGMA, which so far is the best method
  if(sed.return==1){
    dm<-dist.ml(memoirfas)
    dm.ham=dist.hamming(memoirfas)
    treeUPGMA<-upgma(dm)
    treeUPGMA_noSeq<-removeSeqLabel(treeUPGMA)
    allDistances[m+1]= RF.dist(treeUPGMA_noSeq,trueTree)
  }else{
    allDistances[m+1]= NaN
  }

  #Apr 9th
  #Manual distance calculation (v beta1.0)
    #  matdist = manualDist(barcodeLeaves,mu,alpha,nGen)
    #  manualTree =upgma(as.dist(t(matdist)))
    #  manualTree$tip.label= treeUPGMA$tip.label

    #  allDistances[m+2]= RF.dist(removeSeqLabel(manualTree),trueTree)

    allDistances[m+2]=0

    #manualdist from MLfunctinos.R

    matdist_ = manualDistML(barcodeLeaves,mu,alpha,nGen)
    manualTree_ =upgma(as.dist(t(matdist_)))
    manualTree_$tip.label= treeUPGMA$tip.label

    allDistances[m+3]= RF.dist(removeSeqLabel(manualTree_),trueTree)

  #  print("All distances calcualted")


  #delete files

  system(paste("rm ",firstCellFile,sep=""))
  if(sed.return){
    system(paste("rm ",paste(fasIN,".bak",sep=""),sep=""))
    system(paste("rm ",fasIN,sep=""))
  }
  #system(paste("rm ",fasIN,".bak",sep=""))


  return(allDistances)
}
  #END of simulation function
#Apr 9th

Pr_edit <- function (nGen,mu,alpha){
  #probability that at nGen generations there is a mutation
  cumSum =0
  for (nG in 1:nGen){
    cumSum = cumSum + mu * (1-mu)^(nG-1) * alpha
  }
  return(cumSum)
}





#Apr 23
#GIT update for analyzing performance of reconstruction method with simulated data.
#these updates happened before the comparison between binary and tri simulations
# code section for ATOM execution and examples
#registerDoParallel()
#results<-compareDist()
