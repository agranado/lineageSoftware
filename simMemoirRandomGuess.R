

#chat with Mark and grace:
#look for optimal edit rate 

#April 2nd
#Simulate a tree and then calcualte the distance to a randomized version of itself
#returns a vector with the distances of many comparisions (many randomizations)
#Update March 30th, function works well 
simMemoirRandomGuess<-function(nGen,mu,alpha,barcodeLength,nRand){
  #load necessary libraries and functions
  library(phangorn)
  source("simulation2.R")
  pathName= "/Users/alejandrog/MEGA/Caltech/trees/simulation/"
  #Estimate ML trees:
  #ML trees does not significantly increase the performance of lineage reconstruction
  #It does, actually, performa way worse than normal upgma/distance methods. 
  #Moreover, it is computationally expensive, especially when performing parameter analysis. 
  #if you want to estimate ML edit the following line:
  estimML=FALSE
  
  
  
  #clear the variable (since it behaves as global)
  if(!exists("firstCell")){
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
    divideCellRecursive2(firstCell,mu,alpha)
  }
  
  #prints only the barcodes for all leaves
  print(firstCell,"barcode")
  print("Tree simulation completed")
  #save to file as newick tree
  #save the length of branches plus the ID (which so far is a number)
  newickTree<-ToNewick(firstCell)
  #Generate unique ID for writing file to disc (we'll erase it later)
  fileID = toString(runif(1))
  firstCellFile = paste(pathName,"trees/firstCell",fileID,".nwk",sep="")
  write(newickTree,file=firstCellFile)
  #load the tree from the file as a tree structure. 
  trueTree<-read.tree(file=firstCellFile)
  
  
  
  #file is now deleted
  print("True tree read")
  # plot(trueTree,main=paste("True tree ",sep=""))
  randDist = array()
  for (r in 1:nRand){
    randomTrue<- randomTrueTree(trueTree)
    randDist[r]=treedist(randomTrue,trueTree)
  }
  
  
  
  
  system(paste("rm ",firstCellFile,sep=""))
  
  return(randDist)
}   
