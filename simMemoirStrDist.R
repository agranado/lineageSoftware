
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

compareDist <- function(nGen=3,mu=0.4,alpha=2/3,barcodeLength=6,nRepeats=20,methods=c('osa','lv','dl','hamming','lcs','qgram','cosine','jaccard','jw','soundex')){
    library(doParallel)

    results= foreach(i=1:nRepeats) %dopar% simMemoirStrdist(nGen,mu,alpha,barcodeLength,methods)
    results.matrix=do.call(rbind,results)
    #Optional when only interested in the mean
    #apply(results.matrix,2,mean)
}



#April 8th
#Test stringdistance measures using the stringdist R library
#use the same format as before but testing different methods included in the stringdist function
simMemoirStrdist<-function(nGen,mu,alpha,barcodeLength,methods,simulationType='trit'){
  #load necessary libraries and functions
  library(phangorn)
  library(stringdist)
  source("/Users/alejandrog/MEGA/Caltech/trees/GIT/simulation2.R")
  pathName="/Users/alejandrog/MEGA/Caltech/trees/simulation/"


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
    divideCellRecursive2(firstCell,mu,alpha,type=simulationType)
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
  fasIN <-paste(pathName,"fasta/firstCell_",fasID,".fas",sep="")


  write(fastaBarcodes,file=fasIN)
  print("writting fasta file, simulated tree")
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
  convertMemoirToDNA(fasIN)
  #now we can use the phyDat
  memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")
  #for distance based trees
  #from the phangorn tutorial pdf:

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
  dm<-dist.ml(memoirfas)
  dm.ham=dist.hamming(memoirfas)
  treeUPGMA<-upgma(dm)
  treeUPGMA_noSeq<-removeSeqLabel(treeUPGMA)
  allDistances[m+1]= RF.dist(treeUPGMA_noSeq,trueTree)


  #Apr 9th
  #Manual distance calculation (v beta1.0)
  matdist = manualDist(barcodeLeaves,mu,alpha,nGen)
  manualTree =upgma(as.dist(t(matdist)))
  manualTree$tip.label= treeUPGMA$tip.label

  allDistances[m+2]= RF.dist(removeSeqLabel(manualTree),trueTree)

  print("All distances calcualted")


  #delete files

  system(paste("rm ",firstCellFile,sep=""))
  system(paste("rm ",fasIN,sep=""))


  return(allDistances)
}

#Apr 9th

Pr_edit <- function (nGen,mu,alpha){
    #probability that at nGen generations there is a mutation
  cumSum =0
    for (nG in 1:nGen){
      cumSum = cumSum + mu * (1-mu)^(nG-1) * alpha
    }
  return(cumSum)
}

#calcualtion of distance based on the properties of the system.
#should work only for this implementation of memoir
manualDist <- function(barcodeLeaves,mu,alpha,nGen){
   alphabet = array()
   alphabet[1]= "u"
   alphabet[2]= "r"
   alphabet[3]= "x"
   #mu & alpha indicate an explicit probabilistic model
   #Probability of no mutation for nGen cell divisions (mu is rate of edit per cell division)

   #transition probabilities
   #beacuse transitions only happen from u ->  then this is a 1-D vector
   Tran.pr = c(1, alpha,1-alpha)


   #NULL MODEL: probability of observing sites as independent events:
   Pr = array()
   Pr[1] = (1-mu)^nGen
   #probability of nu mutation during nGen-1 devisions and then a mutation in generation nGen times Pr(alphabet[2])
   #then we use the choose to correct for all the order in which this could have happened, (we need still further correction to
   #to account for irreversibility)
   #Pr[2] = choose(nGen,nGen-1)*( 1-mu)^(nGen-1)*mu*alpha
   #corrected for irreversibility:
   Pr[2] = Pr_edit(nGen,mu,alpha)

   #same as before but using (1-alpha)
   #Pr[3] = choose(nGen,nGen-1)*(1-mu)^(nGen-1)*mu*(1-alpha)
   Pr[3] = Pr_edit(nGen,mu,1-alpha)

   PrMatrix  = array(0,dim=c(length(Pr),length(Pr)))
   #calcualte probabilistic model:
   #this just means the pr that those two sites have those characters by random. This is the expected pr for each site
   #it assummes independence but does not tell you how likely they are to come from a common ancestor
   for (p1 in 1:length(alphabet)){
     for (p2 in 1:length(alphabet)){
        PrMatrix[p1,p2] = Pr[p1] * Pr[p2]
     }
   }

   #weights for number of sustitutions
   equalU =0
   oneSust = 1
   twoSust = 2

   nBarcodes = length(barcodeLeaves)
   barcodeLength=nchar(barcodeLeaves[1])
   distMat= array(0,dim =c(nBarcodes,nBarcodes))
   ratioMat = array(0,dim=c(nBarcodes,nBarcodes))
   productMat = array(0,dim=c(nBarcodes,nBarcodes))

    #go through all the elements in the barcode array
    for (i in 1:(nBarcodes-1)){
      for (j in (i+1):nBarcodes){
        barcodeArray1 =strsplit(barcodeLeaves[i],"")[[1]]
        barcodeArray2 =strsplit(barcodeLeaves[j],"")[[1]]
        distSum = 1
        ratio.sum =0
        ratio.product=1
        #for each pairwise comparison
        for (s in 1:barcodeLength){
            #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
            Pr.sust = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])] #this is just the product of both Pr
            Pr.sust.inv = 1/Pr.sust
            #both characters are the same (independently of their identity)
            if(barcodeArray1[s]==barcodeArray2[s]){
                #calcualte probability of both sites
                distSum = distSum - equalU
                #Pr_sust
                #probabilities under assumption of sister cells
                if(barcodeArray1[s]=="u"){
                  #probability of u in the previous generations times pr(no sust) * pr(no sust)
                  Pr_sister=(1-mu)^(nGen-1) * (1-mu)^ 2
                }else if(barcodeArray1[s]=="r"){
                  #if both are r : probability that ancestor is r + pr_a (u) * Pr(sust) ^2
                  # Pr(r_{t-1}) + Pr(u_{t-1}) * Pr(u->r)^2 = Pr(u->r_{t},u->r_{t} | u_{t-1})
                  Pr_sister=Pr_edit(nGen-1,mu,alpha) + (1-mu)^(nGen-1) * (mu *alpha)^2
                }else if(barcodeArray1[s]=="x"){
                  #if both are x : probability that ancestor is r + pr_a (u) * Pr(sust) ^2
                  Pr_sister=Pr_edit(nGen-1,mu,1-alpha) + (1-mu)^(nGen-1) * (mu * (1-alpha))^2
                }
                  #ratio between sister probability and random probability
                  ratio.sum = ratio.sum + Pr_sister/Pr.sust
                  ratio.product = ratio.product* Pr_sister/Pr.sust
            }else{
              #characteres have different sites.
              #is any of the character a u
              b = grepl(alphabet[1],c(barcodeArray1[s],barcodeArray2[s]))

              #one of them is u
              if(length(which(b==FALSE))==1){
                   distSum = distSum + oneSust *Pr.sust
                   #are there any r
                   c = grepl(alphabet[2],c(barcodeArray1[s],barcodeArray2[s]))
                   #there is one r
                   if(length(which(c==TRUE))==1){
                      Pr_sister=(1-mu)^(nGen-1) * (1-mu) * mu * Tran.pr[2]

                   }else {
                     #it is an x:
                     # Pr_sister = Pr(u_{t-1}) * Pr(u->u_{t}) * Pr(u->r_{t})
                     Pr_sister=(1-mu)^(nGen-1) * (1-mu) * mu * Tran.pr[3]

                   }
                   ratio.sum = ratio.sum + Pr_sister/Pr.sust
              #none of them is u AND they are different
              }else if(length(which(b==FALSE))==2){
                  distSum = distSum + twoSust *Pr.sust

                  #Pr_sister= Pr(r,x | u_{t-1})
                  Pr_sister = (1-mu)^(nGen-1) * mu^2 * Tran.pr[2] * Tran.pr[3]
                  ratio.sum  = ratio.sum + Pr_sister/Pr.sust
              }
              #sister probabilities:
            }
            ratio.product = ratio.product* Pr_sister/Pr.sust
        }
    distMat[i,j]= distSum
    ratioMat[i,j]=1/ratio.sum *distSum
    productMat[i,j] = 1/ratio.product
    }
  }
return(ratioMat)
}






runThisScript <- function (){
    #run memoirSim.R once and plot both trees.
    x11()
    par(mfrow=c(1,2))
    source("/Users/alejandrog/MEGA/Caltech/trees/simulation/memoirSim.R")
    ratioMat = manualDist(barcodeLeaves,mu,alpha,nGen );ratioTree = upgma(as.dist(t(ratioMat))); ratioTree$tip.label= treeUPGMA$tip.label;

    plot(ratioTree,main=paste("Manual tree",toString(RF.dist(trueTree,removeSeqLabel(ratioTree))) ))
    plot(treeUPGMA,main=paste("UPGMA dist tree",toString(RF.dist(trueTree,removeSeqLabel(treeUPGMA))) ))

    #execute multiple trees to get statistics
    # results= foreach(i=1:200) %dopar% simMemoirStrdist(3,0.3,alpha,barcodeLength,methods); results.matrix=do.call(rbind,results); meanDist =apply(results.matrix,2,mean)
}




#Apr 23
#GIT update for analyzing performance of reconstruction method with simulated data. 
#these updates happened before the comparison between binary and tri simulations
# code section for ATOM execution and examples
registerDoParallel()
results<-compareDist()
