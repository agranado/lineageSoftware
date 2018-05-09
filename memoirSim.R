#simTreeDist<-function(barcodeLength,nGen,mu,alpha){
      #load all necessary function from the library


      #parameters come from outside the script
      #barcodeLength<-16
      #nGen =5
      #mu= (1/barcodeLength)*3;
      #mu=1/5
      #alpha= 2/3;

      source("/Users/alejandrog/MEGA/Caltech/trees/simulation/simulation2.R")
      #clear the variable (since it behaves as global)
      rm(firstCell)
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
      write(newickTree,file="trees/firstCell.nwk")
      #load the tree from the file as a tree structure. 
      trueTree<-read.tree(file="trees/firstCell.nwk")
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
      varName<-deparse(substitute(firstCell))
      fasIN <-paste("fasta/",varName,".fas",sep="")
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
      library(igraph)
      #as.igraph.Node converts data.tree to igraph object
      firstCell_igraph<-as.igraph.Node(firstCell)
      #rename all vertices in the network using the barcodes. 
      V(firstCell_igraph)$name<-barcodes
      #the object has now all info. 
      #REMEBER: this is the REAL tree, still we need to perform alignment and reconstruction
      
      
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      
      #### FROM here comes the alignment and reconstruction 
      #previously saved fasta file
      #LIBRARY\
      
      library(phangorn)
      #this a new type of object, from the package phangorn "phyDat"
      
      #we can apply lineage reconstruction here so later we can compare with the real tree. 
      #before loading the fasta file we need to convert the characters to DNA, so that it is compatible
      convertMemoirToDNA(fasIN)
      #now we can use the phyDat 
      memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")
      #for distance based trees
      #from the phangorn tutorial pdf: 
      dm<-dist.ml(memoirfas)
      dm.ham=dist.hamming(memoirfas)
      print("distance calculated")
      treeUPGMA<-upgma(dm)
      treeUGPMA.ham <-upgma(dm.ham)
      
      print("UPGMA ready")
      treeNJ<-NJ(dm)
      
      #parsimony score:
      parsimony(treeUPGMA,memoirfas)
      print("parsimony ready")
      #create max likelihood object 
      fitGTR=pml(treeUPGMA,data=memoirfas)#nothing happens to topology here
      print("fit jukes cantor ready")
      #update the model for ML optimization
      #fitGTR is a max Lik tree, with base frequencies. 
      
      #from simulations and modelTest (March 29) G+I do not make big difference
      #Moreover, F81 & HKY seem to be the best models in terms of AIC. (GTR would unnecessarily overcomplicate things)
      #  fitGTR_<-update(fit,k=4,inv=0.1)#from JC to +G+I ## March 29: we no longer add G I 
      
      #print("GTR updated")
      #fitGTR<-optim.pml(fitGTR,model="GTR",optInv=T,optGamma=T,rearragement="NNI",control=pml.control(trace=0),optNni=T)
     
      #HEre is where the ML optimization actually happens
      #fitGTR2<-optim.pml(fitGTR,model="GTR",optInv=T,control=pml.control(trace=0),rearrangement = "stochastic")
	   #March 29
      fit.HKY <- optim.pml(fitGTR, optNni = TRUE, optEdge = TRUE, model = "HKY",optInv=T,optGamma=T)
      fit.F81 <- optim.pml(fitGTR, optNni = TRUE, optEdge = TRUE, model = "F81",optInv=T,optGamma=T)
     # fitGTR <- optim.pml(fitGTR_, model="GTR", optInv=TRUE, optGamma=TRUE,rearrangement = "ratchet", control = pml.control(trace = 0))
      
      print("optim.pml..done")
	  
      #calcualte distance between trees:
      #At this point, we have the real tree, which is a simulation 
      #we also have the tree from the alignment (UPGMA), we can calculate the
      #distance between the trees since both are object from the phylo class
      #we first need to make sure the labels are the same 
      #duplicate the UPGMA tree


      #Here we calcualte the distance between 3 diffrent topologies and the real tree
      
      treeUPGMA_noSeq<-removeSeqLabel(treeUPGMA)
      treeUPGMA.ham_noSeq=removeSeqLabel(treeUGPMA.ham)
      treeNJ_noSeq<-removeSeqLabel(treeNJ)
      #use the maximum likelihood model to calculate the distances. 
      #here we can take the basefrequencies/ rate matrix to create the distances:
      #treeGTR_noSeq<- removeSeqLabel(upgma(dist.ml(memoirfas,bf=fitGTR$bf,Q=fitGTR$Q)))
      #specify only the base frequencies
      #OLD #treeGTR_noSeq<- removeSeqLabel(upgma(dist.ml(memoirfas,bf=fitGTR$bf)))
      
      #take the tree from the fit (we allowed for topology rearrangements)
      #March 29: test two different models (both look good from Model test)
      tree.HKY_noSeq<- removeSeqLabel(fit.HKY$tree)
      tree.F81_noSeq<- removeSeqLabel(fit.F81$tree)
      # NOTE: sometimes Q matrix creates Inf values for distances which crashes the clustering method in UPGMA
      
      #March 28:optim parsimony based on the UPGMA tree
      #this will look for different topologies that have max parsimony
      treePars_noSeq<-removeSeqLabel(optim.parsimony(treeUPGMA,memoirfas))
      
      treePars_noSeq2<-removeSeqLabel(optim.parsimony(treeUPGMA,memoirfas,method="sankoff"))
      #Calculate distance between the trueTree and different methods
      #treedist(treeUPGMA_noSeq,trueTree)
      rfdist=RF.dist(treeUPGMA_noSeq,trueTree)
      rfdist.ham=RF.dist(treeUPGMA.ham_noSeq,trueTree)
      tfdist_NJ=RF.dist(treeNJ_noSeq,trueTree)
      #ML trees
      rfdist_HKY=RF.dist(tree.HKY_noSeq,trueTree)
      rfdist_F81=RF.dist(tree.F81_noSeq,trueTree)

      #not used rfML.HKY = RF.dist(removeSeqLabel(fit.HKY$tree),trueTree)
      # not used rfML.F81 = RF.dist(removeSeqLabel(fit.F81$tree),trueTree)
      rfPars= RF.dist(treePars_noSeq,trueTree)
      rfPars2= RF.dist(treePars_noSeq2,trueTree)
      
      matdist = manualDist(barcodeLeaves,mu,alpha,nGen)
      manualTree =upgma(as.dist(t(matdist)))
      manualTree$tip.label= treeUPGMA$tip.label
      
      
      #plot the trees
      allDistances =c(rfdist,tfdist_NJ,rfdist_HKY,rfdist_F81,rfPars,rfPars2)
      print("All distances calcualted")
    #PLOT
     # x11();
     #  par(mfrow=c(1,2))
     #  plot(treeUPGMA_noSeq,main=toString(rfdist))
     # plot(trueTree)
 #     return(list(trueTree,treeUPGMA_noSeq,rfdist))     
#}