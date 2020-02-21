
# Jan 13th.
# Goal: let's keep everything the same (trit, integrase, mu, alpha, reconstruction, etc)
# Just change the recording mode and add a history component to it


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
