
filePath="/Users/alejandrog/MEGA/Caltech/trees/GraceData/"
treeFileNames = system(paste("ls ", filePath," | grep .csv",sep=""), intern=T)
treeNames = do.call(rbind,strsplit(treeFileNames,"\\."))[,1]

#[1] "Pos21"     "Pos31_1"   "Pos31_2"   "Pos32"     "Pos33"     "Pos34"     "Pos35"    
#[8] "Pos37"     "Pos39"     "Pos40"     "pos1_cont"

treeGen= c(0,0,2,0,0,0,0,0,0,0,0)



#read barcodes 
alphabet= c("u","r","x")
letterFreq=array(0,dim=c(length(alphabet),dim(unique.BC.matrix)[1]));
for (a in 1:length(alphabet)){
  for (d in 1:dim(unique.BC.matrix)[1]){ letterFreq[a,d]=length(which(unique.BC.matrix[d,]==alphabet[a]));}
#average number of unedited barcodes across all cells (we consider only unique barcode profiles)

}
meanFreq=apply(letterFreq,1,mean)




