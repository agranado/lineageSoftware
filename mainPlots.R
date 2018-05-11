#May 9th
#main plots that might be used in the paper
#this script should load the object muVariation or it should be executed afeter running bitVStrit.R


# plotting ----------------------------------------------------------------

a=array()
for(ng in 1:length(generations)){
  a[ng]=simMemoirRandomGuess(generations[ng],mu,alpha,barcodes[1],nRepeats)
}

#this function calculates the area under the curve of the ecdf of the vector. In this case, the cumulative fraction
#of trees (as the plot in Memoir 1.0), that have a certain degree of accuracy (RF.dist)
auc.ecdf<-function(x){
  emp.cdf = ecdf(x)
  return(max(cumsum(emp.cdf(min(x):max(x)))))
}


manual.idx=13
upgam.idx=11
alpha=1/2

distBit=array(0,dim=c(length(barcodes),length(generations)))
distTrit=array(0,dim=c(length(barcodes),length(generations)))

# # # # # # #
pathPlots="/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingMay1st/"
#plot for each generation, RF dist vs BC, overlap blue and black line for trit vs bit
#estimated edit rate 0.4 edits per site per generation
mIdx=3
simType=muVariation[[mIdx]]
pdf(paste(pathPlots,"RFdist_vsBC_empiricMu_compareGenerations.pdf",sep=""))
par(mfrow=c(2,3))
for(ng in 1:length(generations)){
  
  for(bc in 1:length(barcodes)){
    #distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
    distTrit[bc,ng]=apply(simType[['trit']][[bc]][[ng]],2,mean)[13] # 12 -> RF.dist using manualDist + UPGMA
    distBit[bc,ng]=apply(simType[['trit']][[bc]][[ng]],2,mean)[11] #11 -> RF.dist using dist.ml + UPGMA
  }
  
  plot(log(barcodes),distBit[,ng],main=paste("gen = ",toString(generations[ng]),sep=""),
       type="o",ylim=c(0,a[ng]+0.05*a[ng]),ylab="RF dist",xlab="log N scratchpads",
       cex.lab=1.5,cex.axis=1.5,cex.main=2,lwd=1.5)
  lines(log(barcodes),distTrit[,ng],type="o",col="blue",lwd=1.5)
  
  
  abline(h=a[ng],col="red",lwd=1.5)
}
dev.off()

#normalize distances with the random guess and plot heatmaps

distBitNorm=array(0,dim=dim(distBit));for(d in 1:dim(distBit)[2]){distBitNorm[,d]=(a[d]-distBit[,d])/a[d]}

distTritNorm=array(0,dim=dim(distTrit));for(d in 1:dim(distTrit)[2]){distTritNorm[,d]=(a[d]-distTrit[,d])/a[d]}

#rename columns

rownames(distBitNorm)=as.character(barcodes)
colnames(distBitNorm)=as.character(generations)

rownames(distTritNorm)=as.character(barcodes)
colnames(distTritNorm)=as.character(generations)

#plot heatmaps
#normalized distances (by random tree) using the empirical edit rate of 0.4
plotPath= "/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingMay1st/"
pdf(paste(plotPath,"heatmap_Distance_BitVSTrit_muNOTRACE",toString(mus[mIdx]),".pdf",sep=""))
heatmap.2(t(distBitNorm),dendrogram='none', Rowv=FALSE, Colv=FALSE,col=paste("gray",1:99,sep=""),main="Bit",trace="none",xlab="Barcodes",ylab="Generations")
heatmap.2(t(distTritNorm),dendrogram='none', Rowv=FALSE, Colv=FALSE,col=paste("gray",1:99,sep=""),main="Trit",trace="none",xlab="Barcodes",ylab="Generations")

vals=unique(scales::rescale(c(distBitNorm)))
o<-order(vals,decreasing=F)
cols<-scales::col_numeric("Blues",domain=NULL)(vals)
colz<-setNames(data.frame(vals[o],cols[o]),NULL)



p1=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=distBitNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p1,limits=c(0,1))
p2=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=distTritNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p2,limits=c(0,1))
subplot(p1,p2)


dev.off()

#compare the different distance measures: 
distAllCompare = array(0,dim=c(12,length(generations),length(mus)))
par(mfrow=c(2,3))
for(ng in 1:length(generations)){
  for (mIdx in 1:length(mus)){
    simType=muVariation[[mIdx]]
    
    for(d in 1:12){
      #distBitAll[bc,ng,mIdx]=apply(simType[['binary']][[bc]][[ng]],2,mean)[12] #12 is normal UPGMA + as.dist
      distAllCompare[d,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[d]
      #distBitAll[bc,ng, mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[12]
      
    }
  }
  #plot(mus,1-distTritAll[3,ng,]/a[ng],ylim=c(0,1),
  #     type="o",col="blue",ylab="norm dist",xlab="edit rate p/site p/gen",
  #     main=paste("g=",toString(generations[ng]),sep=""),
  #     cex.axis=1.5,cex.lab=1.6,cex.main=2);
  #lines(mus,1-distBitAll[3,ng,]/a[ng],type="o")
  
  #lines(mus,1-distTritAll[5,ng,]/a[ng],col="#3399FF",pch=17,lty=6)
  #lines(mus,1-distBitAll[5,ng,]/a[ng],col="gray",pch=17,lty=6)
  
}



# Apr 30th
# compare bit vs trit using optimal edit rate
pdf(paste(pathPlots,"compareSixmer_forMu_rate.pdf"))
distBitAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))
distTritAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))
par(mfrow=c(2,3))
BCn=2
for(ng in 1:length(generations)){
  for (mIdx in 1:length(mus)){
    simType=muVariation[[mIdx]]
    
    for(bc in 1:length(barcodes)){
      distBitAll[bc,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[11] #12 is normal UPGMA + as.dist
      distTritAll[bc,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[13]
      
      
    }
  }
  plot(mus,1-distTritAll[BCn,ng,]/a[ng],ylim=c(0,1),
       type="o",col="blue",ylab="norm dist",xlab="edit rate p/site p/gen",
       main=paste("g=",toString(generations[ng]),sep=""),
       cex.axis=1.5,cex.lab=1.6,cex.main=2);
  lines(mus,1-distBitAll[BCn,ng,]/a[ng],type="o")
  
  #lines(mus,1-distTritAll[5,ng,]/a[ng],col="#3399FF",pch=17,lty=6)
  #lines(mus,1-distBitAll[5,ng,]/a[ng],col="gray",pch=17,lty=6)
  
}
dev.off()
#take minimun across all edit rates
#optimal rate

optimMatTrit=apply(distTritAll,c(1,2),min)
optimMatBit=apply(distBitAll,c(1,2),min)

pdf(paste(pathPlots,"histograms_tenmer_3gen_muVar.pdf"))
par(mfrow=c(2,2))
ng=1;bc=2
for(m in 1:4){
  medDist = mean(1-muVariation[[m]][['trit']][[bc]][[1]][,13]/a[ng])
  medianDist = median(1-muVariation[[m]][['trit']][[bc]][[1]][,13]/a[ng])
  hist(1-muVariation[[m]][['trit']][[bc]][[ng]][,13]/a[ng],
       main=paste("u=",toString(mus[m]),", <d>=",toString(medDist),", me=",toString(medianDist),sep=""),
       cex.lab=1.5,cex.axis=1.5,cex.main=1.7,
       xlab="norm distance",xlim=c(0,1),ylim=c(0,45))
}
dev.off()



# across generations


pdf(paste(pathPlots,"RFdist_vsBC_optimalRate.pdf",sep=""))
par(mfrow =c(2,3))
barcodeAxis= log(barcodes)
for(ng in 1:length(generations)){
  plot(barcodeAxis,(a[ng]-optimMatTrit[,ng])/a[ng],
       cex.main=2,cex.lab=1.5,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),cex=1.5,main=paste("gen=",toString(generations[ng]),sep=""),
       xlab="log N scratchpads",ylab="Dist norm to random guess",lwd=1.5);
  lines(barcodeAxis,(a[ng]-optimMatBit[,ng])/a[ng],col="black",type="o",lwd=1.5)
  abline(h=a[ng],col="red")
}
dev.off()

pdf(paste(pathPlots,"RFdist_vsBC_optimalRate_vsEmpiricalBIT.pdf",sep=""))
par(mfrow =c(2,3))
barcodeAxis= log(barcodes)
for(ng in 1:length(generations)){
  plot(barcodeAxis,(a[ng]-optimMatBit[,ng])/a[ng],
       cex.main=2,cex.lab=1.5,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),cex=1.5,main=paste("gen=",toString(generations[ng]),sep=""),
       xlab="log N scratchpads",ylab="Dist norm to random guess",lwd=1.5);
  lines(barcodeAxis,distBitNorm[,ng],col="orange",type="o",lwd=1.5)
  abline(h=a[ng],col="red")
}
dev.off()




#heatmaps of RF.dist for BC vs nG
optimMatBitNorm=array(0,dim=dim(optimMatBit));
optimMatTritNorm=array(0,dim=dim(optimMatTrit));
for(o in 1:dim(optimMatBit)[2]){
  optimMatBitNorm[,o]=optimMatBit[,o]/a[o]
  optimMatTritNorm[,o]=optimMatTrit[,o]/a[o]
}

#subplots 
#once the optimal reates were taken into account *by considering the minimum of the distance 
#across all values that were simulated, we can plot the distance in the color axis and consider
#the number of scratchpads and the number of gnerations in the x and y axis. 
#this plot is from c(0,1), 1->perfect reconstruction; 0->random guess

pdf(paste(pathPlots,"heatmaps_blue_BitvsTrit_atOptimalRates.pdf"))
vals=unique(scales::rescale(c(1-optimMatBitNorm)))
o<-order(vals,decreasing=F)
cols<-scales::col_numeric("Blues",domain=NULL)(vals)
colz<-setNames(data.frame(vals[o],cols[o]),NULL)

p1=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=1-optimMatBitNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p1,limits=c(0,1))
p2=plot_ly(y=as.character(log(barcodes)),x=as.character(generations), colorscale = colz,z=1-optimMatTritNorm,type="heatmap",cauto=F,cmin=0,cmax=1) %>% colorbar(p1,limits=c(0,1))

subplot(p1,p2)
dev.off()

x11()
plotList=list()
par(mfrow =c(2,3))
muVar = array(0,dim =c(length(barcodes),length(mus)))
for (ng in 1:length(generations)){
  for (muIdx in 1:length(mus)){
    for (bc in 1:length(barcodes)){
      muVar[bc,muIdx]= apply(muVariation[[muIdx]][['trit']][[bc]][[ng]],2,mean)[12]
    }
  }
  plotList[[ng]]=plot_ly(y=as.character(log(barcodes)),x=as.character(mus),z=muVar,type="heatmap") 
  matplot(log(barcodes),muVar)
}
legend("right",as.character(mus),col=seq_len(length(mus)),cex=0.8,fill=seq_len(length(mus)),title="edit rate")

#matplot(mus,t(apply(muVarMatrix,c(1,2),mean)),type="o",pch=1:4)
#legend("right",distNames,col=seq_len(length(distNames)),cex=0.8,fill=seq_len(length(distNames)),title="Generations")





