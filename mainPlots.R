#May 9th
#main plots that might be used in the paper
#this script should load the object muVariation or it should be executed afeter running bitVStrit.R

barcodes = c(2,4,6,8,10,12,15,17,20,30,40,50,60,70,80,100,200)
#barcodes= c(2,10,50,100)
generations=c(2,3,4,5,6,7,8,9,10)
#generations=c(8)
mus = c(0.999,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.001)
#mus = c(0.99,0.6,0.5,0.4,0.1,0.01)
#mus=c(0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,0.001)
#barcodes = c(6,7)
#generations = c(3,4,5)
nRepeats=72


# plotting ----------------------------------------------------------------
source("simMemoirRandomGuess.R")
 # mu=0.4;alpha=1/2;
 # b=array(0,dim=c(length(generations),nRepeats))
 # a = array()
 # for(ng in 1:length(generations)){
 #   b[ng,]=simMemoirRandomGuess(generations[ng],mu,alpha,barcodes[1],nRepeats)
 #   a[ng] = median(b[ng,])
 # }
a=c(2    ,10,   26,   58,  122,  250,  506, 1018, 2042)


#this function calculates the area under the curve of the ecdf of the vector. In this case, the cumulative fraction
#of trees (as the plot in Memoir 1.0), that have a certain degree of accuracy (RF.dist)
auc.ecdf<-function(x){
  emp.cdf = ecdf(x)
  return(max(cumsum(emp.cdf(min(x):max(x)))))
}


#compare the measures

distMat=array(0,dim=c(length(barcodes),length(generations),3))
mIdx=7
simType=muVariation[[mIdx]]
par(mfrow=c(2,5))
for(ng in 1:length(generations)){

  for(bc in 1:length(barcodes)){
    #distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
    distMat[bc,ng,]=apply(simType[['trit']][[bc]][[ng]],2,mean) # 12 -> RF.dist using manualDist + UPGMA

  }

  matplot(log(barcodes),distMat[,ng,],main=paste("gen = ",toString(generations[ng]),sep=""),
       type="o",ylim=c(0,a[ng]+0.05*a[ng]),ylab="RF dist",xlab="log N scratchpads",
       cex.lab=1.5,cex.axis=1.5,cex.main=2,lwd=1.5)

  abline(h=a[ng],col="red",lwd=1.5)
}

library(RColorBrewer)
darkcols <- brewer.pal(3, "Dark2")
#compare trit performance using 3 different distance measures
#using the Mean/median to compare
distMat=array(0,dim=c(length(barcodes),length(generations),3))
simType=muVariation[[mIdx]]
par(mfrow=c(2,5))
for(ng in 1:length(generations)){

  for(bc in 1:length(barcodes)){
    #distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
    distMat[bc,ng,]=apply(simType[['trit']][[bc]][[ng]],2,mean) # 12 -> RF.dist using manualDist + UPGMA

  }

  matplot(log10(barcodes),distMat[,ng,],main=paste("gen = ",toString(generations[ng]),sep=""),
          type="o",ylim=c(0,a[ng]+0.05*a[ng]),ylab="RF dist",xlab="log N scratchpads",
          cex.lab=1.5,cex.axis=1.5,cex.main=2,lwd=1.5,col = darkcols)


  legend("right",as.character(c(1,2,3)),col=darkcols,cex=0.8,fill=darkcols,title="edit rate")
  abline(h=a[ng],col="red",lwd=1.5)

}
#
#
#
#Proportion of perfect trees
bcIdx=5;gIdx=1
distMat=array(0,dim=c(length(barcodes),length(generations),3))
#mIdx=7
simType=muVariation[[mIdx]]
par(mfrow=c(2,5))
for(ng in 1:length(generations)){

  for(bc in 1:length(barcodes)){
    #distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
    distMat[bc,ng,]=apply(simType[['trit']][[bc]][[ng]]==0,2,sum)/dim(simType[['trit']][[bcIdx]][[1]]==0)[gIdx] # 12 -> RF.dist using manualDist + UPGMA

  }

  matplot(log10(barcodes),distMat[,ng,],main=paste("gen = ",toString(generations[ng]),sep=""),
          type="o",ylim=c(0,1),ylab="Prop perfect trees",xlab="log N scratchpads",
          cex.lab=1.5,cex.axis=1.5,cex.main=2,lwd=1.5,col = darkcols)


  legend("right",as.character(c(1,2,3)),col=darkcols,cex=0.8,fill=darkcols,title="edit rate")
  abline(h=a[ng],col="red",lwd=1.5)

}


manual.idx=3
upgam.idx=1
alpha=1/2
pdf(paste(pathPlots,"RFdist_vsBC_empiricMu_compareGenerations.pdf",sep=""))


distBit=array(0,dim=c(length(barcodes),length(generations)))
distTrit=array(0,dim=c(length(barcodes),length(generations)))

# # # # # # #
pathPlots="/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingMay15/"
#plot for each generation, RF dist vs BC, overlap blue and black line for trit vs bit
#estimated edit rate 0.4 edits per site per generation

#mIdx=5
simType=muVariation[[mIdx]]
par(mfrow=c(2,5))
for(ng in 1:length(generations)){

  for(bc in 1:length(barcodes)){
    #distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[11]
    distTrit[bc,ng]=apply(simType[['trit']][[bc]][[ng]],2,mean)[manual.idx] # 12 -> RF.dist using manualDist + UPGMA
    distBit[bc,ng]=apply(simType[['binary']][[bc]][[ng]],2,mean)[upgam.idx] #11 -> RF.dist using dist.ml + UPGMA
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
plotPath= "/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingMay15/"

library(plotly)
library(ggplot2)
library(gplots)
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
distAllCompare = array(0,dim=c(3,length(generations),length(mus)))
par(mfrow=c(2,3))
for(ng in 1:length(generations)){
  for (mIdx in 1:length(mus)){
    simType=muVariation[[mIdx]]

    for(d in 1:3){
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



# May 24th
# compare bit vs trit using optimal edit rate
pdf(paste(pathPlots,"compare10mer_forMu_rate_may24_MEAN.pdf"))
#paper
BCn=5
distBitAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))
distTritAll=array(0,dim=c(length(barcodes),length(generations),length(mus)))
par(mfrow=c(2,5))
for(ng in 1:length(generations)){
  for (mIdx in 1:length(mus)){
    simType=muVariation[[mIdx]]

    for(bc in 1:length(barcodes)){
      distBitAll[bc,ng,mIdx]=apply(simType[['binary']][[bc]][[ng]],2,mean)[1] #12 is normal UPGMA + as.dist
      distTritAll[bc,ng,mIdx]=apply(simType[['trit']][[bc]][[ng]],2,mean)[3]


    }
  }
  plot(mus,1-distTritAll[BCn,ng,]/a[ng],ylim=c(0,1),
       type="o",col="blue",ylab="Reconstructability",xlab="edit rate [p/site p/gen]",
       main=paste("g=",toString(generations[ng]),sep=""),
       cex.axis=1.5,cex.lab=1.8,cex.main=2,lwd=2,cex.axis=1.5);
  lines(mus,1-distBitAll[BCn,ng,]/a[ng],type="o",lwd=2)

  #lines(mus,1-distTritAll[5,ng,]/a[ng],col="#3399FF",pch=17,lty=6)
  #lines(mus,1-distBitAll[5,ng,]/a[ng],col="gray",pch=17,lty=6)

}
dev.off()
#take minimun across all edit rates
#optimal rate





pdf(paste(pathPlots,"ecdf_tenmer_3gen_muVar_BIT_May24th.pdf"))


optimMatTrit=apply(distTritAll,c(1,2),min)
optimMatBit=apply(distBitAll,c(1,2),min)

par(mfrow=c(2,2))
ng=1;bc=5
for(m in 1:11){
  medDist = mean(1-muVariation[[m]][['binary']][[bc]][[1]][,1]/a[ng])
  medianDist = median(1-muVariation[[m]][['binary']][[bc]][[1]][,1]/a[ng])
  plot(ecdf(1-muVariation[[m]][['binary']][[bc]][[ng]][,1]/a[ng]),
       main=paste("u=",toString(mus[m]),", <d>=",toString(medDist),", me=",toString(medianDist),sep=""),
       cex.lab=1.5,cex.axis=1.5,cex.main=1.7,
      xlab="norm distance",xlim=c(0,1),ylim=c(0,1))
}
dev.off()



# across generations

#paper (optimal rate bit vs trit)
pdf(paste(pathPlots,"RFdist_vsBC_optimalRate_May24th_MEAN.pdf",sep=""))
par(mfrow =c(2,5))
barcodeAxis= log10(barcodes)

for(ng in 1:length(generations)){
  plot(barcodeAxis,(a[ng]-optimMatTrit[,ng])/a[ng],
       cex.main=2,cex.lab=1.8,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),main=paste("gen=",toString(generations[ng]),sep="",log="x",xaxt="n"),
       xlab="N recording units",ylab="Dist norm to random guess",lwd=2);
  lines(barcodeAxis,(a[ng]-optimMatBit[,ng])/a[ng],col="black",type="o",lwd=2)
  abline(h=a[ng],col="red")
}

at.y <- outer(1:9, 10^(0:4))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1)


dev.off()

#paper (empirical vs optimal rate)
pdf(paste(pathPlots,"RFdist_vsBC_optimalRate_vsEmpiricalTRIT_MEAN.pdf",sep=""))
par(mfrow =c(2,5))
barcodeAxis= log10(barcodes)
for(ng in 1:length(generations)){ #generations to plot length(generations) : to plot all
  plot(barcodeAxis,(a[ng]-optimMatTrit[,ng])/a[ng],
       cex.main=2,cex.lab=1.8,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),cex=1.5,main=paste("gen=",toString(generations[ng]),sep=""),
       xlab="log N units",ylab="Reconstructability",lwd=2);
  lines(barcodeAxis,distTritNorm[,ng],col="orange",type="o",lwd=2)
  abline(h=a[ng],col="red")
}
dev.off()


pdf(paste(pathPlots,"RFdist_vsBC_optimalRate_vsEmpiricalTRIT.pdf",sep=""))
par(mfrow =c(2,2))
barcodeAxis= log(barcodes)
for(ng in 1:length(1:2)){ #generations to plot length(generations) : to plot all
  plot(barcodeAxis,(a[ng]-optimMatTrit[,ng])/a[ng],
       cex.main=2,cex.lab=1.5,cex.axis=1.5,type="o",
       col="blue",ylim=c(0,1),cex=1.5,main=paste("gen=",toString(generations[ng]),sep=""),
       xlab="log N units",ylab="Reconstructability",lwd=2);
  lines(barcodeAxis,distTritNorm[,ng],col="orange",type="o",lwd=2)
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

#paper plot
# HHMI plots / retreat plots
#how many generations you can do with N barcodes
par(mfrow=c(2,2))
optimMatNorm= array(0, dim=dim(optimMatTrit) ); for(i in 1:dim(optimMatTrit)[2]){ optimMatNorm[,i] = (a[i]-optimMatTrit[,i])/a[i] }
optimMatBitNorm= array(0, dim=dim(optimMatTrit) ); for(i in 1:dim(optimMatBit)[2]){ optimMatBitNorm[,i] = (a[i]-optimMatBit[,i])/a[i] }
t=c(0.80,0.85,0.90,0.95)
#ylims=c(100,100,100,200
ylims=c(60,60,60,100)
legend.y.pos = c(90,90,90,180)
for(tt in 1:length(t)){


    threshold=t[tt]
    #bc.per.gen.bit = array();for(i in 1:dim(optimMatBitNorm)[2]){ bc.per.gen.bit[i] = min( which(optimMatBitNorm[,i]>threshold) ) }
    bc.per.gen.bit = array();for(i in 1:dim(optimMatNorm)[2]){ bc.per.gen.bit[i] = min( which(optimMatNorm[,i]>t[3]) ) }
    bc.per.gen.bit2 = array();for(i in 1:dim(optimMatNorm)[2]){ bc.per.gen.bit2[i] = min( which(optimMatNorm[,i]>t[1]) ) }
    bc.per.gen.trit = array();for(i in 1:dim(optimMatNorm)[2]){ bc.per.gen.trit[i] = min( which(optimMatNorm[,i]>threshold) ) }

    plot(generations,barcodes[bc.per.gen.trit],type="l",lwd=2.5,
         xlim=c(2,10),ylim=c(2,ylims[tt]),cex.lab=1.5,cex.axis=1.2,ylab="Array length",xlab="Generations",
         main=paste(toString(t[tt])," % reconstructability",sep=""))
    grid (10,10, lty = 6, col = "cornsilk2")
    lines(generations,barcodes[bc.per.gen.bit],type="l",col="black",lwd=2)
    lines(generations,barcodes[bc.per.gen.trit],type="l",col="blue",lwd=2)
    lines(generations,barcodes[bc.per.gen.bit2],type="l",col="gray",lwd=2)
    #legend(3,legend.y.pos[tt],c("Trit","Bit"),fill=c("blue","black"))
}
