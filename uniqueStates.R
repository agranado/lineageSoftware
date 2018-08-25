#compare bit vs trit
#Apr 23th
#This script runs variations of nGen and barcodeLenght for bit and trit memoir simulations.
#The idea is to compare how well the trit is doing vs the bit.
#runs in parallel.

rm(list=ls())
source("simMemoirStrDist3.R")
source("simulation2.R")
source("MLfunctions.R")
library(doParallel)
os=system("cat ../os.txt",intern = T)
if(os=="mac"){
  registerDoParallel(cores=8)
}else if(os=="linux"){
  registerDoParallel(cores=36)
}

#SET PARAMETERS
barcodes = c(2,4,6,8,10)#,15,20,25,50,100)
barcodes= c(2,4,6,10,15,50,100)
generations=c(3,4,5,6,7,8,9)
generations=c(3,4,5,6,7,8,9)
mus = c(0.9999,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0001)
#mus = c(0.7)#,0.6)#,0.5,0.4,0.2,0.1)
mus =c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.001)
#barcodes = c(6,7)
#generations = c(7,8)
nRepeats=72*2


types=c('binary','trit')
#types=c('trit')

muVariation=list()
for(m in 1:length(mus)){
  simType=list()
  mu=mus[m]
  for(st in 1:length(types)){
    simulationType=types[st]
    barcodeData=list()
    for(bc in 1:length(barcodes)){
      barcodeLength=barcodes[bc]
      genData=list()
      for(ng in 1:length(generations)){
        nGen=generations[ng]
        
        results= foreach(i=1:nRepeats) %dopar% simMemoirBarcodes(nGen=nGen,alpha=1/2,mu=mu,simulationType=simulationType,barcodeLength = barcodeLength)
        results.mat=do.call(rbind,results)
        genData[[ng]]=results.mat
        
        print(paste("sim: g=",toString(nGen)," ",simulationType," mu",toString(mu)," BC=",toString(barcodeLength),sep=""))
      }
      barcodeData[[bc]]=genData
    }
    simType[[simulationType]]=barcodeData
  }
  muVariation[[m]]=simType
}

x11();
#plot average number of unique states for each parameter set: nGen, mu, BC=10
par(mfrow=c(2,3))
for(mIdx in 6:8){
  #mIdx=5
  bc=4
  ng=1
  
  nRep=dim(muVariation[[4]][['binary']][[3]][[1]])[1]
  
  bit.states=array()
  trit.states=array()
  bit.sdev = array()
  trit.sdev = array()
  for(ng in 1:length(generations)){
    bit.unique=array(); 
    trit.unique=array(); 
    for( i in 1:nRep){
      bit.unique[i]=length(unique(muVariation[[mIdx]][['binary']][[bc]][[ng]][i,]))
      trit.unique[i]=length(unique(muVariation[[mIdx]][['trit']][[bc]][[ng]][i,]))
    }
    bit.states[ng]=mean(bit.unique)
    trit.states[ng]=mean(trit.unique)
    
    trit.sdev[ng]=sd(trit.unique)
    bit.sdev[ng]=sd(bit.unique)
  }
  



  plot(generations,bit.states,ylim=c(1,250),type="o",xlab="time (cell divisions)",ylab="Number of unique states",
       lwd=2,cex.lab=1.7,cex.ax=1.7)
  arrows(generations, bit.states-bit.sdev, generations, bit.states+bit.sdev, length=0.05, angle=90, code=3)
  
  lines(generations,trit.states,type="o",col="blue",lwd=2)
  arrows(generations, trit.states-trit.sdev, generations, trit.states+trit.sdev, length=0.05, angle=90, code=3)
}

#plot separate panels, all lines for bit and all lines for trit
x11();
#plot average number of unique states for each parameter set: nGen, mu, BC=10
par(mfrow=c(2,2))
for(simulationType in c('binary','trit')){
  if(simulationType=='binary'){col.plot='black'}else{col.plot='blue'}
  for(mIdx in 1:length(mus)){
    #mIdx=5
    bc=4
    ng=1
    
    nRep=dim(muVariation[[4]][['binary']][[3]][[1]])[1]
    
    bit.states=array()
    trit.states=array()
    bit.sdev = array()
    trit.sdev = array()
    for(ng in 1:length(generations)){
      bit.unique=array(); 
      trit.unique=array(); 
      for( i in 1:nRep){
        bit.unique[i]=length(unique(muVariation[[mIdx]][[simulationType]][[bc]][[ng]][i,]))
      }
      bit.states[ng]=mean(bit.unique)
      bit.sdev[ng]=sd(bit.unique)
    }
    if(mIdx==1){
      plot(generations,bit.states,ylim=c(1,2^9),type="o",xlab="time (cell divisions)",ylab="Number of unique states",lwd=2,cex.lab=1.5,cex.ax=1.5)
      arrows(generations, bit.states-bit.sdev, generations, bit.states+bit.sdev, length=0.05, angle=90, code=3)
    }else{
      
      lines(generations,bit.states,type="o",col=col.plot,lwd=2)
      arrows(generations, bit.states-bit.sdev, generations, bit.states+bit.sdev, length=0.05, angle=90, code=3)
    }
  }
}



#unique sampling from distribution of all states. 
pathPlots="/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/fig1/"

bc.length=10
n.cells=10:1000
pr.bit=array()
pr.trit=array()


pdf(paste(pathPlots,"Pr_NO_redundantSet.pdf",sep=""))
pr.bit[1] = 1
pr.trit[1] = 1
par(mfrow=c(2,2))
for (i in 2:length(n.cells)){

    pr.bit[i] = (2^bc.length-(i-1))/2^bc.length   * pr.bit[i-1]
    pr.trit[i]= (3^bc.length-(i-1))/3^bc.length * pr.trit[i-1]
 
}

plot(pr.trit,type="l",col="blue",lwd=2,ylab="Probability of no redundancy",xlab="N cells",
     xlim=c(0,250),cex.lab=1.7,cex.axis=1.7);
lines(pr.bit,type="l",col="black",lwd=2)

dev.off()



nRepeats=100
n.unique=array()
n.unique2=array()
max.N=5
n.cells=10000

clone.distribution.bit = array(0,dim=c(n.cells,nRepeats,max.N))
clone.distribution.trit = array(0,dim=c(n.cells,nRepeats,max.N))

for(i in 1:n.cells){

  for(j in 1:nRepeats){
    #generate the samples
    this.sample.bit = sample(1:2^10,i,replace=T);
    this.sample.trit = sample(1:3^10,i,replace=T)
    #count the number of sinlges, pairs, triaads, etc
    this.sample.bit.table = table(this.sample.bit)
    this.sample.trit.table = table(this.sample.trit)
                                   
    for(z in 1:max.N){
      clone.distribution.bit[i,j,z] = (sum(this.sample.bit.table==z)/i) * (1/z)^(z-1)
      clone.distribution.trit[i,j,z] = (sum(this.sample.trit.table==z)/i) * (1/z)^(z-1)

    }
  }
}

#plot the sum of the values (as the total probability of clonality)
pdf(paste(pathPlots,"variationBit_Trit.pdf",sep=""))


par(mfrow=c(2,2))

clone.conf.bit = array()
clone.conf.trit = array()
for(i in 1:n.cells){
  clone.conf.bit[i]=sum(apply(clone.distribution.bit[i,,],2,mean))
  clone.conf.trit[i]=sum(apply(clone.distribution.trit[i,,],2,mean))
  
}

plot(clone.conf.bit,lwd=2,type="l",cex.lab=1.7,cex.axis=1.7,ylab="clonal confidence",xlab="N cells",ylim=c(0.5,1),
     main="Total probability of clonality",log="x")
lines(clone.conf.trit,lwd=2,col="blue")
abline(h=0.95,lwd=3,lty=2,col="red")


#plot the fraction of unique states

clone.conf.bit = array()
clone.conf.trit = array()
clone.conf.bit.sd = array()
clone.conf.trit.sd = array()


for(i in 1:n.cells){
  clone.conf.bit[i]=mean(clone.distribution.bit[i,,1])
  clone.conf.trit[i]=mean(clone.distribution.trit[i,,1])
  
  clone.conf.bit.sd[i]=sd(clone.distribution.bit[i,,1])
  clone.conf.trit.sd[i]=sd(clone.distribution.trit[i,,1])
  
  
}

plot(clone.conf.bit,lwd=2,type="l",cex.lab=1.7,cex.axis=1.7,
     log="x",ylab="clonal confidence",xlab="N cells",ylim=c(0.5,1),
     main="Fraction of distinct states")
lines(clone.conf.trit,lwd=2,col="blue")
abline(h=0.95,lwd=3,lty=2,col="red")

dev.off()
#direct calculation 

#####################

c1<-data.frame(1:10000,clone.conf.bit)
 c1$type="bit"
 c2<-data.frame(1:10000,clone.conf.trit)
 c2$type="trit"

 names(c1)<-c("N","prob","type")
 names(c2)<-c("N","prob","type")
 
 A=rbind(c1,c2)
 
 sd.test.max=c(clone.conf.bit+clone.conf.bit.sd,clone.conf.trit+clone.conf.trit.sd)
 sd.test.min=c(clone.conf.bit-clone.conf.bit.sd,clone.conf.trit-clone.conf.trit.sd)
 
xlabel<- "N cells" #This is just a way of getting the "500" as a subscript (See graph above)
ylabel <- "Clonal probability"


p<-ggplot(data=A, aes(x=A$N, y=prob, colour=type,linetype=type)) + geom_point() + geom_line();
p<-p+geom_ribbon(aes(ymin=sd.test.min, ymax=sd.test.max), alpha=0.1); 
p<-p+scale_x_log10()
p


#ggplot2 allows you to add pieces to a plot as you go along. 
#I've rewritten the original code to make this more obvious
#First, we enter the data. Note that this odes not plot anything yet.
#fill=type says to set the fill separately for each A$type. 
#Likewise, linteype=type says to use a different line type for each A$type. Same for color=type, if we added that.
p <- ggplot(data=A, aes(x=A$N, y=prob, ymin=sd.test.min, ymax=sd.test.min, fill=type, linetype=type))

#To actually plot this, we'd type "plot(p)" or just "p". Right now, that will give an error since we haven't said what to do with the data yet!
#In the following, geom_line() says to draw a line graph using these data. 
p <- p + geom_line()
#If you want to fill in the ribbons, add
p <- p + geom_ribbon(alpha=.5) #obviously, you can set the alpha as you like
# This next line puts us on a log-log scale
p <- p + scale_x_log10() 
p <- p + xlab(xlabel) + ylab(ylabel)
#now plot
p







N=2^10;k=0;m=100;
Pr_0=array()
for(i in 1:10000){
  m=i
  Pr_0[i]=(1-choose(m,k)* (N-1)^(m-k)/N^m) *N
}

