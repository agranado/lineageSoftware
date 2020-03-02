

pathPlots="/Users/alejandrog/MEGA/Caltech/trees/simulation/plots/meetingJul11/"

#Meeting July 11
#Shows Reconstructability as a function of the number of generations. 
#Also, boxplots for time series. 
#Show that 40 units is the minumum to reconstruct 8 generaitons with >0.8 reconstructability

random.dist=c(10  ,26  ,58 ,122 ,250 ,506)
random.dist = c(2    ,9,   26,   58,  122,  250  ,506 ,1018 ,2042)


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
nrep = nRepeats

types=c('binary','trit')
#Plots supplementary and subgroup meeting. 
#get optimal reconstructability as distribution for all simulations
#previously I was showing only the mean

df.dist = list() #list of data frames
for(g in 1:length(generations)){       
  
      optim.ng=g
  
      all.dist = array()
      all.categ = array()
      all.axis = array()
      ml.idx = 3# which index corresponds to the ML distance
      for(o in 1:length(barcodes)){
        optim.bc=o #10-mer
        
          
        #Find the optimal edit rate for this combination of generations and barcode, 
        mIdx=which(distTritAll[optim.bc,optim.ng,]==min(distTritAll[optim.bc,optim.ng,]),arr.ind=T)
        mIdx=mIdx[1]
        
        simType=muVariation[[mIdx]]
        
        optim.dist = simType[['trit']][[optim.bc]][[optim.ng]]
        optim.dist=(random.dist[g]-optim.dist)/random.dist[g]
        all.dist = c(all.dist,optim.dist[,ml.idx])
        all.categ = c(all.categ,rep(toString(barcodes[o]),nrep)  )
        all.axis = c(all.axis,rep(barcodes[o],nrep))
      }
      #df.dist=data.frame(all.dist[2:length(all.dist)],all.categ[2:length(all.categ)],all.axis[2:length(all.axis)])
      
      df.dist[[g]]=data.frame(all.dist,all.categ,all.axis)
}

library(HH)

#all blue boxplots: 
  #####
  #####
boxplot_series<-function(x,y,main.title){
  all.axis=x
  all.dist=y
  df.dist =data.frame(all.axis,all.dist)
  #convert data to factor and log scale
  df.dist$all.axis<- factor(log10(as.numeric(df.dist$all.axis)))
  position(df.dist$all.axis)<-as.numeric(levels(df.dist$all.axis))
  #main call to boxplot function 
  bwplot(all.dist ~ all.axis,data =df.dist,col="blue", scales=list(cex=1.2,x=list(limits=c(0, 2.1), 
      at=position(df.dist$all.axis),labels=round(10^position(df.dist$all.axis))   )),  
      box.width=0.08,panel=panel.bwplot.intermediate.hh,colr=1,do.out=F,
      par.settings=list(box.rectangle=list(col="salmon",fill="blue",
                        alpha=0.3),box.umbrella=list(col="salmon",alpha=0.4)),
      main=main.title,ylab=list("Reconstructability", cex=1.5), 
      xlab=list("Recording units", cex=1.5))
}

  ###
  ###
  ###

#save pdfs

for(g in 1:length(generations)){
  pdf(paste(pathPlots,paste("boxplots_",toString(generations[g]),"G.pdf",sep="" ) ,sep=""),width=5,height = 5)

  boxplot_series(df.dist[[g]]$all.axis,df.dist[[g]]$all.dist,paste(toString(generations[g])," generations",sep=""))
  Sys.sleep(3)
  dev.off()
  

}






#plot the normalized reconstructability 
######
###### Function
plot_reconstruction_matrix<-function(df.dist,log="x",legend.x=40,lwd=2,legend.y=0.4,generations=c(2  ,3,  4,  5  ,6  ,7  ,8)){
   
  all.means=matrix(0,length(generations),length(unique(df.dist[[1]]$all.axis))-1)
  
    cols=brewer.pal(n=9,name="Blues")
    for(g in 1:length(generations)){
      data=df.dist[[g]]
      a=aggregate(data$all.dist , by=list(data$all.axis, data$all.categ) , mean)
      if(g==1){
      plot(sort(a[,1],index.return=T)$x,a[sort(a[,1],index.return=T)$ix,3], type="l" , 
           col=cols[g+2] , lwd=lwd,ylim=c(0,1),log=log,cex.lab=1.5,cex.axis=1.2,
           ylab="Reconstructability",xlab="Recording units")
      }else{ lines(sort(a[,1],index.return=T)$x,a[sort(a[,1],index.return=T)$ix,3],lwd=lwd ,type="l" , col=cols[g])}
      
      
      all.means[g,]=a[sort(a[,1],index.return=T)$ix,3]
    }
    legend(legend.x,legend.y,c("2","3","4","5","6","7","8","9","10"),fill=cols[3:length(cols)],title="Generations",pt.cex=2)
  
    return(all.means)
}


pdf(paste(pathPlots,paste("reconstructionLines_10G.pdf",sep="" ) ,sep=""))

par(mfrow=c(2,2))
legend.y=0.75
plot_reconstruction_matrix(df.dist,"",legend.x=70,legend.y=legend.y,generations = generations)
plot_reconstruction_matrix(df.dist,"x",legend.x=30,legend.y=legend.y,generations=generations)

dev.off()
#########
#########

### plot each recording array as a line
### and see how reconstructability drops as the number of G increases
#add space to the right to make the legend

pdf(paste(pathPlots,paste("compareBarcodes_AllG.pdf",sep="" ) ,sep=""))


#par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
par(mfrow=c(2,2))
matplot(generations,all.means[,c(1,3,5,6,8,9,10)],type="l",
        col=brewer.pal(n=9,name="Purples")[2:9],lwd=2,lty=1,
        cex.axis=1.2,cex.lab=1.5,ylab="Reconstructability",xlab="Generations")

legend("topright",inset=c(-0.8,0),as.character(barcodes[c(1,3,5,6,8,9,10)]),
       fill=brewer.pal(n=9,name="Purples")[2:9],title="Units",pt.cex=2)

dev.off()

# tutorial plotting multiple lines using ggplot (maybe next time)

library("reshape2")
library("ggplot2")

test_data <-
  data.frame(
    var0 = 100 + c(0, cumsum(runif(49, -20, 20))),
    var1 = 150 + c(0, cumsum(runif(49, -10, 10))),
    date = seq(as.Date("2002-01-01"), by="1 month", length.out=100)
  )

test_data_long <- melt(test_data, id="date")  # convert to long format

ggplot(data=test_data_long,
       aes(x=date, y=value, colour=variable)) +
  geom_line()


# plot all boxplots in the same trend. 


#old way to make boxplot 


#convert the axis to log scale

