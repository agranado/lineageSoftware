# 6-Sep-2018


#take a distance matrix and calculate the minimum distance for each cell,
#we are going to assume that this is the sister and we will infer the most likely commmon ancenstor

#the object matdist_ is upper diagonal


  # a= matdist_ + t(matdist_)

  #a already has names

  #> a[1:2,1:2]
#                   19_ggagaaagta 30_gtgaaaaaaa
#19_ggagaaagta      0.000000      1.159583
#30_gtgaaaaaaa      1.159583      0.000000

# diagonal of a is 0 , so we make it NaN to fin dthe min values.

  a[a==0]=NaN
  #make it full matrix, symmetric
  aa= a + t(a)

  #for each cell, get the cell with the minimum distance
  sis1=colnames(aa)[apply(aa,1,which.min)]
  sis2=names(apply(aa,1,which.min))


  #extract the sequence part of the names (they do have a number at the beggining)
  m<-regexpr("[gta]+",sis1);
  sis1=regmatches(sis1,m)

  m2<-regexpr("[gta]+", sis2)
  sis2=regmatches(sis2,m2)
 #convert to MEMOIR notation
  a=gsub("g","u",sis1); a=gsub("t","x",a); a=gsub("a","r",a)
  b=gsub("g","u",sis2); b=gsub("t","x",b);b= gsub("a","r",b)

  sis1 = a
  sis2 = b
  for (i in 1:length(sis1)){ #for each pair of cells
    barcodeArray1 =strsplit(a[i],"")[[1]]
    barcodeArray2 =strsplit(a[i],"")[[1]]


    for (j in 1:length(barcodeArray1))


    for(a in 1:length(alphabet)){
       sumPr = sumPr+ Tran.pr[a,s1] * Tran.pr[a,s2] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
    }


  }



  manualDistML <- function(a,b,mu,alpha,nGen){
    alphabet = array()
    alphabet[1]= "u"
    alphabet[2]= "r"
    alphabet[3]= "x"
    #mu & alpha indicate an explicit probabilistic model
    #Probability of no mutation for nGen cell divisions (mu is rate of edit per cell division)

    #transition probabilities

    #this is of the form i,j   Pr(i->j | i)
    Tran.pr = t(matrix(c(1-mu, mu*alpha,mu*(1-alpha),0,1,0,0,0,1),nrow=3,ncol=3))


    #NULL MODEL: probability of observing sites as independent events:
    Pr = array()
    Pr[1] = (1-mu)^nGen

    Pr[2] = Pr_edit(nGen,mu,alpha)

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

    Pr_s0_array = array()
    for (i in 1:length(alphabet)){
      Pr_s0_array[i]= Pr_s0(alphabet[i], mu, alpha, nGen-1)
    }
    #go through all the elements in the barcode array
    for (i in 1:(nBarcodes-1)){
      barcodeArray1 =strsplit(barcodeLeaves[i],"")[[1]]
      for (j in (i+1):nBarcodes){

        barcodeArray2 =strsplit(barcodeLeaves[j],"")[[1]]
        distSum = 1
        ratio.sum =0
        ratio.product=1
        #for each pairwise comparison
        Pr.s1_s2=array()
        Pr.ind=array()
        for (s in 1:barcodeLength){
          Pr.ind[s] = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])]
          #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
          # Sum{i=1}{3} Pr(s1|S0=i)Pr(s2|s0=i)Ps(S0=i)
          s1=which(alphabet==barcodeArray1[s])
          s2=which(alphabet==barcodeArray2[s])

          #sum over all possible S0 states
          sumPr=0
          for(a in 1:length(alphabet)){
             sumPr = sumPr+ Tran.pr[a,s1] * Tran.pr[a,s2] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
          }
          Pr.s1_s2[s] =sumPr
        }
          ratio.product=sum(log(Pr.s1_s2))/sum(log(Pr.ind))

        distMat[i,j]= distSum
        ratioMat[i,j]=1/ratio.sum *distSum
        productMat[i,j] = ratio.product

      }
    }
    return(productMat)
  }
