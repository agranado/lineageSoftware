
#script uses the probabilistic distances/model to infer the most likely ancestor for sister cells.


ml.ancenstor <-function(bc1,bc2,mu,alpha_,nGen){





}

Pr_edit <- function (nGen,mu,alpha){
  #probability that at nGen generations there is a mutation
  cumSum =0
  if(nGen>0)
    for (nG in 1:nGen)
      cumSum = cumSum + mu * (1-mu)^(nG-1) * alpha
  else
    cumSum =0
    
  return(cumSum)
}

Pr_noedit<-function(nGen,mu){
  (1-mu)^nGen
}

#this function gives the pr that a site is in each of the thre states for a given number of generations.
#works fine as May 8th

Pr_s0<-function(a,mu,alpha,nGen){
  if(a=="u")
    pr = Pr_noedit(nGen,mu)
  else if(a=="r")
    pr =Pr_edit(nGen,mu,alpha)
  else
    pr =Pr_edit(nGen,mu,1-alpha)

  return(pr)
}


manualDistML <- function(barcodeLeaves,mu,alpha,nGen){
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
        #NEW section for W
        if(barcodeArray1[s]=="w" | barcodeArray2[s]=="w"){ #any of them is W


          #if the first one is w
          if(barcodeArray1[s]=="w" & barcodeArray2[s]=="w"){
            sumPr=1;Pr.ind[s] =1
          }else if(barcodeArray1[s]=="w"){
              Pr.ind[s]=Pr[which(alphabet ==barcodeArray2[s])] # Pr(non-w char) * Pr(w)\\

              s2=which(alphabet==barcodeArray2[s])

              #sum over all possible S0 states
              sumPr=0
              for(a in 1:length(alphabet)){
                 sumPr = sumPr+ Tran.pr[a,s2] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
              }


          }else if(barcodeArray2[s]=="w"){
              Pr.ind[s]=Pr[which(alphabet ==barcodeArray1[s])]


              s1=which(alphabet==barcodeArray1[s])

              sumPr=0
              for(a in 1:length(alphabet)){
                 sumPr = sumPr+ Tran.pr[a,s1] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
              }

          }

        }else{

          #NORMAL

            Pr.ind[s] = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])]
            #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
            # Sum{i=1}{3} Pr(s1|S0=i)Pr(s2|s0=i)Ps(S0=i)
            s1=which(alphabet==barcodeArray1[s])
            s2=which(alphabet==barcodeArray2[s])

            #sum over all possible S0 states
            sumPr=0
            for(a in 1:length(alphabet)){
               sumPr = sumPr+ Tran.pr[a,s1] * Tran.pr[a,s2] * Pr_s0_array[a]
                #   All possible ways in which the sisters s1, s2 could be generated from a previous state.
                #   The sum is all over possible states of S0, weighted by their proabilities.
                #   Pr_{g-1}(S_0) x Pr_{g}(S_0->S_1) x P_{g}(S_0->S_2)  # Pr_s0(alphabet[a],mu,alpha,nGen-1)
            }
        }
        Pr.s1_s2[s] =sumPr

      } #end barcode loop
        ratio.product=sum(log(Pr.s1_s2))/sum(log(Pr.ind))

      distMat[i,j]= distSum
      ratioMat[i,j]=1/ratio.sum *distSum
      productMat[i,j] = ratio.product

    }
  }
  return(productMat)
}


manualDistML__ <- function(barcodeLeaves,mu,alpha,nGen){
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
      } #end barcode loop
        ratio.product=sum(log(Pr.s1_s2))/sum(log(Pr.ind))

      distMat[i,j]= distSum
      ratioMat[i,j]=1/ratio.sum *distSum
      productMat[i,j] = ratio.product

    }
  }
  return(productMat)
}
