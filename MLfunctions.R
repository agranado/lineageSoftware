
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


# function for cousing distance

cousinPr<-function(s0, s1, mu,alpha, nGen){

    if(s0 =="u"){
      p = Pr_s0(s1, mu, alpha,2)
    }else if(s0==s1){
      p = 1
    }else{
      p = 0
    }

    return(p)
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
  #if(alpha==0){ alpha = 1e-6} else if (alpha==1){alpha = 0.999999}
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



all.probs<-function(mu,alpha,nGen,alphabet,  cousin = F){
  #this is of the form i,j   Pr(i->j | i)
  Tran.pr = t(matrix(c(1-mu, mu*alpha,mu*(1-alpha),0,1,0,0,0,1),nrow=3,ncol=3))
  #may29 NOTE: here we need a list of Tran.pr matrices, one per site

  # if we calculate cousin distance, the parent is g-2 instead of g-1
  if(cousin) parent_time = 2 else parent_time = 1;

  #NULL MODEL: probability of observing sites as independent events:
  Pr = array()
  Pr[1] = (1-mu)^nGen

  Pr[2] = Pr_edit(nGen,mu,alpha)

  Pr[3] = Pr_edit(nGen,mu,1-alpha)

  PrMatrix  = array(0,dim=c(length(Pr),length(Pr)))
  #calcualte probabilistic model: P(x,y)
  #this just means the pr that those two sites have those characters by random. This is the expected pr for each site
  #it assummes independence but does not tell you how likely they are to come from a common ancestor
  for (p1 in 1:length(alphabet)){
    for (p2 in 1:length(alphabet)){
      PrMatrix[p1,p2] = Pr[p1] * Pr[p2]
    }
  }

  #Probability for a cel in the previous generation to have a given state
  Pr_s0_array = array()
  for (i in 1:length(alphabet)){
    Pr_s0_array[i]= Pr_s0(alphabet[i], mu, alpha, nGen-parent_time)
  }

  return(list(PrMatrix,Pr_s0_array,Tran.pr))
}

#June 19th: account for large deletions
#will modify both barcodes if large deletions are found
removeDels<-function(barcodeArray1="",barcodeArray2=""){
  #4 behaviours
    #non has - >nothing happens
    #1 has -> copy from 2
    #2 has -> copy from 1
    #both have -> copy from each other = massive xxxxxx

  large_del = "xxxx"
  aa= regexpr(barcodeArray1,pattern = 'xxxx+') #4 or more consecutive x
  #if find large deletion:
  if(aa[1]>0){
    substr(barcodeArray1,aa[1],aa[1]+attr(aa,"match.length")-1) <- substr(barcodeArray2,aa[1],aa[1]+attr(aa,"match.length")-1)
  }

  #at this point, barcodeArray1 is already modified to have whatever barcodeArray2 has:
  # if barcodeArray1 has del, now it has contect of array2
  # if barcodeArray2 ALSO has del, not barcodeArray1 has del (from 2) AND both will have del
  #
  aa= regexpr(barcodeArray2,pattern = 'xxxx+') #4 or more consecutive x
  #if find large deletion:
  if(aa[1]>0){
    substr(barcodeArray2,aa[1],aa[1]+attr(aa,"match.length")-1) <- substr(barcodeArray1,aa[1],aa[1]+attr(aa,"match.length")-1)
  }



  return(c(barcodeArray1,barcodeArray2))
}



manualDistML_2 <- function(barcodeLeaves,mu,alpha,nGen){

  removeLargeDels = F

  alphabet = array()
  alphabet[1]= "u"
  alphabet[2]= "r"
  alphabet[3]= "x"
  #mu & alpha indicate an explicit probabilistic model
  #Probability of no mutation for nGen cell divisions (mu is rate of edit per cell division)

  #transition probabilities

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

      barcodeLeaves1 = barcodeLeaves[i]
      barcodeLeaves2 = barcodeLeaves[j]

      if(removeLargeDels){

            #New function that will look for large deletions (more than 4 consecutive x)
            #The function will replace the xxxx to match whatever state the other cell has such that
            #those sites have minimum effect on the distance calculation
            #IF large_del in barcodeArray1 THEN copy what barcodeArray2 has
            maskedBarcodes = revomeDels(barcodeLeaves1,barcodeLeaves2) #barcodeArray2 is not modified
            barcodeLeaves1 = maskedBarcodes[1]
            barcodeLeaves2 = maskedBarcodes[2]

            #END section: large scale deletion
      }
      #after deletion check we can split now the barcode into a character array
      #and continue with the normal function
      barcodeArray1 =strsplit(barcodeLeaves1,"")[[1]]
      barcodeArray2 =strsplit(barcodeLeaves2,"")[[1]]


      distSum = 1
      ratio.sum =0
      ratio.product=1
      #for each pairwise comparison
      Pr.s1_s2=array()
      Pr.ind=array()
      for (s in 1:barcodeLength){ #may29 NOTE here is where we go site by site
        #new edit:
        #mu and alpha are vectors, so for each element we estimate the transition matrix and P(x,y)
        res.list =all.probs(as.numeric(mu[s]),as.numeric(alpha[s]),nGen,alphabet )
        PrMatrix = res.list[[1]]
        Pr_s0_array = res.list[[2]] #Pr of different letters in previous generation
        Tran.pr = res.list[[3]]
        #

        # PrMatrix = PrMatrix_list[[s]] #NOTE
        Pr.ind[s] = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])]
        #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
        # Sum{i=1}{3} Pr(s1|S0=i)Pr(s2|s0=i)Ps(S0=i)
        s1=which(alphabet==barcodeArray1[s])
        s2=which(alphabet==barcodeArray2[s])

        #sum over all possible S0 states
        sumPr=0
        #may29 NOTE here get a Tran.pr
        #Tran.pr = Tran.pr_list[[s]]
        for(a in 1:length(alphabet)){
           sumPr = sumPr+ Tran.pr[a,s1] * Tran.pr[a,s2] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
        }
        Pr.s1_s2[s] =sumPr
      } #end barcode loop
      #end may 29 NOTE



        ratio.product=sum(log(Pr.s1_s2))/sum(log(Pr.ind))
        # NOTE: remove the following line and everything will be OK
        # ratio.product=sum(log(Pr.s1_s2)) #/sum(log(Pr.ind))


      distMat[i,j]= distSum
      ratioMat[i,j]=1/ratio.sum *distSum
      productMat[i,j] = ratio.product

    }
  }
  return(productMat)
}


#The idea is to generate a new distance matrix
#based on the cousin distance Sep 4th

cousinDistML <- function(barcodeLeaves,mu,alpha,nGen){

  removeLargeDels = F

  alphabet = array()
  alphabet[1]= "u"
  alphabet[2]= "r"
  alphabet[3]= "x"
  #mu & alpha indicate an explicit probabilistic model
  #Probability of no mutation for nGen cell divisions (mu is rate of edit per cell division)

  #transition probabilities

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

      barcodeLeaves1 = barcodeLeaves[i]
      barcodeLeaves2 = barcodeLeaves[j]

      #after deletion check we can split now the barcode into a character array
      #and continue with the normal function
      barcodeArray1 =strsplit(barcodeLeaves1,"")[[1]]
      barcodeArray2 =strsplit(barcodeLeaves2,"")[[1]]


      distSum = 1
      ratio.sum =0
      ratio.product=1
      #for each pairwise comparison
      Pr.s1_s2=array()
      Pr.s1_s2__ = array()

      Pr.ind=array()
      for (s in 1:barcodeLength){ #may29 NOTE here is where we go site by site
        #new edit:
        #mu and alpha are vectors, so for each element we estimate the transition matrix and P(x,y)
        res.list =all.probs(as.numeric(mu[s]),as.numeric(alpha[s]),nGen,alphabet, cousin  = T)
        PrMatrix = res.list[[1]]
        Pr_s0_array = res.list[[2]] #Pr of different letters in previous generation
        Tran.pr = res.list[[3]]
        #

        # PrMatrix = PrMatrix_list[[s]] #NOTE
        Pr.ind[s] = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])]
        #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
        # Sum{i=1}{3} Pr(s1|S0=i)Pr(s2|s0=i)Ps(S0=i)
        s1=which(alphabet==barcodeArray1[s])
        s2=which(alphabet==barcodeArray2[s])

        #sum over all possible S0 states
        sumPr=0
        sumPr2=0


        for(a in 1:length(alphabet)){

            sumPr = sumPr + Pr_s0_array[a]  * cousinPr(alphabet[a], alphabet[s1], mu[s],alpha[s], nGen) * cousinPr(alphabet[a], alphabet[s2], mu[s],alpha[s], nGen)# -Tran.pr[a,s1] * Tran.pr[a,s2]

        #  for(b in 1:length(alphabet)){
        #    sumPr = sumPr + Pr_s0_array[a]  * Tran.pr[a,b] * Tran.pr[b,s1] * Tran.pr[b,s2]
        #  }
        #  sumPr2 = sumPr+ Tran.pr[a,s1] * Tran.pr[a,s2] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
        }


        Pr.s1_s2[s] =sumPr
      #  Pr.s1_s2__[s] = sumPr2
      } #end barcode loop
      #end may 29 NOTE


      #  ratio.product=sum(log(Pr.s1_s2))/sum(log(Pr.ind))
        ratio.product = sum(log(Pr.s1_s2))

      distMat[i,j]= distSum
      ratioMat[i,j]=1/ratio.sum *distSum
      productMat[i,j] = ratio.product

    }
  }
  return(productMat)
}
