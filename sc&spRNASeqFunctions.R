######################################
#     PACKAGES
######################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages('dplyr')
install.packages('Seurat')
install.packages("tidyverse")
#
######################################
#    LIBRARIES
######################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(dplyr)
library(Seurat)
library(tidyverse)
#
###########################################################################################################
#     Functions for 2-sample scRNASeq analysis
#     All functions require Seurat objects
#     All functions use only the "counts" slot and are independent of Seurat data transformations
###########################################################################################################
# Combined Chi2/weighted-t test (calls Chi2 test and IterWghtTest)
# To use Chi2 as a prefilter, set max.chi value at the prefilter threshold (must be <1)
# To use maximum value from Chi and weighted t-test as the final p-value keep default max.chi=1
# output is produced only for those genes that meet minimum requirements for both Chi2 and weighted t-test
# Use one of 4 icc= options: "i" weighted t-test with iterative icc (default option)
#                            "A" weighted t-test with ANOVA ICC (faster than "i" but may be less accurate)
#                            0 weighted t-test for aggregate counts
#                            1 common (unweighted) t-test
DGE.2samples <- function(object,features=NULL,ident.1=NULL,ident.2=NULL,
                         fc.thr=1,min.pct=0,max.chi=1,max.pval=1,min.count=30,
                         icc="i",df.correction=FALSE) {
  if (max.chi<1){                                                                                                  # runs Chi2Test when max.chi<1
    Chi2<-Chi2Test(object,features,ident.1,ident.2,fc.thr=1,min.pct,max.pval=max.chi,min.count)                    # Chi2 test
    features<-rownames(Chi2)                                                                                       # feature found by Chi2 test
    output<-IterWghtTtest(object,features,ident.1,ident.2,fc.thr,min.pct,max.pval,min.count,icc=icc,df.correction) # runs IterWghtTest
  }
  else{
    Chi2<-Chi2Test(object,features,ident.1,ident.2,fc.thr=1,min.pct,max.pval=1,min.count)
    IWT<-IterWghtTtest(object,features,ident.1,ident.2,fc.thr,min.pct,max.pval=1,min.count,icc=icc,df.correction)
    features<-intersect(rownames(Chi2),rownames(IWT))
    output<-IWT[features,]
    output[features,"wtd.t.p.value"]<-output[,2]
    output[features,"Chi2.p.value"]<-Chi2[features,2]
    for (gene in features) {output[gene,2]<-max(IWT[gene,2],Chi2[gene,2])}
    output<-output[output[,2]<=max.pval,]
  }
  if (max.chi<1){                                                                                                  # adds Chi2 p-values when max.chi<1
    output.features<-rownames(output)
    output[["Chi2.p.value"]]<-Chi2[output.features,]$p.value}
  return(output)
}
#Chi2 test
Chi2Test <- function(object,features=NULL,ident.1=NULL,ident.2=NULL,fc.thr=1,min.pct=0,max.pval=1,min.count=30) {
  if (is.null(features)){
    GL <- rownames(object)
  } else {
    GL <- as.vector(features)
  }
  if (is.null(ident.1)|is.null(ident.2)) {stop("Two identities in a Seurat object must be defined using ident.1= and ident.2=")}
  output<-data.frame(log2FC=numeric(),p.value=numeric(),col3=numeric(),col4=numeric())                               # output dataframe
  colnames(output)<-c("log2FC","p.value",str_c("Counts/Cell.",ident.1),str_c("Counts/Cell.",ident.2))
  object.1 <- subset(object, idents = ident.1)                                                                       # object subsetting
  object.2 <- subset(object, idents = ident.2)
  Ci.1 <- as.matrix(object.1[["RNA"]]$counts)                                                                        # count matrices
  Ci.2 <- as.matrix(object.2[["RNA"]]$counts)
  Nc.1 <- ncol(Ci.1)                                                                                                 # number of cells
  Nc.2 <- ncol(Ci.2)
  AC.1 <- rowSums(Ci.1)                                                                                              # aggregate counts /gene
  AC.2 <- rowSums(Ci.2)
  TC.1 <- sum(AC.1)                                                                                                  # total counts /sample
  TC.2 <- sum(AC.2)
  rowi<-0
  print("Conducting Chi2 testing:")
  pb <- txtProgressBar(min = 0, max = length(GL), initial = 0, style = 3)                                             # initialize progress bar
  for (rownum in c(1:length(GL))) {                                                                                  # main loop, gene by gene
    if ((AC.1[GL[rownum]] >= min.count | AC.2[GL[rownum]] >= min.count)&                                             # checking for min.count and min.pct
        ((sum(Ci.1[GL[rownum],]!=0)/Nc.1>min.pct)|(sum(Ci.2[GL[rownum],]!=0)/Nc.2>min.pct))){
      ContTable <- matrix(c(TC.1 - AC.1[GL[rownum]], TC.2 - AC.2[GL[rownum]],                                        # count table for chi2 test
                            AC.1[GL[rownum]], AC.2[GL[rownum]]), nrow = 2, ncol = 2)
      fc <- as.numeric((AC.1[GL[rownum]]/TC.1)/(AC.2[GL[rownum]]/TC.2))                                              # fold change (FC)
      if(!is.na(fc)){
        if (fc >= fc.thr | fc <= 1/fc.thr) {                                                                         # checking FC threshold
          chisquare <- as.numeric(chisq.test(ContTable)$p.value)                                                     # chi2 test
          if (chisquare <= max.pval) {                                                                               # checking max.pval threshold
            rowi<-rowi+1
            output[[rowi,1]] <- as.numeric(format(log(fc,2), digits=4, scientific=FALSE))                            # output assembly
            output[[rowi,2]] <- as.numeric(format(chisquare, digits=4, scientific=FALSE))
            output[[rowi,3]] <- as.numeric(format(AC.1[GL[rownum]]/Nc.1, digits=4, scientific=FALSE))
            output[[rowi,4]] <- as.numeric(format(AC.2[GL[rownum]]/Nc.2, digits=4, scientific=FALSE))
            row.names(output)[rowi]<-GL[rownum]
          }  
        }
      }
    }
    setTxtProgressBar(pb, rownum)                                                                                    # progress bar update
  }
  close(pb)                                                                                                          # close progress bar
  return(output)
}
#Weighted t-test with iterative weight calculation
IterWghtTtest <- function(object,features=NULL,ident.1=NULL,ident.2=NULL,fc.thr=1,min.pct=0,max.pval=1,min.count=30,
                          icc="i",df.correction=FALSE) {
  if (is.null(features)){
    GeneList <- rownames(object)
  } else {
    GeneList <- as.vector(features)
  }
  if (is.null(ident.1)|is.null(ident.2)) {stop("Two identities in a Seurat object must be defined using ident.1= and ident.2=")}
  output<-data.frame(log2FC=numeric(),p.value=numeric(), col3=numeric(),col4=numeric())                   # output dataframe
  colnames(output)<-c("log2FC","p.value", str_c("Counts/Cell.",ident.1),str_c("Counts/Cell.",ident.2))   # output dataframe
  object.1 <- subset(object, idents = ident.1)                                                           # ident.1 object
  Ci.1 <- as.matrix(object.1[["RNA"]]$counts)                                                            # ident.1 count matrix
  Ni.1 <- colSums(Ci.1)                                                                                  # ident.1 Ni vector
  Nc.1 <- ncol(Ci.1)                                                                                     # ident.1 number of cells
  Xi.1 <- Ci.1                                                                                           # normalized counts initiation     
  for (i in c(1:nrow(Ci.1))) {Xi.1[i,]=Ci.1[i,]/Ni.1}                                                    # ident.1 normalized counts
  object.2 <- subset(object, idents = ident.2)
  Ci.2 <- as.matrix(object.2[["RNA"]]$counts)
  Ni.2 <- colSums(Ci.2)
  Nc.2 <- ncol(Ci.2)
  Xi.2 <- Ci.2
  for (i in c(1:nrow(Ci.2))) {Xi.2[i,]=Ci.2[i,]/Ni.2}                                      
  rowi<-0                                                                                                # output row counter
  print("Conducting t-test:")
  pb <- txtProgressBar(min = 0, max = length(GeneList), initial = 0, style = 3)                          # initialize progress bar
  for (rownum in c(1:length(GeneList))) {                                                                # Main test loop
    AC.1 <- sum(Ci.1[GeneList[rownum],])                                                                 # ident.1 aggregated counts
    AC.2 <- sum(Ci.2[GeneList[rownum],])                                                                 # ident.2 aggregated counts
    Xi1<-Xi.1[GeneList[rownum],]                                                                         # ident.1 normalized counts
    Xi2<-Xi.2[GeneList[rownum],]                                                                         # ident.2 normalized counts
    if ((AC.1 >= min.count | AC.2 >= min.count)&                                                         # checking for minumum aggregated counts
        ((sum(Xi1!=0)/Nc.1>min.pct)|(sum(Xi2!=0)/Nc.2>min.pct))){                                        # checking for minumum expression
      wi.1<-ICCWeight(h=Ci.1[GeneList[rownum],],n=Ni.1,icc=icc)                                          # ident.1 weights
      wi.2<-ICCWeight(h=Ci.2[GeneList[rownum],],n=Ni.2,icc=icc)                                          # ident.2 weights
      fc <- sum(Xi1*wi.1)/sum(Xi2*wi.2)                                                                  # fold change
      if (!is.na(fc)){                                                                                   # removing 0/0 division
        if ((fc >= fc.thr | fc <= 1/fc.thr)&                                                             # checking FC threshold
            (sum(Xi1!=0)>=3|sum(Xi2!=0)>=3)) {                                                           # checking for at least 3 nonzero values in ident.1 or ident.2                                                      
          if(df.correction) {wTtest <- as.numeric(alt.wttest2(Xi1,Xi2,wi.1,wi.2))}                       # weighted t-test with effective df
          else {wTtest <- as.numeric(alt.wttest(Xi1,Xi2,wi.1,wi.2))}                                     # weighted t-test without df correction
          if (wTtest <= max.pval) {                                                                      # checking p-value threshold    
            rowi<-rowi+1                                                                                 # output row counter
            output[[rowi,1]] <- as.numeric(format(log(fc,2), digits=4, scientific=FALSE))                # output
            output[[rowi,2]] <- as.numeric(format(wTtest, digits=4, scientific=FALSE))
            output[[rowi,3]] <- as.numeric(format(AC.1/Nc.1, digits=4, scientific=FALSE))
            output[[rowi,4]] <- as.numeric(format(AC.2/Nc.2, digits=4, scientific=FALSE))
            row.names(output)[rowi]<-GeneList[rownum]
          }  
        } 
      }
    }
    setTxtProgressBar(pb, rownum)                                                                        # progress bar update
  }
  close(pb)                                                                                              # close progress bar
  return(output)
}
alt.wttest <- function(x1, x2, w1, w2) {
  # alternative weighted t-test based on Margolin-Leikin variance estimator
  stopifnot(length(x1)==length(w1) && length(x2)==length(w2))
  w1 = w1/sum(w1)
  m1 = sum(x1*w1)
  vm1 = sum(w1^2*(x1-m1)^2) / (1-sum(w1^2)) # unbiased when w ~ 1/s^2 (and sum(w)=1)
  w2 = w2/sum(w2)
  m2 = sum(x2*w2)
  vm2 = sum(w2^2*(x2-m2)^2) / (1-sum(w2^2)) # unbiased when w ~ 1/s^2 (and sum(w)=1)
  s12 = sqrt(vm1 + vm2)
  t = (m1 - m2) / s12
  df = s12^4 / (vm1^2/(length(x1)-1) + vm2^2/(length(x2)-1)) # not sure it's right ...
  p = 2*pt(-abs(t), df=df)
  return(p)
}
alt.wttest2 <- function(x1, x2, w1, w2) {
  # alternative weighted t-test based on Margolin-Leikin variance estimator and effective degrees of freedom
  stopifnot(length(x1)==length(w1) && length(x2)==length(w2))
  w1 = w1/sum(w1)
  n1 = 1/sum(w1^2)
  m1 = sum(x1*w1)
  vm1 = sum(w1^2*(x1-m1)^2) / (1-1/n1) # unbiased when w ~ 1/s^2 (and sum(w)=1)
  w2 = w2/sum(w2)
  n2 = 1/sum(w2^2)
  m2 = sum(x2*w2)
  vm2 = sum(w2^2*(x2-m2)^2) / (1-1/n2) # unbiased when w ~ 1/s^2 (and sum(w)=1)
  s12 = sqrt(vm1 + vm2)
  t = (m1 - m2) / s12
  df = s12^4 / (vm1^2/(n1-1) + vm2^2/(n2-1)) # not sure it's right ...
  p = 2*pt(-abs(t), df=df)
  return(p)
}
#Calculation of weights based on ICC
ICC.AN<-function(h,n){
  N=sum(n)                                                    #ANOVA ICC calculation
  k=length(n)
  n0=(1/(k-1))*(N-sum(n^2/N))
  MSw=(1/(N-k))*(sum(h)-sum(h^2/n))
  MSb=(1/(k-1))*(sum(h^2/n)-(1/N)*(sum(h))^2)
  if ((MSb+(n0-1)*MSw)==0) {ICC=0}                            #setting ICC=0 when ANOVA denominator = 0
  else ICC=(MSb-MSw)/(MSb+(n0-1)*MSw)
  if(ICC < 0) ICC=0                                           #resetting negative ICC to 0
  if(ICC > 1) ICC=1                                           #resetting ICC>1 to ICC=1
  return(ICC)
}
ICC.iter<-function(h,n){                                      #Iterative ICC calculation
  x=h/n
  w0<-n/sum(n)                                              # initial weights 
  x0<-sum(x*w0)                                             # initial weighted average count
  VarT0<-x0*(1-x0)/sum(n)                             # initial variance @ icc=0
  VarE0<-sum(w0^2*(x-x0)^2)/(1-sum(w0^2))                   # initial measured variance with w0 weights
  if (VarE0<=VarT0){icc=0}    
  else{
    f <- function(icc,h=h,n=n) {
      x = h/n                                                   #normalized counts
      wprop = n/(1 + icc*(n-1))                                 #proportional weights
      w <- wprop/sum(wprop)                                     #normalized weights
      x1 = sum(x*w)                                             #weighted average
      VarT<-x1*(1-x1)/(sum(wprop))                              #VarT 
      VarE<-sum(w^2*(x-x1)^2)/(1-sum(w^2))                      #VarE
      VarE-VarT                                                 #VarE-VarT
    }
    ur = uniroot(f, 0:1, n=n, h=h, check.conv=T, extendInt="downX", tol = 1e-4/max(n))
    icc = ur$root
    if(icc > 1) icc=1                                           #resetting ICC>1 to ICC=1
  }
  return(icc)
}
ICCWeight <- function(h,n,icc="i") {                     #Calculation of weights based on ICC
  Nc<-length(n)
  if(length(h)!=Nc) {stop("Unequal lengths of Ni and Nzi vectors")} 
  if(Nc<3) {stop("At least 3 cells are required in each ident")}
  if (sum(h!=0)<3){w<-rep(1/Nc,Nc)} #exit and return equal weights if the number cells/samples with nonzero counts is less than 3
  else {
    if(icc=="i"){icc=ICC.iter(h,n)} else {
      if(icc=="A"){icc=ICC.AN(h,n)} else {
        if (icc==0) {icc=0} else {
          if (icc==1) {icc=1} else {
            stop("Invalid icc, must be icc = 'i', 'A', 0, or 1")
          }
        }  
      }
    } 
    wprop = n/(1 + icc*(n-1))                                 #proportional weights
    w <- wprop/sum(wprop)
  }
  return(w)
}
###########################################################################################################
#     Functions for multi-sample scRNASeq analysis
#
#     Recommended usage for scRNASeq: results<-DGE.MultiSample(object, samples.1=c(...), samples.2=c(...))
#     Recommended usage for spRNASeq: results<-DGE.MultiSample(object, samples.1=c(...), samples.2=c(...), icc=0)
#     Recommended usage for all assay: reassemble Seurat object(s) from raw counts after QC by using CreateSeuratObject()
# 
#     Requirements: At least 6 active identities within the object corresponding to different samples
#                   No other active identities
#                   At least 3 samples in samples.1 and 3 samples in samples.2
#                   A single layer for counts in the Seurat object (use JoinLayers function when merging data)
#                   Raw data without any transformations must be stored in the counts slot of Seurat object to ensure
#                   reliable results
#                   Active assay in the Seurat object must be "RNA" (default assay when using CreateSeuratObject())
#
#     Optional parameters:
#     features=c("A","B",...) vector of gene names to be analysed, default - all genes in object
#     min.pct= (default=0.03) limits analysis to only those genes for which at least one of the 2 sample groups contains
#     at least 3 samples with the fraction of cells expressing the gene larger than min.pct
#     fc.thr= (default=1) limits analysis to only those genes for which FC>fc.thr or (1/FC)>fc.thr; FC is the fold-change
#     max.pval= (default=1) limits output only to the genes with p-value <= max.pval
#     t.test= (defalut=FALSE) when TRUE performs common (unweighted) t-test for average counts
#     icc= (default="i", other options "A",0,1); icc="i" calculates ICC iteratively for count averaging 
#          icc="A" uses ANOVA ICC; icc=0 performs count aggregation (use for spRNASeq); icc=1 performs unweighted count averaging
#     df.correction (default=FALSE) when TRUE uses effective degrees of freedom (1/sum(w^2)) in the weighted t-test
#
#     All functions use only the "counts" slot and are independent of Seurat data transformations
###########################################################################################################
# Count averaging and variance calculation
# Usage: modified.object<-CntAv(object)
CntAv<-function(object,features=NULL,icc="i"){
  if (is.null(features)){                                                                # check for genes to be analyzed
    GeneList <- rownames(object)
  } else {
    GeneList <- as.vector(features)
  }
  object@misc$AV.data<-list()                                                            # creating AV.data list in misc slot of the Seurat object 
  Samples<-levels(object@active.ident)                                                   # retrieving sample identities from the object
  Ns<-length(Samples)                                                                    # number of samples in the object
  if(Ns<6) {stop("At least 6 samples with different active identities are required")}    # check for sample number
  Sstats<-matrix(nrow=3,ncol=Ns)
  colnames(Sstats)<-Samples
  rownames(Sstats)<-c("N.cells","N.counts","Counts/cell")
  Ng<-length(GeneList)                                                                   # number of genes
  print("Averaging counts and calculating variance:")
  pb <- txtProgressBar(min = 0, max = Ns * Ng, initial = 0, style = 3)                   # initialize progress bar
  for (j in c(1:Ns)){                                                                    # main loop for calculating weighted average and variance sample by sample
    sample=Samples[j]
    OBJ<-subset(object, idents=sample)                                                   # sample subsetting
    Ci<-as.matrix(OBJ[["RNA"]]$counts)                                                   # count matrix
    Ni<-colSums(Ci)                                                                      # total counts vector (Ni)
    Nc<-ncol(Ci)                                                                         # number of cells
    Sstats["N.cells",j]<-Nc
    Sstats["N.counts",j]<-sum(Ni)
    Sstats["Counts/cell",j]<-Sstats["N.counts",j]/Sstats["N.cells",j]
    AV<-vector()                                                                         # vector initialization
    VAR<-vector()
    PCT<-vector()
    for (i in c(1:Ng)){                                                                  # gene by gene calculation subloop 
      w<-ICCWeight(h=Ci[GeneList[i],],n=Ni,icc=icc)                                                 # statistical weights
      x<-as.vector(100*Ci[GeneList[i],]/Ni)                                              # normalized genes counts
      AV[i]<-sum(x*w)                                                                    # weighted average gene count (% UMI)
      VAR[i]<-if(AV[i]!=0){sum(w^2*(x-AV[i])^2)/(1-sum(w^2))} else{1/(sum(Ni)^2)}        # variance of the weighted average or 1/sum(Ni)^2 when all x=0
      PCT[i]<-sum(Ci[GeneList[i],]!=0)/Nc                                                # fraction of cells expressing the gene
      setTxtProgressBar(pb, i + j * Ng)                                                  # progress bar update
    }
    AVdat<-cbind(AV,VAR,PCT)                                                             # output matrix assembly
    rownames(AVdat)<-GeneList
    object@misc$AV.data[[sample]]<-as.data.frame(AVdat)                                  # writing output matrix into AV.data list within misc slot
  }
  close(pb)                                                                              # close progress bar
  object@misc$AV.data[["Sstats"]]<-as.data.frame(Sstats)
  return(object)
}
#
# Sample matrix generation
# Usage: sample.matrix<-SampleMatrix(object, samples.1=c(...),samples.2=c(...))
# generates matrix of weighted average counts, variances, and expression fractions for subsequent analysis
SampleMatrix<-function(object,samples.1=NULL,samples.2=NULL){
  if(is.null(samples.1)|is.null(samples.2)) {stop("Define sample name vectors using samples.1= and samples.2=")}
  N.1=length(samples.1)                                                                  # number of samples in group 1
  N.2=length(samples.2)                                                                  # number of samples in group 2
  if(N.1<3|N.2<3) {stop("At least 3 samples per group are required")}                    # checking for at least 3 samples per group
  Ng=nrow(object@misc$AV.data[[samples.1[1]]])                                           # number of genes
  if (is.null(Ng)) {stop("Cannot find average data matrix in the object. Please run CntAv()")}
  AV.1<-matrix(0,nrow=Ng,ncol=N.1)                                                       # initialization of sample matrices
  colnames(AV.1)=samples.1
  rownames(AV.1)=rownames(object@misc$AV.data[[samples.1[1]]])
  VAR.1<-AV.1
  PCT.1<-AV.1
  AV.2<-matrix(0,nrow=Ng,ncol=N.2)
  colnames(AV.2)=samples.2
  rownames(AV.2)=rownames(AV.1)
  VAR.2<-AV.2
  PCT.2<-AV.2
  print("Generating sample matrix:")
  pb <- txtProgressBar(min = 0, max = N.1 + N.2, initial = 0, style = 3)                 # initialize progress bar
  for (i in c(1:N.1)) {                                                                  # reading object@misc$AV.data created by CntAv()
    AV.1[,i]<-object@misc$AV.data[[samples.1[i]]]$AV                                     # into sample matrices for groups 1 and 2
    VAR.1[,i]<-object@misc$AV.data[[samples.1[i]]]$VAR
    PCT.1[,i]<-object@misc$AV.data[[samples.1[i]]]$PCT
    setTxtProgressBar(pb, i)                                                             # progress bar update
  }
  for (i in c(1:N.2)) {
    AV.2[,i]<-object@misc$AV.data[[samples.2[i]]]$AV
    VAR.2[,i]<-object@misc$AV.data[[samples.2[i]]]$VAR
    PCT.2[,i]<-object@misc$AV.data[[samples.2[i]]]$PCT
    setTxtProgressBar(pb, i + N.1)                                                       # progress bar update
  }
  close(pb)                                                                              # close progress bar
  Sstats<-object@misc$AV.data[["Sstats"]]
  samples.matrix<-list(AV.1,VAR.1,AV.2,VAR.2,PCT.1,PCT.2,Sstats=Sstats)                                # assembling a single list of sample matrices for
  return(samples.matrix)                                                                 # for subsequent analysis
}
#
# Weight calculation for average,variance vectors
# used by WT.MultiSample() function
IterVar<-function(Av,Va){                                                                
  Nv<-length(Va)                                                                         # lengths of average and variance vectors
  if(length(Av)!=Nv) {stop("Unequal lengths of average and variance vectors")} 
  if(Nv<3) {stop("At least 3 samples are required")}
  x<-as.vector(Av)                                                                     # average counts vector
  Vari<-as.vector(Va)                                                                  # variance vector  
  Nc<-length(x)                                                                        # number of remaining samples
  VarT0=1/sum(1/Vari)
  w0=(1/Vari)/sum(1/Vari)
  x0=sum(x*w0)
  VarE0<-sum(w0^2*(x-x0)^2)/(1-sum(w0^2))
  if(VarE0<=VarT0) {Vp=0}
  else{
    f <- function(Vp,x=x,Vari=Vari) {
      wprop<-1/(Vari+Vp)                                                            #proportional weights
      VarT<-1/sum(wprop)      
      w <- wprop/sum(wprop)                                                            #normalized weights
      x1 = sum(x*w)                                                                    #weighted average
      VarE<-sum(w^2*(x-x1)^2)/(1-sum(w^2))                                             #VarE
      VarE-VarT                                                                
    }
    Vp1 = max(Vari)
    ur = uniroot(f, c(0,Vp1),x=x,Vari=Vari,check.conv=T,extendInt="downX", tol = min(Vari)*1e-4)
    Vp = ur$root
  }
  if(Vp<0){Vp=0}
  w=(1/(Vari+Vp))/sum(1/(Vari+Vp))
  return(w)
}
#
# Weighted t-test for multiple samples  
WT.MultiSample <- function(sample.matrix,features=NULL,t.test=FALSE ,min.pct=0.03,fc.thr=1,max.pval=1,df.correction=FALSE) {
  AV1<-sample.matrix[[1]]                                                            # read average expression matrix
  var1<-sample.matrix[[2]]                                                           # read variance matrix
  AV2<-sample.matrix[[3]]
  var2<-sample.matrix[[4]]
  pct1<-sample.matrix[[5]]                                                           # read expression fraction matrix
  pct2<-sample.matrix[[6]]
  Sstats<-sample.matrix[["Sstats"]]
  if (is.null(features)){                                                            # check for genes to be analyzed
    GeneList <- rownames(AV1)
  } else {
    GeneList <- as.vector(features)
  }
  data1<-as.matrix(AV1)                                                              # setting up data matrices
  data2<-as.matrix(AV2)                                                              # 2 line below set up output dataframe
  output<-data.frame(p.value=numeric(),log2FC=numeric(),Av1=numeric(),Sd1=numeric(),C1=numeric(),Av2=numeric(),Sd2=numeric(),C2=numeric())
  colnames(output)<-c("log2FC","p.value","Wtd.%UMI.1","Sd.%UMI.1","Av.min.pct.1","Wtd.%UMI.2","Sd.%UMI.2","Av.min.pct.2")
  rowi<-0                                                                            # output row counter
  print("Conducting t-test:")
  pb <- txtProgressBar(min = 0, max = length(GeneList), initial = 0, style = 3)      # initialize progress bar
  for (rownum in c(1:length(GeneList))) {                                            # gene by gene DGE calculation loop
    Xi1<-data1[GeneList[rownum],]                                                    # data vectors for groups 1 and 2  
    Xi2<-data2[GeneList[rownum],]
    VARi1<-var1[GeneList[rownum],]                                                   # variance vectors for groups 1 and 2
    VARi2<-var2[GeneList[rownum],]
    PCTi1<-pct1[GeneList[rownum],]                                                   # expression fraction vectors for groups 1 and 2 
    PCTi2<-pct2[GeneList[rownum],]
    N1<-length(Xi1)                                                                  # number of samples in group 1 included in the analysis
    N2<-length(Xi2)                                                                  # number of samples in group 2 included in the analysis
    if (N1 >= 3 & N2 >= 3 & (mean(PCTi1)>=min.pct|mean(PCTi2)>=min.pct)){            # selecting genes meeting minimum sample and expression requirements
      wi.1<-if(t.test==TRUE) {replicate(N1,1/N1)} else {IterVar(Xi1,VARi1)}          # calculating weights for groups 1 and 2
      wi.2<-if(t.test==TRUE) {replicate(N2,1/N2)} else {IterVar(Xi2,VARi2)} 
      Xi1av<-sum(Xi1*wi.1)                                                           # weighted average expression for groups 1 and 2
      Xi2av<-sum(Xi2*wi.2)
      Xi1sd<-sqrt(sum(wi.1^2*(Xi1-Xi1av)^2)/(1-sum(wi.1^2)))                         # standard deviations for the weighted averages
      Xi2sd<-sqrt(sum(wi.2^2*(Xi2-Xi2av)^2)/(1-sum(wi.2^2)))
      fc <- Xi1av/Xi2av                                                              # fold change (FC)
      if (!is.na(fc)){                                                               # removing 0/0 division
        if (fc >= fc.thr | fc <= 1/fc.thr) {                                         # checking FC threshold
          if(df.correction) {wTtest <- as.numeric(alt.wttest2(Xi1,Xi2,wi.1,wi.2))}   # weighted t-test with effective df
          else {wTtest <- as.numeric(alt.wttest(Xi1,Xi2,wi.1,wi.2))}                 # weighted t-test without df correction
          if (wTtest <= max.pval) {                                                  # output filtering for p-value
            rowi<-rowi+1
            output[[rowi,1]] <- as.numeric(format(log(fc,2), digits=4, scientific=FALSE))  #output
            output[[rowi,2]] <- as.numeric(format(wTtest, digits=4, scientific=FALSE))
            output[[rowi,3]] <- as.numeric(format(Xi1av, digits=4, scientific=FALSE))
            output[[rowi,4]] <- as.numeric(format(Xi1sd, digits=4, scientific=FALSE))
            output[[rowi,5]] <- as.numeric(format(mean(PCTi1), digits=4, scientific=FALSE))
            output[[rowi,6]] <- as.numeric(format(Xi2av, digits=4, scientific=FALSE))
            output[[rowi,7]] <- as.numeric(format(Xi2sd, digits=4, scientific=FALSE))
            output[[rowi,8]] <- as.numeric(format(mean(PCTi2), digits=4, scientific=FALSE))
            row.names(output)[rowi]<-GeneList[rownum]
          }  
        } 
      }
    }
    setTxtProgressBar(pb, rownum)                                                     # progress bar update
  }
  close(pb)                                                                           # close progress bar
  list(DGE=output,Sstats=Sstats)
}
#
#Simplified wrapper function for one-step analysis
DGE.MultiSample<-function(object,samples.1=NULL,samples.2=NULL,features=NULL,t.test=FALSE, min.pct=0.03,fc.thr=1,max.pval=1,
                          icc="i",df.correction=FALSE) {
  if(is.null(samples.1)|is.null(samples.2)) {stop("Define sample name vectors using samples.1= and samples.2=")}
  N.1=length(samples.1)                                                              # number of samples in each group
  N.2=length(samples.2)
  if(N.1<3|N.2<3) {stop("At least 3 samples per group are required")}                # early check for sample numbers
  object.m<-CntAv(object,features,icc=icc)                                                   # Count averaging within each sample
  SM<-SampleMatrix(object.m,samples.1,samples.2)                                     # Assembly of sample matrices
  output<-WT.MultiSample(SM,features,t.test,min.pct,fc.thr,max.pval,df.correction)                        # Analysis of sample matrices
  param<-as.character(c(t.test,min.pct,fc.thr,max.pval,icc,df.correction))
  names(param)<-c("t.test","min.pct","fc.thr","max.pval","icc","df.correction")
  list(DGE=output$DGE,Sstats=output$Sstats,parameters=param)
}
###########################################
