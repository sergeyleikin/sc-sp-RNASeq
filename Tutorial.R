#'*WELCOME TO THE TUTORIAL*
#
#'*UNDER CONSTRUCTION!*
#
# Written and maintained by Sergey Leikin, leikins@mail.nih.gov
#
#'*BRIEF GENERAL INTRODUCTION*
# 
# This tutorial contains 3 sections, the titles of which speak for themselves: 
# 1. BASIC TUTORIAL, lines 22 - 310 (highly recommended for all method users). Initial Release.
# 2. ADVANCED TUTORIAL, lines XXX - XXX (may be of interest to advanced users). Under Construction!
# 3. NOTES AND COMMENTS FOR ADVANCED USERS (may be of interest to mathematicians, physicists, and statisticians). Under Construction!
#
# Recommended folders are based on Windows notations since I use a Windows workstation. Please adjust accordingly for Mac and other operating systems
#
# Seurat v.5.4.0 or later is required for full tutorial code to work. Earlier Seurat versions are incompatible with code for import and processing of 
# Visium and VisiumHD data due to a change in the coordinate system implemented starting from Seurat 5.4.0.
#
# Our method is designed for analysis scRNASeq data and spatial RNASeq based on full transcriptome sequencing (including Visium and similar assays).
# It is NOT designed for spatial data based on Xenium, Vizgen, and similar RNA imaging assays, analysis of which requires a different approach.
#
#'*BASIC TUTORIAL*
#'
#'Chapter 1: Analysis of high-quality scRNASeq data (plenty of cells and transcipt counts in each cell group of interest), lines XXX - XXX
#'Chapter 2: Analysis of moderate-quality scRNASeq data (limited number of cells and/or transcript counts, e.g., low expression genes in high-quality data),
#'           lines XXX - XXX
#'Chapter 3: Analysis of Visium and VisiumHD data (not tested for other spatial RNASeq assays so far. Under construction!), lines XXX - XXX
#'Chapter 4: Multiple comparison correction basics (False Discovery Rate, FDR)
#'
#'*BASIC TUTORIAL. Chapter 1: scRNASeq*
#'
# I recommend RStudio for users like me, who are not programmers or R experts
# R experts, please do not expect professional or optimized R code. Instead, my goal is code that is easy to write and understand.
#
##################################################################################
#'*Section 1.0. Getting started with two-sample analysis*#
#
# Create a folder C:\scRNASeq
# Dowload Neuron.QC.rds and sc-spRNASeqFunctions_v2.R files into this folder
# Run the following code:
setwd("C:/scRNASeq") 
# to set C:\scRNASeq as your working directory
# Open the sc-spRNASeqFunctions_v2.R file and run the entire code in the file to make sure that necessary libraries and functions are loaded.
# Install missing packages and libraries are directed by RStudio warnings.
# Load a pre-assembled Seurat object containing 10X Genomics scRNASeq data for E18 mouse neurons (see our paper) 
Neuron.QC=readRDS("Neuron.QC.rds") 
# Neuron.QC contains raw (unnormalized) scRNASeq counts of glutamatergic (active identity "GLUT") and GABAergic neurons (active identity "GABA")
#
##################################################################################
#'*Section 1.1. Two-sample analysis: Identifying and labeling cell groups/clusters*
#
# There are multiple approaches to separating cells within each scRNASeq dataset into groups referred to as identities within Seurat. Please 
# see Seurat instructions for approaches recommended by Seurat. In my opinion, however, the most reliable approach is based
# on marker genes (see our paper for the discussion of Seurat clustering drawbacks for differential expression analysis). For instance, the first
# step to select cells based of expression of markers genes may be to normalize Neuron.QC using the Relative Counts (RC) normalization method.
Neuron.QCN=NormalizeData(Neuron.QC, normalization.method = "RC")
# A possible next step is to subset the object based on expression of the marker genes: 
GABA.QCN=subset(Neuron.QCN, subset = Gad1>0 & Gad2>0 & Slc6a1>0)  
GLUT.QCN=subset(Neuron.QCN, subset = Gls>0 & Grin2b>0 & Slc17a7>0)
# Note that the initial quality control of data (minimum counts and features per cell ans well as fraction of mitochondrial genes may be performed
# in exactly the same way, see Seurat instructions)
#     Because there are few cells that express all 6 marker genes in Neuron.QCN, the two datasets we generated contain a small number of the same cells.
# To eliminate these "ambiguous cells", we can use the following code (there are many ways to do it)
glutcells=setdiff(colnames(GLUT.QCN[["RNA"]]$counts),colnames(GABA.QCN[["RNA"]]$counts))
GABAcells=setdiff(colnames(GABA.QCN[["RNA"]]$counts),colnames(GLUT.QCN[["RNA"]]$counts))
GLUT.QC=subset(Neuron.QC, cells = glutcells)
GABA.QC=subset(Neuron.QC, cells = GABAcells)
# Examination of the objects we created reveals "GLUT" and "GABA" identities in both GLUT.QCN and GABA.QCN, just "GLUT" identity in GLUT.QC and just "GABA" 
# identity in GABA.QC. At this point, the two latter datasets can be merged. However, let as reassign cell identities in each of them
Idents(GLUT.QC)="glutcells"
Idents(GABA.QC)="GABAcells"
# This step ensures unique single cell identity within each object. It is important only if the identities in the initial objects are the same.
# We can now merge the two objects:
Neuron.QC2=JoinLayers(merge(GLUT.QC, y=GABA.QC))
#'*IMPORTANT*: 1. Always use the JoinLayers() function because our data analysis functions require a single counts layer
#              2. To speed up data analysis, the final object should contain no assays other than "RNA" and no layers other than "counts". Remove
#                 all other assays and layers. If needed, rename the default/active assay to "RNA" (see Seurat's RenameAssays() function).
#              3. Do not normalize the data once our object is prepared for our differential gene expression analysis. This will just slow
#                 down the analysis. Our functions will use only the counts layer and perform their own normalization anyway.
#
##################################################################################
#'*Section 1.2. Two-sample analysis: Differential gene expression (DGE) analysis*
#
# The object Neuron.QC2 we just created is ready for the analysis. The results of the analysis are produced by a single function
#
DGE.results=DGE.2samples(Neuron.QC2, ident.1 = "glutcells", ident.2 = "GABAcells")
#
# This code performed DGE analysis between "glutcells" and "GABAcells" using default settings. Before we consider optional analysis parameters and
# their optimal settings, let's discuss the analysis output and quality control parameters reported in it.
#
##################################################################################
#'*Section 1.3. Two-sample analysis: DGE analysis output*
#
# Take a look at the head of the DGE.result dataframe generated by the analysis:
head(DGE.results)
#
# The row names in this dataframe are the names of thes genes that were analyzed. The dataframe contains 11 columns:
# log2FC, p.value, Counts/Cell.1, Counts/Cell.2, N.+_cells.1, N.+_cells.2, Counts.1, Counts.2, Sum_w^2.1, Sum_w^2.2, and Chi2.p.value
# Here ".1" and ".2" indicate identities 1 and 2, respectively. 
#'*log2FC*       is the logarithm base 2 of the ratio of the weighted mean relative counts in the two identities. log2FC=0 means no DGE.
#'               log2FC=1 means 2-fold upregulation of the gene in identity 1 vs. identity 2. log2FC=n means 2^n-fold upregulation.
#'*p.value*      is the DGE p-value (probability of the observed counts if the gene expression were the same in the two identities)
#'*Chi2.p.value* is the DGE p-value calculated by the chi squared test, which can be used to reduce false positive findings caused by sampling bias. 
#                The chi squared test is optional. Please see our paper on when and how to use these p-values. 
#'*Quality Control Parameters*:
#'*Counts/cell*  average number of gene counts per cell in each identity
#'*N.+_cells*    number of cells in each identity that express the gene. The accuracy of weighted averaging decreases at smaller values of
#                this parameter. Common (non-weighted) averaging and the non-weighted t-test provide more reliable (although more conservative) estimates 
#                of log2FC and p.value when this parameter is below 5-10 in both identities.
#'*Counts*       aggregate count of the gene in each identity; must be large (we recommend ~30 or more) in at least one of the two
#'               identities for log2FC and p.value to be accurate.
#'*Sum_w^2*      sum of cell weight squares in each identity. Sum_w^2>0.4 indicates that just one or two cells in the identity dominate the weighted 
#                average, which is not good. Large values of N.+_cells x Sum_w^2 indicate very uneven distribution of gene counts, which not good either.
#
#'*DGE results are reported only for genes meeting quality control requirements defined by optional parameters of the DGE.2samples() function.*
# 
##################################################################################
#'*Section 1.4. Two-sample analysis: DGE analysis parameters*
#
#'*DGE.2samples() has the following optional parameters that may be used for quality control and filtering of the results*
#
#'*features*     list of genes to be analyzed (default features=NULL, all genes), e.g., features=c("Col1a1","Col1a2") cause analysis of just Col1a1 and Col1a2.
#'*fc.thr*       fold-change threshold filter for the analysis and reporting (default fc.thr=1, no filter), e.g., only genes upregulated or downregulated more  
#                than 2-fold are analyzed and reported at fc.thr=2.
#'*min.pct*      minimum fraction of cells that must express the gene in at least one of the two identities (default min.pct=0, no filter)
#'*min.count*    minimum aggregate count in at least one of the two identities for the gene to be analyzed and reported (default min.count=30).
#'*sum.w2*       the value of Sum_w^2 above which weighted averaging switches to non-weighted one (default sum.w2=0.4).
#'*min.cells*    minimum number of cells in at least one of the identities for the gene to be analyzed and reported (default min.cells=10)
#'*icc*          method of weighted data averaging (default icc="i", weighted averaging with iterative intracluster correlation coefficient [ICC]). Other possible 
#                values: icc="A" (averaging with ANOVA ICC, see Supplement 1 and tutorial section for advanced users, not recommended for basic users), 
#                icc=1 (averaging with ICC=1, same as common non-weighted averaging), icc=0 (averaging with ICC=0, same as aggregate count averaging, useful for VisiumHD)  
#'*df.correction* (default df.correction=FALSE), see NOTES AND COMMENTS FOR ADVANCED USERS.
#
##################################################################################
#'*Section 1.5. Two-sample analysis: Basic workarounds*
#
#'*Small number of cells expressing genes of interest*
# This is a common scRNASeq analysis problem for low expression genes regardless of the number of cells per identity. For instance,
DGE.results.mc5=DGE.2samples(Neuron.QC2, ident.1 = "glutcells", ident.2 = "GABAcells", min.cells=5, min.count = 10)
# contains results for almost almost 2,000 more genes than DGE.Results calculated by DGE.2samples() with default min.cells=10 and min.count=30
diff.genes=setdiff(rownames(DGE.results.mc5),rownames(DGE.results))
# identifies all these "extra" low expression genes
DGE.results.dg=DGE.results.mc5[diff.genes,]
# is the dataframe containing just these low expression genes. Note that 183 of these genes are upregulated or downregulated more than 16-fold:
DGE.results.LE=DGE.results.dg[abs(DGE.results.dg$log2FC)>4,]
# but is this a real DGE? Some of these 183 findings could result from weighted averaging inaccuracy at low min.cell and min.count settings.
#   To test which of these 183 findings are most likely to be real, we reanalyze these genes by common non-weighted averaging and more conservative t-test
DGE.results.LEgenes=DGE.2samples(Neuron.QC2, features=rownames(DGE.results.LE), ident.1 = "glutcells" ,ident.2 = "GABAcells", min.cells=5, min.count = 10, icc=1)
# We next examine which of the genes have weighted log2FC and non-weighted log2FC within 30% of each other and therefore are not significantly affected by cell weights:
DGE.results.LEgenes.a=DGE.results.LEgenes[(abs(DGE.results.LE$log2FC*DGE.results.LEgenes$log2FC)>0&
                                             (abs(DGE.results.LE$log2FC-DGE.results.LEgenes$log2FC)<0.37|
                                                (abs(DGE.results.LE$log2FC)==Inf&abs(DGE.results.LEgenes$log2FC)==Inf))),]
# which turn out to be 174 of the 183 genes. Filtering for the p-value of the more conservative t-test
DGE.results.LEgenes.b=DGE.results.LEgenes.a[DGE.results.LEgenes.a$p.value<0.0001,]
# identifies 20 low expression genes with DGE p-value below 0.0001, most of which are indeed likely to have significant DGE. 
#     This strategy may be very useful when missing an potentially important low expression gene could affect the study. Of course, all such genes need to be validated
# by independent experiments, but one needs to know the candidate genes for validation...
#
#'*Small number of cells (<25-50) in one or both identities* 
# This is another scRNASeq problem, which is very common for analysis of rare cell types. It can be addressed by using the same strategy as we outlined above. The examples
# are under construction. Please stay tuned.
#
##################################################################################
#'*Section 1.6. Multiple sample analysis: Basic concepts*
#
# The main pitfall of 2-sample analysis is that individual cells are not really independent biological replicates in the context of most studies. This leads
# to artificially low p-values and therefore to potentially significant false findings. The only proper way to address this problem is to use as many real
# biological replicates (e.g., cell extracts from different animals) as one can afford (ideally 5 or more per genotype or treatment, but at least 3). We refer to 
# analysis of such data as multiple sample analysis. Please see our paper for more details on why it is important and more reliable than 2-sample analysis above.
# Here I describe only how to perform such analysis using our analysis functions.
#
#'*The key to multiple sample analysis is to assign a separate, individual identity to cells from each biological replicate. All replicates must be merged into*
#'*a single Seurat object, but their cells must have distinct identities* 
#
##################################################################################
#'*Section 1.7. Multiple sample analysis: Assembling Seurat object*
#
# Let's create a mock multiple sample experiment. First, we randomly split the 1079 glutamatergic neurons in GLUT.QC into 3 "biological replicates"
# using glutcells barcodes we extracted from data in Section 1.1
glutcells.r1=sample(glutcells,350)
temp=setdiff(glutcells,glutcells.r1)
glutcells.r2=sample(temp,350)
glutcells.r3=setdiff(temp,glutcells.r2)
Untr.r1=subset(GLUT.QC,cells=glutcells.r1)
Idents(Untr.r1)="untreated.1"
Untr.r2=subset(GLUT.QC,cells=glutcells.r2)
Idents(Untr.r2)="untreated.2"
Untr.r3=subset(GLUT.QC,cells=glutcells.r3)
Idents(Untr.r3)="untreated.3"
# This code created 3 "untreated" biological replicates of glutamatergic neurons with distinct identities.
# We now create 3 "treated" replicates by 50% upregulation of 1000 randomly selected genes
genes=rownames(GLUT.QC[["RNA"]]$counts)
modgenes=sample(genes,1000)
Tr.r1.mat=as.matrix(Untr.r1[["RNA"]]$counts)
Tr.r1.mat[modgenes,]=round(1.5*Tr.r1.mat[modgenes,])
Tr.r1=CreateSeuratObject(Tr.r1.mat)
Idents(Tr.r1)="treated.1"
Tr.r2.mat=as.matrix(Untr.r2[["RNA"]]$counts)
Tr.r2.mat[modgenes,]=round(1.5*Tr.r2.mat[modgenes,])
Tr.r2=CreateSeuratObject(Tr.r2.mat)
Idents(Tr.r2)="treated.2"
Tr.r3.mat=as.matrix(Untr.r3[["RNA"]]$counts)
Tr.r3.mat[modgenes,]=round(1.5*Tr.r3.mat[modgenes,])
Tr.r3=CreateSeuratObject(Tr.r3.mat)
Idents(Tr.r3)="treated.3"
# Now we can assemble a single Seurat object with 3 "untreated" replicates and 3 "treated" replicates
GLUT.MS=JoinLayers(merge(Untr.r1, y=c(Untr.r2,Untr.r3,Tr.r1,Tr.r2,Tr.r3)))
# This object contains 6 unique identities and is ready for multiple sample DGE analysis. Object from real multiple sample experiments should be assembled 
# in the same way (of course, without randomly splitting cells and mock "treatment").
#
##################################################################################
#'*Section 1.7. Multiple sample analysis: Simplest approach, slow parameter optimization*
#
# The easiest way to analyze multiple sample experiment is
DGE.multi.results=DGE.MultiSample(GLUT.MS, samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"))
# The output of this function is a list containing 3 items: DGE - dataframe of the analysis results, Sstats - sample statistics, parameters - parameters 
# used for the analysis. Let's take a look at the genes with DGE p-value below 0.05
DGE.ms.results=DGE.multi.results$DGE
DGE.ms.p05=DGE.ms.results[DGE.ms.results$p.value<0.05,]
DE.genes=rownames(DGE.ms.p05)
DE.genes.check=length(intersect(DE.genes,modgenes))
# The analysis found DGE in 529 of the 1000 genes we modified and DGE in 462 genes we did not modify. This is as good as we could expect because:
# 1. Modification of very low expression genes, which were not analyzed, would not be detected. Only 579 of the modified genes were analyzed
DE.genes.check2=length(intersect(rownames(DGE.ms.results),modgenes)) 
# 529 of these 579 were discovered with p<0.05, which is excellent considering a large fraction of low expression genes among the 579. 
# 2. Our simplistic gene modification caused a change in total gene count per cell, modifying relative expression of all other genes, which could result 
# in DGE of genes other than the ones we explicitly changed detectable with p<0.05. To test this, we can filter the results based on log2FC because the 
# "unintended" modification should lead to |log2FC|< 0.1
DGE.ms.p05.log2FC=DGE.ms.p05[abs(DGE.ms.p05$log2FC)>0.1,]
# This filter removed all 462 genes we did not explicitly modify and kept all 529 genes we did modify.
#     Let us now check the performance performance of DGE.MultiSample() with t.test parameter set as t.test=TRUE
DGE.multi.results=DGE.MultiSample(GLUT.MS, samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"), t.test=TRUE)
DGE.ms.results=DGE.multi.results$DGE
DGE.ms.p05=DGE.ms.results[DGE.ms.results$p.value<0.05,]
DE.genes=rownames(DGE.ms.p05)
DE.genes.check=length(intersect(DE.genes,modgenes))
DE.genes.check2=length(intersect(rownames(DGE.ms.results),modgenes))
# This analysis discovered DGE in 533 of the 573 analyzed genes we modified and DGE in 462 we did not explicitly modify. It performed better that the
# analysis with t.test=FALSE because having just 3 samples per sample group is not enough for accurate calculation of the weight of each sample. The
# t.test=TRUE setting replaces weighted averaging of samples and weighted t-test between the sample groups with common (non-weighted) averaging and t.test. 
# When the samples within groups have similar numbers of cells and similar sequencing depths, like in our example, the non-weighted averaging and t.test are
# expected to perform slightly better than the weighted ones.  
#
##################################################################################
#'*Section 1.8. Multiple sample analysis: DGE output*
#
# Before we proceed to other parameters of multiple sample DGE analysis, let us consider its output more closely
head(DGE.ms.p05)
# Compared to 2-sample analysis, this dataframe is "missing" the Chi2.p.value column and has 3 extra columns Wtd%UMI, Sd.%UMI, and Av.cell.fr
# for each of the two sample groups samples.1 and samples.2. The Chi2.p.value column is missing because because the chi squared test does not make
# sense for multiple sample analysis. Wtd.%UMI is the weigthed average expression of the gene across each sample group measured as % of all UMIs.
# Sd.%UMI is the standard deviation of Wtd.%UMI. Av.cell.fr. is the fraction of cells expressing the gene averaged over each sample group. Other
# columns are the same as in 2-sample analysis, except Counts and N.+_cells are reported for each sample in the group, separated by ":". 
#
##################################################################################
#'*Section 1.9. Multiple sample analysis: DGE analysis parameters*
#
#'*DGE.MultiSample() has the following optional parameters that may be used for quality control and filtering of the results*
#
#'*features*     same as in DGE.2samples (default features=NULL, all genes)
#'*fc.thr*       same as in DGE.2samples (default fc.thr=1, no filter)
#'*min.pct*      similar to DGE.2samples, except this is the threshold for Av.cell.fr (default min.pct=0.03, Av.cell.fr>0.03)
#'*min.count*    same as in DGE.2samples, except summed over all samples (default min.count=10). Does not filter results at default min.cells;
#'               may need to be increased at min.cells<10
#'*sum.w2*       similar to DGE.2samples but used for switching to non-weighted averaging for both cells within a sample and samples
#'               (default sum.w2=0.4).
#'*min.cells*    minimum number of cells in at least 3 samples within one of the two sample groups for the gene to be analyzed and reported
#'               (default min.cells=10)
#'*icc*          method of weighted data averaging within each sample (default icc="i", other possible values same as in DGE.2samples).   
#'*df.correction* (default df.correction=FALSE), see NOTES AND COMMENTS FOR ADVANCED USERS.
#'*t.test*       averaging and t-test method when comparing sample groups (default t.test=FALSE, alternative t.test=TRUE). At t.test=FALSE weighted
#'               averaging of the samples and weighted t-test are performed. At t.test=TRUE, common (non-weighted) averaging and t.test are performed
#
##################################################################################
#'*Section 1.10. Multiple sample analysis: Faster parameter optimization*
#
# DGE.MultiSample() is a wrapper function that requires significant memory and may perform slowly for very large datasets (see Supplement 1, Fig. S1.6
# in our paper). When parameter optimization (e.g., comparison of the results for t.test=FALSE and t.test=TRUE settings) is needed, a better approach is
# to perform the analysis in 3 steps:
# Step 1:
GLUT.MS=CntAv(GLUT.MS)
# At this step, weighted averaging is performed within each identity as specified by the icc= setting. The results are saved into the misc slot of GLUT.MS.
# Step 2:
GLUT.SampleMatrix=SampleMatrix(GLUT.MS,samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"))
# generates various matrices as described in the Readme file.
# Step 3: 
DGE.multi.results=WT.MultiSample(GLUT.SampleMatrix)
# generates the same results as DGE.MultiSample(). DGE.MultiSample() simply wraps these 3 steps together
#     The benefit of this 3-step approach are:
# 1. Adding weighted avarage values, their variances, and other parameters for each identity in the object to object's misc slot.
# 2. For very large datasets, step 1 may take hours. It has only 1 parameter, icc. Under most circumstances the default icc="i" is appropriate for scRNASeq.
#    Therefore, this step is required just once even when all other DGE.MultiSample() parameters need to be optimized. Steps 2 and 3 are much faster and shouldn't
#    take more than a few seconds to couple of minutes on a modern computer. The last step based on the Wt.MultiSample() function is the one that utilizes all the 
#    DGE.MultiSample() parameters except icc. Therefore, parameter optimization (next section) is much faster by using just the step 3.
#
##################################################################################
#'*Section 1.11. Multiple sample analysis: Parameter optimization*
#
#'*Small number of cells expressing genes of interest*
# Very low expression genes pose a problem for multiple sample analysis just like for the 2-sample one. For instance,
DGE.multi.results.mc5=WT.MultiSample(GLUT.SampleMatrix, min.pct=0, min.cells = 5, t.test=TRUE)
DGE.multi.results.p05.lg2FC01=DGE.multi.results.mc5$DGE[(DGE.multi.results.mc5$DGE$p.value<0.05&abs(DGE.multi.results.mc5$DGE$log2FC)>0.1),]
# finds 550 rather than 533 genes with DGE among the genes we modified
length(intersect(rownames(DGE.multi.results.p05.lg2FC01),modgenes))
# We can recover even more modified genes by further reducing min.cells and min.count as long as we keep t.test=TRUE. However, the returns diminish and statistics
# becomes truly questionable with such small numbers of cells and counts.
#     In my opinion, t.test=TRUE setting is the only reasonable one for min.cells<10. However, this is a tricky issue. I am still thinking and testing different
# options for optimizing multiple sample analysis parameters to better recover low expression genes without compromising the analysis accuracy too much. I am also
# not sure yet how to provide simple recommendations for t.test=FALSE or t.test=TRUE in general, when the number of samples per group is small (e.g., 3 or 4).
# So far, I am finding the choice to be anything but simple. Please stay tuned for updates. Please provide your suggestions as well.
#'
#'*BASIC TUTORIAL. Chapter 2: spatial RNASeq*
#'
# UNDER CONSTRUCTION!
#
#'
#'*BASIC TUTORIAL. Chapter 3: Multiple comparison corrections and False Discovery Rate (FDR)*
#'
# UNDER CONSTRUCTION!
#