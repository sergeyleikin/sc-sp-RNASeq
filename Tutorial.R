#'*UNBIASED ANALYSIS OF sc- and sp-RNASeq DATA IN SEURAT PACKAGE FOR R*
# 
#'*Tutorial written and maintained by Sergey Leikin, leikins@mail.nih.gov*
#
# Based on the manuscript and R functions referred to and deposited @ https://github.com/sergeyleikin/sc-sp-RNASeq 
#
#'Chapters and sections marked *UNDER CONSTRUCTION* have not been written yet. Please stay tuned for weekly/monthly updates.
#
# Seurat v.5.4.0 or later is required for full tutorial code to work. Earlier Seurat versions are incompatible with code for import and processing of 
# Visium and VisiumHD data due to a change in the coordinate system implemented starting from Seurat 5.4.0.
#
# Our method is designed for analysis scRNASeq data and spatial RNASeq based on full transcriptome sequencing (including Visium and similar assays).
# It is NOT designed for spatial data based on Xenium, Vizgen, and similar RNA imaging assays, analysis of which requires a different approach.
#
#################################################################################
#################################################################################
#'*TOC*
#'
#'*Lines       Content*
#' 52 - XXX    BASIC TUTORIAL  @ALL-USERS 
#' 54 - 291       Chapter 1. Single cell RNASeq
#' 60 - 73            Section 1.0. Getting started
#' 74 - 100           Section 1.1. Cell identities and assay requirements
#' 101 - 111          Section 1.2. 2-sample DGE analysis with DGE.2samples()
#' 112 - 132          Section 1.3. Understanding DGE.2samples() output
#' 133 - 145          Section 1.4. Understanding DGE.2samples() parameters
#' 146 - 199          Section 1.5. Seurat objects for multiple sample DGE analysis
#' 200 - 205          Section 1.6. Multiple sample DGE analysis with DGE.MultiSample()
#' 206 - 241          Section 1.7. Understanding DGE.MultiSample() output
#' 242 - 258          Section 1.8. DGE.MultiSample() parameters
#' 259 - 291          Section 1.9. 3-step multiple sample DGE analysis (recommended)
#' 292 - 302      Chapter 2. Quality control and parameter optimization
#' 303 - 342          Section 2.1  Optimizing 2-sample DGE analysis
#' 343 - 437          Section 2.2  Optimizing multiple sample DGE analysis
#' XX - XXX       Chapter 3. Spatial RNASeq *(UNDER CONSTRUCTION)*
#' XX - XXX           Section 3.1  Low resolution assays (Visium)
#' XX - XXX           Section 3.2  High resolution assays (VisiumHD)
#' XX - XXX           Section 3.3  Optimizing spRNASeq analysis
#' XX - XXX       Chapter 4. Multiple comparison corrections *(UNDER CONSTRUCTION)*
#' XX - XXX           Section 4.1  Family-wise error rate
#' XX - XXX           Section 4.2  False discovery rate
#' 
#' XX - XXX     ADVANCED TUTORIAL *(UNDER CONSTRUCTION)*
#' XX - XXX           Chapter 5. Weighted t-test vs. non-weighted t-test
#' XX - XXX           Chapter 6. Degrees of freedom in weighted t-test
#'
#' XX - XXX     PERSONAL THOUGHTS AND COMMENTS *(UNDER CONSTRUCTION)*
#' XX - XXX           Data analysis in physics vs. biology
#' XX - XXX           Occam's razor, maximum entropy methods, and logical maxims 
#' XX - XXX           Uses and abuses of generative AI in transcriptomic analysis
#
#################################################################################
#################################################################################
#'*BASIC TUTORIAL. Chapter 1: Single cell RNASeq*
#'
# I recommend RStudio for users like me, who are not programmers or R experts.
# If you are an R expert, please do not expect professional or optimized R code. Instead, my goal is code that is easy to write and understand.
#
##################################################################################
#'*Section 1.0. Getting started*
#
# Create a folder C:\scRNASeq (Windows OS notations, adjust accordingly for other OSs)
# Dowload Neuron.QC.rds and sc-spRNASeqFunctions_v2.R files into this folder
# Run the following code:
setwd("C:/scRNASeq") 
# to set C:\scRNASeq as your working directory
# Open the sc-spRNASeqFunctions_v2.R file and run the entire code in the file to make sure that necessary libraries and functions are loaded.
# Install missing packages and libraries as directed by RStudio warnings.
# Load a pre-assembled Seurat object containing 10X Genomics scRNASeq data for E18 mouse neurons (see our paper) 
Neuron.QC=readRDS("Neuron.QC.rds") 
# Neuron.QC contains raw scRNASeq counts of glutamatergic (active identity "GLUT") and GABAergic neurons (active identity "GABA")
#
##################################################################################
#'*Section 1.1. Cell identities and assay requirements*
#
# There are multiple approaches to separating cells within each scRNASeq dataset into groups (identities within Seurat) for DGE analysis. Please 
# see Seurat instructions for approaches recommended by Seurat. In my opinion, however, the most reliable approach is based
# on marker genes (see our paper for the discussion of Seurat clustering drawbacks for DGE analysis).
#
# Step 1: Normalize data using the Relative Counts (RC) normalization method. DO NOT use log-normalization (default in Seurat) 
Neuron.QCN=NormalizeData(Neuron.QC, normalization.method = "RC")
# Step 2: Quality control (QC) based on fraction of mitochondrial transcripts, number of counts, number of features, etc (Neuron.QC is QC pre-processed).
# Step 3: Subset cells based on relative expression of marker genes (the normalized expression threshold does not have to be zero):
GABA.QC=subset(Neuron.QCN, subset = Gad1>0 & Gad2>0 & Slc6a1>0 & (Gls==0 | Grin2b==0 | Slc17a7==0))  
GLUT.QC=subset(Neuron.QCN, subset = Gls>0 & Grin2b>0 & Slc17a7>0 & (Gad1==0 | Gad2==0 | Slc6a1==0))
# Step 4: Remove normalized data layers and objects that are no longer needed and only increase memory demand and slow down analysis
GABA.QC[["RNA"]]$data<-NULL
GLUT.QC[["RNA"]]$data<-NULL
rm(Neuron.QCN)
# Step 5: Set unique identities of the cells in GABA.QC and GLUT.QC and merge all cells into a single Seurat object (rebuilt Neuron.QC) 
Idents(GLUT.QC)="glutcells"
Idents(GABA.QC)="GABAcells"
Neuron.QC=JoinLayers(merge(GLUT.QC, y=GABA.QC))
#'*IMPORTANT*: 1. Always use the JoinLayers() function because our data analysis functions require a single counts layer
#              2. The active assay must be "RNA" for our functions to work (see Seurat's RenameAssays() function).
#'            *3. Use only original raw counts and NEVER counts pre-processed with SCTransform() or other functions.* 
#              4. SCTransform creates a new "SCT" assay and makes it active while preserving "RNA" assay. If your object has been SCT pre-processed,
#                 you can just make the "RNA" assay active and then remove the SCT assay (see Seurat instructions).
#
##################################################################################
#'*Section 1.2. 2-sample DGE analysis with DGE.2samples()*
#
# The object Neuron.QC we just reassembled is ready for the analysis. The results of the analysis are produced by a single function
DGE.results=DGE.2samples(Neuron.QC, ident.1 = "glutcells", ident.2 = "GABAcells")
# This code performed DGE analysis between "glutcells" and "GABAcells" using default settings. You may encounter situations in which you may want to group
# several identities together and compare the groups. In the following mock example, we compare all Neuron.QC cells to themselves
NULL.results=DGE.2samples(Neuron.QC, ident.1 = c("glutcells","GABAcells"), ident.2 = c("glutcells","GABAcells")) 
# You may combine as many identities as you wish, just keep in mind that it is not the best approach to multiple sample analysis.
# Before we discuss optional DGE.2samples() parameters and optimal settings, let's discuss the analysis output.
#
##################################################################################
#'*Section 1.3. Understanding DGE.2samples() output*
#
# Take a look at the structure of the DGE.result dataframe generated by the analysis:
head(DGE.results)
# The row names in this dataframe are the names of the genes that were analyzed. The dataframe contains 11 columns:
# log2FC, p.value, Counts/Cell.1, Counts/Cell.2, N.+_cells.1, N.+_cells.2, Counts.1, Counts.2, Sum_w^2.1, Sum_w^2.2, and Chi2.p.value
# In these columns, ".1" and ".2" indicate identities 1 and 2, respectively. 
#'*log2FC*       logarithm base 2 of the ratio of the weighted mean relative counts in the two identities. log2FC=0 means no DGE.
#'               log2FC=1 means 2-fold upregulation of the gene in identity 1 vs. identity 2. log2FC=n means 2^n-fold upregulation.
#'*p.value*      DGE p-value (probability of the observed counts under the "null hypothesis" of the same gene expression in the two identities).
#'*Chi2.p.value* Chi squared test p-value, which may be useful for identifying sampling bias effects in two-sample tests (please see manuscript).
#'*Counts/cell*  average number of gene counts per cell in each identity, *used for quality control, see Section 2.1*.
#'*N.+_cells*    number of cells in each identity that express the gene, *used for quality control, see Section 2.1*. 
#'*Counts*       aggregate count of the gene in each identity, *used for quality control, see Section 2.1*.
#'*Sum_w^2*      sum of cell weight squares in each identity, *used for advanced quality control, see ADVANCED TUTORIAL*.
#
#'*IMPORTANT:*
# Please note that DGE.results reports only 13,070 genes out of 21,091 genes with detected counts. Results for the "missing" genes have been 
# filtered out as unreliable based on the default settings of DGE.2samples() parameters (see Section 2.1 on optimizing these settings) 
# 
##################################################################################
#'*Section 1.4. Understanding DGE.2samples() parameters*
#
#'*features*     list of genes to be analyzed (default: features=NULL [all genes]), e.g., features=c("A","B","C") will cause analysis of just genes A, B, and C.
#'*fc.thr*       fold-change threshold filter (default: fc.thr=1 [no filter]), e.g., fc.thr=2 will cause analysis of just the genes with abs(log2FC)>1.
#'*min.pct*      minimum fraction of cells that must express the gene in at least one of the two identities (default: min.pct=0, [no filter]).
#'*min.count*    minimum aggregate count in at least one of the two identities (default: min.count=30).
#'*sum.w2*       the value of Sum_w^2 above which weighted averaging switches to non-weighted one (default: sum.w2=0.4).
#'*min.cells*    minimum number of cells in at least one of the identities for the gene to be analyzed and reported (default: min.cells=10).
#'*icc*          method of calculating cell weights (default: icc="i" [iterative weights]); icc="A" [ANOVA-based weights], icc=1 [equal weights], and icc=0 
#'               [weights proportional to total counts] are also allowed icc values (see Supplement 1 of the manuscript).
#'*df.correction* (default, df.correction=FALSE), see tutorial for advanced users.
#
##################################################################################
#'*Section 1.5. Seurat objects for multiple sample DGE analysis*
#
# The main pitfall of 2-sample analysis is that individual cells are not independent biological replicates in the context of most studies. This may cause 
# significant false positive findings. The traditional solution for this problem is the "pseudobulk" analysis of multiple (at least 3) samples per genotype/treatment.
# Multiple shortcomings of the pseudobulk approach (including analysis bias) are discussed in our manuscript. Below is the tutorial for an alternative, unbiased 
#'approach proposed in our manuscript. *IMPORTANT:* Our implementation of the latter approach requires all data being assembled into a single Seurat object.
#
# To demonstrate how to assemble a Seurat object meeting requirements for multiple sample analysis, in this section we create a mock multiple sample experiment and
# assemble the corresponding Seurat object. First, we randomly split the 1079 glutamatergic neurons in GLUT.QC into 3 "biological replicates"
glutcells=colnames(GLUT.QC[["RNA"]]$counts)
glutcells.r1=sample(glutcells,350)
glutcells.r2=sample(setdiff(glutcells,glutcells.r1),350)
glutcells.r3=setdiff(setdiff(glutcells,glutcells.r1),glutcells.r2)
Untr.r1=subset(GLUT.QC,cells=glutcells.r1)
Idents(Untr.r1)="untreated.1"
Untr.r2=subset(GLUT.QC,cells=glutcells.r2)
Idents(Untr.r2)="untreated.2"
Untr.r3=subset(GLUT.QC,cells=glutcells.r3)
Idents(Untr.r3)="untreated.3"
# This code created 3 "untreated" biological replicates of glutamatergic neurons with distinct identities.
# We now create 3 "treated" replicates by 50% upregulation of 100 randomly selected genes that have between 0.1 and 1 counts per cell
# Step 1: identifying genes with 0.1 to 1 counts per cell
GLUT.counts=as.matrix(GLUT.QC[["RNA"]]$counts)
GLUT.counts.01_1=GLUT.counts[(rowSums(GLUT.counts)/ncol(GLUT.counts)>0.1 & rowSums(GLUT.counts)/ncol(GLUT.counts)<1),]
genes=rownames(GLUT.counts.01_1)
# Step 2: randomly selecting 100 of these genes 
modgenes=sample(genes,100)
# Step 3: generating "treated" samples
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
# We generated 6 Seurat objects for glutamatergic neurons from 6 mock experiments (3 untreated and 3 treated). Each group of glutamatergic neuron has a unique identity 
# We can now assemble a single Seurat object with 3 "untreated" replicates and 3 "treated" replicates that meets the requirements for our multiple sample analysis
GLUT.MS=JoinLayers(merge(Untr.r1, y=c(Untr.r2,Untr.r3,Tr.r1,Tr.r2,Tr.r3)))
#    Real multiple sample experiments should be assembled into a single Seurat object in the same way. As long as all relevant groups of cells have unique active identities,  
# the Seurat object may contain hundreds or thousands of identities. The only requirement is that all cells within the same cell group (e.g., cell type) in a sample
# (i.e., biological replicate) have the same identity. For instance, if we were to generate 3 untreated and 3 treated groups of GABAergic neurons, we would give unique
# identities to all of them and then assemble a single Seurat object with both glutamatergic and GABAergic neurons (6 samples, 2 identities per sample). This would enable
# us to compare treated and untreated glutamatergic neurons as well as treated and untreated GABAegric neurons. However, if we would want to compare treated and untreated 
# glutamatergic+GABAergic neurons as a single group, we would need to change the identities of GABAergic neurons within each sample to match the corresponding identities 
# of glutamatergic neurons (so that they are grouped together). Alternatively, one may use multiple identifiers for each cell in Seurat and set active identity to different 
# identifiers depending on how one wants to group the cells (see Seurat instructions). The multiple sample DGE analysis is based on selecting relevant active identities 
# from a single Seurat object.
#
##################################################################################
#'*Section 1.6. Multiple sample DGE analysis with DGE.MultiSample()*
#
# The simplest (but not optimal) way to analyze multiple sample experiment is
DGE.multi.results=DGE.MultiSample(GLUT.MS, samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"))
#
##################################################################################
#'*Section 1.7. Understanding DGE.MultiSample() output*
#
# The output of DGE.MultiSample() is not a dataframe but a list containing 3 items: DGE, Sstats, and parameters 
# Item 1: DGE.multi.results$DGE is the dataframe of the analysis results
DGE.ms.results=DGE.multi.results$DGE
# Item 2: DGE.multi.results$Sstats is the statistics of the number of cells, total transcript counts, and average counts per cell in each sample
View(DGE.multi.results$Sstats)
# Item 3: DGE.multi.results$parameters is the analysis parameters
View(as.matrix(DGE.multi.results$parameters))

#'*DGE.MultiSample() performance*
# First let's make a smaller dataframe showing only the results with DGE p-value below 0.05
DGE.ms.p05=DGE.ms.results[DGE.ms.results$p.value<0.05,]
# This dataframe contains 98 of the 100 genes we modified and no false positive findings: 
DE.genes=rownames(DGE.ms.p05)
DE.genes.check=length(intersect(DE.genes,modgenes))
# The two "missing" genes have borderline p-values between 0.05 and 0.08:
m.genes=setdiff(modgenes,DE.genes)
View(DGE.ms.results[m.genes,])
# The DGE.MultiSample() function performed very well despite fairly low expression of genes we chose to modify.   
#
#'*DGE.MultiSample() output columns*
#
#'*log2FC*       same as 2-sample analysis.
#'*p.value*      same as 2-sample analysis.
#'*Wtd%UMI*      gene expression averaged across biological replicates (% of UMI counts per cell).
#'*Sd.%UMI*      standard deviation of Wtd%UMI.
#'*Av.cell.fr*   fraction of cells expressing the gene averaged over biological replicates.
#'*Counts/cell*  number of gene counts per cell averaged over biological replicates, *used for quality control, see Section 2.1*.
#'*N.+_cells*    number of cells in each sample that express the gene, separated by ":", *used for quality control, see Section 2.1*. 
#'*Counts*       aggregate count of the gene in sample, separated by ":", *used for quality control, see Section 2.1*.
#'*Sum_w^2*      sum of sample weight squares, *used for quality control, see t.test setting optimization*.
#
# Note that this dataframe is "missing" the Chi2.p.value column because the chi squared test is not useful multiple sample analysis.  
#
##################################################################################
#'*Section 1.8. DGE.MultiSample() parameters*
#
#'*DGE.MultiSample() has the following optional parameters that may be used for quality control and filtering of the results*
#
#'*features*     same as in DGE.2samples (default: features=NULL [all genes]).
#'*fc.thr*       same as in DGE.2samples (default: fc.thr=1 [no filter]).
#'*min.pct*      similar to DGE.2samples, except this is the threshold for Av.cell.fr (default: min.pct=0.03 [Av.cell.fr>0.03]).
#'*min.count*    same as in DGE.2samples, except summed over all samples (default: min.count=10). Does not filter results at default min.cells;
#'               may need to be increased at min.cells<10.
#'*sum.w2*       similar to DGE.2samples, used for switching to non-weighted averaging of cells within each sample and switching to non-weighted sample
#'               averaging (default: sum.w2=0.4).
#'*min.cells*    minimum number of cells in at least 3 samples within one of the two sample groups (default: min.cells=10).
#'*icc*          method for calculating cell weights within each sample (default: icc="i"), same as in DGE.2samples.   
#'*df.correction* (default: df.correction=FALSE), see ADVANCED TUTORIAL.
#'*t.test*       method for comparing sample groups (default: t.test=FALSE [weighted t-test used]; alternative: t.test=TRUE [non-weighted t-test used]). 
#
##################################################################################
#'*Section 1.9. 3-step multiple sample DGE analysis*
#
# DGE.MultiSample() is a wrapper function for 3 analysis steps. It makes the analysis a bit simpler but requires much more memory and may perform slowly
# for large datasets (see Supplement 1, Fig. S1.6 in our paper). In general, I recommend going through the analysis steps one by one instead of using DGE.MultiSample():
# The parameters and output will be exactly the same, but parameter optimization will be a lot faster (seconds to minutes instead of hours for large datasets) and your
# analysis will be better documented. I illustrate the steps below using the GLUT.MS Seurat object we created for the DGE.MultiSample() analysis.
#
#'*Step 1:* Weighted averaging of counts within each active identity of the Seurat object. Calculates weighted average counts, their variances and several other parameters.  
GLUT.MS=CntAv(GLUT.MS)
# I recommend using the same Seurat object name for the output and CntAv() argument. In this case, CntAV() simply adds its output to the misc slot of the Seurat object 
# instead of creating a new object. This step may take from minutes to hours, depending on the number of identities and cells in your object. I strongly recommend saving
# the resulting Seurat object, which you can then analyze in many different ways with each analysis taking less a minute on a good computer.  
saveRDS(GLUT.MS,"GLUT.MS.rds")
#'*IMPORTANT:* The only two parameters of this function are features= (genes to be analyzed) and icc= (averaging method). At this step, I recommend averaging all genes
#'(default). I recommend using the default icc="i" for analyzing scRNASeq data and icc=0 for analyzing VisiumHD data (see our manuscript). The value of icc defined in
#'DGE.MultiSample() is used just for this first step.  
#
#'*Step 2*: Assembly of sample matrix for further analysis. This is the step that defines which identities are used for DGE analysis
GLUT.SampleMatrix=SampleMatrix(GLUT.MS,samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"))
# The only two parameters of this function are samples.1= and samples.2=. At least 3 identities are required for each of the two groups of samples. The Seurat object
# used by this function must be the one generated by CntAv() function. The structure of the output matrix is described in the Readme.md file on out Github site.
#
#'*Step 3*: Analysis of DGE between the two sample groups in the sample matrix generated at Step 2. 
DGE.multi.results=WT.MultiSample(GLUT.SampleMatrix)
# The WT.MultiSample() function has the same optional parameters as the DGE.MultiSample() function, except for icc= (not used at this step). The selection of genes for 
# for the analysis can be made at this step. The output produced by WT.MultiSample() is identical to the output of DGE.MultiSample() when the same analysis parameters
# are used. 
#'*IMPORTANT:* DGE.MultiSample() will simplify your life only when: (a) You are sure that you will performs only one DGE analysis with your Seurat object (e.g., for
# just one cell type), (b) You are sure that you know the analysis parameters you want to use, and (c) You are sure that you will not want to make some adjustments 
# to analyze or reanalyze low expression genes. In all other cases, the 3-step analysis will save you a lot of time and headaches.
#
#####################################################################################
#####################################################################################
#'*Chapter 2. Quality control and parameter optimization*
# 
# Accuracy of any statistical analysis relies on having enough replicates and enough data points for each replicate
# This Chapter provides basic recipes for balancing the analysis accuracy with practical realities of scRNASeq (high cost per sample, rare cell types, low expression genes). 
# Main analysis parameters that can be adjusted to achieve the desired balance are: 
#'*min.pct*, *min.cells*, and *min.count*
#'they define minimum thresholds for *Av.cell.fr*, *N.+_cells*, and *Counts*, respectively. I find nonzero min.pct setting to be useful only for multiple sample experiments,
#'in which it ensures consistency between samples containing very different numbers of cells. Therefore, *Av.cell.fr* column is included only in the multiple 
#'sample analysis output, while *N.+_cells* and *Counts* are reported as quality control parameters for 2-sample analysis as well. 
#'
#'Additional settings that may affect/improve the balance between the analysis accuracy and data limitations in multiple sample experiments are:
#' *icc*, *t.test* (multiple experiments only), and  *Sum_w^2* (advanced users only). 
# 
##################################################################################
#'*Section 2.1 Optimizing 2-sample DGE analysis*
#
#'To increase analysis accuracy (at the cost of not being able to analyse lower expression genes), I recommend performing the analysis at default *min.cells=10* and
#'and *min.count=30* followed by filtering the results based on the values of *N.+_cells* and *Counts* (using the larger of the 2 reported values). The results for genes with
#'larger *N.+_cells* and *Counts* are more accurate, yet the additional accuracy benefit rapidly decreases at *N.+_cells > 30-50* and *Counts > 100*.  
#
#'To include more lower expression genes into the analysis, I recommend running DGE.2samples() twice, once with *min.cells=3*, *min.count=0*, and default *icc="i"* and once with   
#' *min.cells=3*, *min.count=0*, and *icc="1"*, and follow the logic of the example below:
#
DGE.results.mc3=DGE.2samples(Neuron.QC, ident.1 = "glutcells", ident.2 = "GABAcells", min.cells=3, min.count = 0)
DGE.results.mc3.icc1=DGE.2samples(Neuron.QC, ident.1 = "glutcells", ident.2 = "GABAcells", min.cells=3, min.count = 0, icc=1)
# The multiple R warnings are for the chi squared test, which cannot be used for genes with low counts anyway. 
#
# Compared to the default parameters, this analysis includes over 4,000 additional lower expression genes
low.exp.genes=setdiff(rownames(DGE.results.mc3),rownames(DGE.results))
# The results just for these genes are:
DGE.results.mc3.LE=DGE.results.mc3[low.exp.genes,]
DGE.results.mc3.icc1.LE=DGE.results.mc3.icc1[low.exp.genes,]
# Some of these results may be unreliable because of reduced accuracy of cell weights calculations based on iterative icc (icc="i") setting. The results for cells with
# similar weights are expected to be more accurate. To select the corresponding genes, we keep only the results with Log2FC within ~30% at icc="i" and icc=1 (equal weights)
DGE.results.LE=DGE.results.mc3.icc1.LE[(abs(DGE.results.mc3.icc1.LE$log2FC*DGE.results.mc3.LE$log2FC)>0&
                                              (abs(DGE.results.mc3.icc1.LE$log2FC-DGE.results.mc3.LE$log2FC)<0.37|
                                                 (abs(DGE.results.mc3.icc1.LE$log2FC)==Inf&abs(DGE.results.mc3.LE$log2FC)==Inf))),]
# Here we selected the results from the DGE.results.mc3.icc1.LE dataframe to err on the side of caution (using more conservative statistical test with icc=1). 
# Among the remaining 2,819 low expression genes
DGE.results.LE.p01=DGE.results.LE[DGE.results.LE$p.value<0.01,]
# 303 have p<0.01. All of them are upregulated/dowregulated more than 2-fold and
DGE.results.LE.p01.log2FC4=DGE.results.LE.p01[abs(DGE.results.LE.p01$log2FC)>4,]
# 194 are up/downregulated more than 16-fold. And,
DGE.results.LE.p0001.log2FC4=DGE.results.LE.p01.log2FC4[DGE.results.LE.p01.log2FC4$p.value<0.0001,]
# 20 of the latter have p<0.0001.
#
# The default analysis missed all these genes because of more stringent accuracy requirements. The DGE in these genes may be real or may be an artifact of insufficient 
# transcript sampling (too few cells and too few counts). Depending on questions being addressed by experiments, one may want to at least know about such genes and 
# select the ones for validation by additional experiments (using not only on log2FC and p-values but also on N.+_cells and Counts values). 

##################################################################################
#'*Section 2.2 Optimizing multiple sample DGE analysis*
#
# While the idea behind optimizing multiple sample analysis parameters is similar, there are several important distinctions. 
#
#'*2.2.1 Utilizing and optimizing min.pct setting*
#
#'The default *min.pct=0.03* setting ensures that at least 3% of all cells across all samples in one of the two groups express the gene. This is different from requiring that
#'that at least 3 samples in one of the two sample groups have *> min.cells* cells expressing the gene. The *min.pct* setting is therefore meaningful and should be used. The
#'specific value for this setting depends on the goals and specific details of the study. A general recommendation may be counterproductive. 
#
#'*2.2.2 Utilizing and optimizing min.count setting*
#'
#'The default *min.count=10* setting has no effect on analysis at *min.cells > 3* because 3 x *min.cells > 3* > 10. At default *min.cells=10*, at least one of the two sample 
#'groups will have at least 3 x 10 = 30 counts. This parameter is provided mostly for compatibility with DGE.2samples() function, enabling one to compare multiple sample and 
#'2-sample analysis of the same data with the same parameters (An exercise I recommend to see pitfalls of 2-sample analysis when data are available from multiple samples).
#'
#'*2.2.3 Utilizing and optimizing t.test setting*
#'
#'The default *t.test=FALSE* setting is optimal for many (at least5+) samples per group when each sample contains many (>10-30) cells expressing genes of interest. The other 
#'*t.test=TRUE* setting is optimal when the number of samples per group is small (3-4) and all samples have about the same number of cells. The reality is usually between these 
#'two limits. Then, it is worthwhile to perform analysis with both settings. Similar results for both settings are most reliable. Different results indicate that one of the tests
#'is not accurate. In that case, I recommend erring on the side of caution and choosing the more conservative results, the results of which need to be validated anyway. 
#'
#'*2.2.4 Identifying DGE in low expression genes at non-default icc and min.cells settings*
#'
#'The defaullt *icc="i"* setting is optimal for genes with moderate to high expression (Counts/cell>0.1) detected in many cells (>10-30 cells per sample). The *icc=1* setting is not 
#'appropriate for multiple sample analysis. The "icc=0" setting is optimal for very low expression genes (Counts/cell<0.05-0.1). With this in mind I recommend running CntAv() once
#'with default icc="i" and once with icc=0 followed by the analysis illustrated below.
#
# To illustrate analysis of low expression genes, let's assemble a new mock multiple sample experiment, in which we modify 100 randomly selected with Counts/cell < 0.3 and total counts per 
# sample group > 20
GLUT.counts.under01=GLUT.counts[(rowSums(GLUT.counts)/ncol(GLUT.counts)<0.3&rowSums(GLUT.counts)>20),]
genes=rownames(GLUT.counts.under01)
set.seed(5)
modgenes=sample(genes,100)
Tr.r1.mat=as.matrix(Untr.r1[["RNA"]]$counts)
Tr.r1.mat[modgenes,]=round(2*Tr.r1.mat[modgenes,])
Tr.r1=CreateSeuratObject(Tr.r1.mat)
Idents(Tr.r1)="treated.1"
Tr.r2.mat=as.matrix(Untr.r2[["RNA"]]$counts)
Tr.r2.mat[modgenes,]=round(2*Tr.r2.mat[modgenes,])
Tr.r2=CreateSeuratObject(Tr.r2.mat)
Idents(Tr.r2)="treated.2"
Tr.r3.mat=as.matrix(Untr.r3[["RNA"]]$counts)
Tr.r3.mat[modgenes,]=round(2*Tr.r3.mat[modgenes,])
Tr.r3=CreateSeuratObject(Tr.r3.mat)
Idents(Tr.r3)="treated.3"
GLUT.MS=JoinLayers(merge(Untr.r1, y=c(Untr.r2,Untr.r3,Tr.r1,Tr.r2,Tr.r3)))
#
# Default CntAv()
GLUT.MS.i=CntAv(GLUT.MS)
GLUT.SampleMatrix=SampleMatrix(GLUT.MS.i,samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"))
DGE.multi.results=WT.MultiSample(GLUT.SampleMatrix)
DGE.ms.results=DGE.multi.results$DGE
DGE.ms.p05=DGE.ms.results[DGE.ms.results$p.value<0.05,]
# finds just 43 of the 100 modified genes, partly because 22 of the 100 modified genes have been excluded by default minimum expression requirements.
length(intersect(modgenes,rownames(DGE.ms.results)))
# We can now reanalyze to the data for low expression genes 
GLUT.MS.0=CntAv(GLUT.MS, icc = 0)
GLUT.SampleMatrix.0=SampleMatrix(GLUT.MS.0,samples.1 = c("untreated.1","untreated.2","untreated.3"), samples.2 = c("treated.1","treated.2","treated.3"))
DGE.multi.results.0=WT.MultiSample(GLUT.SampleMatrix, t.test=TRUE, min.pct=0, min.cells=3, min.count = 0)
DGE.ms.results.0=DGE.multi.results.0$DGE
# where we pushed all parameters to their limit.
# Next, we merge the default analysis results at Counts/cell>=0.1 with the second analysis results at Counts/cell<0.1
DGE.ms.results.filt=DGE.ms.results[(DGE.ms.results$`Counts/cell.1`>=0.1|DGE.ms.results$`Counts/cell.2`>=0.1),]
DGE.ms.results.filt.0=DGE.ms.results.0[(DGE.ms.results.0$`Counts/cell.1`<0.1&DGE.ms.results.0$`Counts/cell.2`<0.1),]
DGE.ms.results.comb=rbind(DGE.ms.results.filt,DGE.ms.results.filt.0)
length(intersect(modgenes,rownames(DGE.ms.results)))
# In this analysis, we find 67 of the 100 genes we modified with p<0.05
DGE.ms.results.comb.p05=DGE.ms.results.comb[DGE.ms.results.comb$p.value<0.05,]
# and 82 of genes we modifies with borderline p<0.1
DGE.ms.results.comb.p01=DGE.ms.results.comb[DGE.ms.results.comb$p.value<0.1,]
# This analysis excluded only 3 of the genes we modifies based on minimum expression requirements.
length(intersect(modgenes,rownames(DGE.ms.results.comb)))
# Overall, the results are about as good as they get for analysis of low expression genes with Counts/cell<0.3
#
# Importantly, this example illustrates only one out of many possible ways to combine analysis of higher expression genes based on weighted averaging (icc="i"), and weighted
# t.test (t.test=FALSE) with analysis of lower expression genes based on aggregate counts (icc=0) and non-weighted t.test (t.test=TRUE). The threshold for transitioning
# between the two approaches as well as the parameters within each of the approaches may be modified. For instance, our mock experiment contains only 3 samples per group with fairly
# even number of cells per sample. In this case the t.test=TRUE setting yet icc="i" actually makes more sense for higher expression genes. If we do that 
DGE.multi.results.t=WT.MultiSample(GLUT.SampleMatrix, t.test=TRUE)
DGE.ms.results.t=DGE.multi.results.t$DGE
DGE.ms.results.filt.t=DGE.ms.results.t[(DGE.ms.results.t$`Counts/cell.1`>=0.1|DGE.ms.results.t$`Counts/cell.2`>=0.1),]
DGE.ms.results.comb.t=rbind(DGE.ms.results.filt.t,DGE.ms.results.filt.0)
# we recover 71 of the 100 modified genes with p<0.05 
DGE.ms.results.comb.t.p05=DGE.ms.results.comb.t[DGE.ms.results.comb.t$p.value<0.05,]
# and 84 of the 100 modified genes with p<0.1
DGE.ms.results.comb.t.p01=DGE.ms.results.comb.t[DGE.ms.results.comb.t$p.value<0.1,]
# which is a better result because the more accurate statistical test has been used. 
#
# I couldn't come up with a single recipe for the best approach so far. I think that the approach should really depend on the study goals and methods for validating the results for
# low expression genes, for which accurate statistics is impossible due to low cell and count numbers. The tricks illustrated above just help to expand the list of candidate genes 
# using less accurate yet still unbiased statistical analysis.
#
#######################################################################################
