**R code for analysis of single cell and spatial transcriptomics data**

sc&spRNASeqFunctions.R file contains the following simple Seurat-compatible functions for differential gene expression (DGE) expression analysis:
1. **DGE.2samples()** performs analysis of DGE between 2 identities in a Seurat object using weighted averaging and wieghted-t-test/chi-squared-test combination.
2. **IterWghtTest()** same as DGE.2samples but performs only weighted averaging and weighted t-test (called by DGE.2Samples).
3. **Chi2Test()** compares 2 identities based on chi-squared test for aggregated counts (called by DGE.2samples).
4. **DGE.Multisample()** performs analysis of DGE between 2 groups of samples (at least 3 samples in each group), uses a single integrated Seurat object with at least 6 identities, performs weighted averaging and weighted t-test.

#####################################################

Both DGE.2samples() and DGE.Multisample() can be used for scRNASeq and spRNASeq. However, different parameters are recommended for different assays (see function documentation).

#####################################################

All functions in the file must be loaded for the all 4 basic functions listed above to work.
The Seurat object must contain an active "RNA" assay with a single counts layer containing all data to be analayzed (see Preparing Seurat Objects).

######################################################

**Additional functions included in the sc&spRNASeqFunctions.R file and called by the 4 basic functions above:**

alt.wttest() performs weighted t-test using optimized variance estimator (resolves known weights::wtd.t.test() issues).

alt.wttest2() same as altwttest but with an additional correction for effective degrees of freedom.

ICC.AN() calculates ANOVA intracluster correlation coefficient (ICC).

ICC.iter() calculates iterative ICC providing a more accurate match bewteen estimated and expected variances of the data compared to ANOVA ICC.

ICCWeight() calculates statistical weights of cells based on ICC values.

CntAv() performs wiehgted averaging or aggregation of counts for multisample DGE analysis, depending on parameter setting.

SampleMatrix() assembles a matrix weighted average counts and variances for multisample DGE analysis.

IterVar() calculates iterative ICC for multisample analysis.

WT.MultiSample() performs weighted t-test analysis of DGE using SampleMatrix() output.
 
