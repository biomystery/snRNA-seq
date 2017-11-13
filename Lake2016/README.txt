##########################################################################################
### Clustering and Classification.R
###
### Version 1.1
### Last updated, Mar 2016
### Rizi Ai
### riai@ucsd.edu
### University of California San Diego
##########################################################################################
Required packages:

R version >= 3.1
library(randomForest)
library(gplots)
##########################################################################################
Input:

gene_matrix.dat  # input gene expression matrix (rows: genes; column: samples)
expression_type  # "1" (default) for log2 expression (eg. log2TPM, log2FPKM, log2Counts); 
                 # "2" for expression without taking log2 transformation
##########################################################################################
Output:

log  # log file
01.Report.txt                      # Summary and statisitics 
02.Expression-variation_stat.txt   # Gene expression variation (eg. CV2)
03.Fitting_stat.pdf                # Fitting over-dispersed genes to an inverse distribution
04.Over-dispersed_gene-list.txt    # Over-dispersed genes
05.Selected-features.txt           # Genes selected during feature selection
06.Classification_stat.pdf         # Plots of mis-classification error rates and number 
                                   # of samples in each cluster during cross-validation
07.Cluster1_sample-list.txt        # Samples associated with cluster 1
07.Cluster2_sample-list.txt        # Samples associated with cluster 2
08.All_genes_stat.txt              # Statistics of all genes
09.Cluster1_AllGenes_NewInput.dat  # Gene expression matrix with all genes in cluster 1
                                   # This matrix can be used for the next round of 
                                   # Clustering and Classification analysis
09.Cluster2_AllGenes_NewInput.dat  # Gene expression matrix with all genes in cluster 2
                                   # This matrix can be used for the next round of 
                                   # Clustering and Classification analysis
10.Heatmap.pdf                     # Heatmap of the clusters using 2-fold DEGs
##########################################################################################


