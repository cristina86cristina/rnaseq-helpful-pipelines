# rnaseq-helpful-pipelines

This is the pipeline I am using for RNAseq analysis for human samples. Mapping is done with Kallisto. 
The following scripts will help with:
- summarise results at the gene level
- calculate cpm (old)
- calculate total number of reads mapping to protein coding genes
- calculate relationship between cpm and gene size
- extract tpm 
- write a file/sample with tpm info
- write a file/sample with a few other metrics 

Second script is for the QC - this is more interactive, people can change as needed (and use samples needed for the analysis)
- transform tpm <= 0.001 
- log2 tpm
- density plot 
- pca plots
