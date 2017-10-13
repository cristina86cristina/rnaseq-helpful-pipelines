# rnaseq-helpful-pipelines

This is the pipeline I am using for RNAseq analysis for human samples. Mapping is done with Kallisto. 

##processing_rnaseq.r 
The following scripts will help with:
- summarise results at the gene level
- calculate cpm (old)
- calculate total number of reads mapping to protein coding genes
- calculate relationship between cpm and gene size
- extract tpm 
- write a file/sample with tpm info
- write a file/sample with a few other metrics 
Inside the folder with kallisto results, from the terminal: R CMD BATCH processing_rnaseq.r






##annotation.r


This script is to annotate several files where the Id are Ensebl genes and we want gene names, biotype and gene description. Initially written for SARTools output. 
