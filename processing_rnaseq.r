#First part RNAseq 
#- read kallisto output
#- gene level data 
#- calculate cpm (old) - we are now extracting tpm
#- write files for each sample
#- summary of read counts


library(base)
library(BiocParallel)
library(RCurl)
library(tximport)
library(readr)
library(biomaRt)



accessions <- list.dirs(full.names=FALSE)[-1]

mart <- biomaRt::useMart(biomart = "ensembl", dataset =  "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "description", "transcript_biotype"), mart = mart)
t2g$target_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep=".") # append version number to the transcript ID
t2g[,c("ensembl_transcript_id","transcript_version")] <- list(NULL) # delete the ensembl transcript ID and transcript version columns
t2g <- dplyr::rename( t2g, gene_symbol = external_gene_name, full_name = description, biotype = transcript_biotype )
t2g<-t2g[,c(ncol(t2g),1:(ncol(t2g)-1))]
#Let's use tximport to summarize results into genes
kallisto.dir<-paste0(accessions)
kallisto.files<-file.path(kallisto.dir,"abundance.tsv")
names(kallisto.files)<- accessions
tx.kallisto <- tximport(kallisto.files, type = "kallisto", tx2gene = t2g, 
                        reader = read_tsv, countsFromAbundance ="no")

#Select only protein coding genes



#load packages
library(dplyr)
library(sm)
#install biomart
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#library(biomaRt)
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gb <- getBM(attributes=c("ensembl_gene_id","gene_biotype"), mart=mart)
detach("package:biomaRt", unload=TRUE) #unload biomart because it creates problems with dplyr
gb_coding<-subset(gb, gb$gene_biotype=="protein_coding")
genes<-gb_coding$ensembl_gene_id
counts<-as.data.frame(tx.kallisto$counts[row.names(tx.kallisto$counts) %in% genes, ])
tpm <- as.data.frame(tx.kallisto$abundance[row.names(tx.kallisto$abundance) %in% genes, ])
len <- as.data.frame(tx.kallisto$len[row.names(tx.kallisto$len) %in% genes, ])





#Let's divide the count for the total read counts - we then split the count file and write a new file for each sample


ids<-rownames(counts)

total_counts<-apply(counts,2,sum)
counts_divided<-sweep(counts, 2, total_counts, `/`)
cpm<-counts_divided*1000000
write.csv(total_counts,"total_counts_mapped.csv")

#let's check gene expression - gene size relationship 

cpm[cpm <= 0.001] <- 0.001

corre <- 0
for (i in 1:length(cpm)){
  corre[i]<-cor(log2(cpm[i]),log2(len[i]))
  print(corre[i])
}

corre_df<-data.frame(sample_id=names(cpm),cor_value=corre)
write.csv(corre_df,"correlation_matrix.csv")


lapply(names(tpm), function(x) write.table(tpm[[x]],paste(x,"tpm.txt",sep=""),row.names=ids,quote=F,col.names=F,sep="\t"))


