################################################################################################
####### file_name: annotation.r
####### file_description: this script annotates several files with a column of Ensembl gene ids into external gene name + description
###### Written with Sartools results in mind
##################################################################################################
#######InPUT FILE
###### Id stuff1 stuff2 
###### ENSGxx 2 3
##################################################################################################



#This function annotates Sartools results - 
library(biomaRt)
annotation_meow = function(Species,datafile,fileout){  ##meow meow
  
  #Read of the output of SARTools
  data2=read.delim(toString(datafile), header =T, sep ="\t") ##change csv if you have that
  
  #annotation
  bmdataset = toString(Species)
  mart=useMart(biomart="ensembl", dataset= bmdataset)
  ann <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), filters="ensembl_gene_id",values=data2$Id,mart = mart)
  data2$ensembl_gene_id<-data2$Id
  data2$Id<-NULL
  data_ann<-merge(data2,ann,by="ensembl_gene_id")
  fileout<- paste("Annotated",datafile,sep="")
  write.csv(data_ann,fileout,row.names=FALSE)
}

fileNames <- Sys.glob("*.txt")
for (fileName in fileNames) {
  annotation_meow("hsapiens_gene_ensembl",fileName) #change here if you have another species
  
}

