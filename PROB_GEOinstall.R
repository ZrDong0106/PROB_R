#install.packages("BiocManager")
#install.packages("forcats")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readr")
#install.packages("tidyr")
#install.packages("survminer")
#BiocManager::install("GEOquery")
PROB_GEOinstall=function(){
  Sys.setenv(VROOM_CONNECTION_SIZE=1e8);
library(GEOquery)
library(Biobase)
my_id="GSE7390";
gset=getGEO(my_id,GSEMatrix =TRUE, getGPL=FALSE);

data2=gset[["GSE7390_series_matrix.txt.gz"]]@assayData[["exprs"]];


#------Download GPL
options('download.file.method.GEOquery'='auto');
options('GEOquery.inmemory.gpl'=FALSE);
gpl96=getGEO("GPL96");
probe2symbol=Table(gpl96)[,c(1,11)];
probe96=Table(gpl96)[,1];
Gene_symbol96=Table(gpl96)[,11];
row.names(probe2symbol)=probe2symbol[,1];
#row.names(data2)=probe2symbol[row.names(data2),2];

#----Mark gene symbols
data2=as.data.frame(data2);
ID=row.names(data2);
gene_symbol=probe2symbol[ID,2];
data2=cbind(data2,ID,gene_symbol);
test2=aggregate(data2[,1:198],by=list(data2$gene_symbol),FUN=mean);

test2=test2[-1,];
row.names(test2)=test2[,1];
test2=test2[,-1];

grade=gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["grade:ch1"]];
grade[120]=grade[127]="1";
grade=as.numeric(grade);
test2=rbind(test2,grade);
write.csv(test2,"Gene_GSE7390.csv");
}