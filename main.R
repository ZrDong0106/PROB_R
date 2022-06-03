# Package names
packages=c("trend", "pheatmap", "OmnipathR", "tidyr", "gprofiler2", "minerva", "reshape2", "ggplot2"
           ,"Biobase","GEOquery","monomvn","igraph","survival","Brq");

# Install packages not yet installed
installed_packages=packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


source("BL_to_csv.R");
source("KM_analysis.R");
source("Locate_Key_Genes.R");
source("ODE_Bayesian_Lasso.R");
source("PROB_GEOinstall.R");
source("Progression_Inference.R");
source("Time_course.R");
source("trans_cytoscape.R");

#---Read GSE7390
Gene_Data=read.csv("Gene_GSE7390.csv");
row.names(Gene_Data)=Gene_Data[,1];
Gene_Data=Gene_Data[,-1];
Gene_Data=data.frame(Gene_Data);
#---Progression Inference-------------

PI=Progression_Inference(Gene_Data);
save(PI,file="GSE7390_PI.Rdata");

pseudo_series=PI$Ordered_Data;
pseudo_time=PI$Sampled_Time;
#plot(pseudo_time,pseudo_series[1,],type="l");

#---Identify TCGs-------------------------
library(trend)
ngenes=nrow(pseudo_series)-1;
TCG_mark=rep(FALSE,ngenes);
pval_trend=rep(0,ngenes);
for(i in 1:ngenes){
  trial=as.numeric(pseudo_series[i,]);
  res.t=mk.test(trial); #use Mann-Kendall test to identify temporal trend of each gene
  if(res.t$pvalg<0.05) TCG_mark[i]=TRUE;
  pval_trend[i]=res.t$pvalg;
}

TCGid=order(pval_trend,decreasing=FALSE)[1:100];
save(TCGid,file="TCGid.Rdata");


grade2=pseudo_series[ngenes+1,];
TCG_series=rbind(pseudo_series[TCGid,],grade2);
#plot(pseudo_time,TCG_series["FOXM1",],type="l");
save(TCG_series,file="TCG_series.Rdata");
TCG_names=row.names(pseudo_series)[TCGid];

#---Draw heatmap
library(pheatmap)
dev.off()
Data=TCG_series[-nrow(TCG_series),];
M=apply(Data,1,max)
A=apply(Data,1,mean)
S=apply(Data,1,sd)
D=(Data-A)/S
dev.new()
pheatmap(D,cluster_row=T, cluster_cols=F, clustering_distance_rows='euclidean',clustering_method = "ward.D", color = colorRampPalette(c("CornflowerBlue", "white", "firebrick3"))(8), fontsize=9, fontsize_row=6,labRow=NA, show_colnames = FALSE,show_rownames=FALSE)

#----Specify prior network
library(OmnipathR)
library(tidyr)
library(gprofiler2)
iai_all=import_all_interactions();
w1_all=which((iai_all$source_genesymbol %in% TCG_names) & (iai_all$target_genesymbol %in% TCG_names));
ntcg=length(TCG_names);
prior_PPI =matrix(rep(0,ntcg*ntcg),nrow=ntcg,ncol=ntcg);
for(i in w1_all){
  sel=iai_all[i,];
  prior_PPI[which(TCG_names==sel$target_genesymbol),which(TCG_names==sel$source_genesymbol)]=1;
}


library(minerva);
Gene_TCGs=t(Gene_Data[TCGid,]);
MI=mine(Gene_TCGs)$MIC;
MI95=quantile(MI,0.95);
prior=(MI>MI95);
for(i in 1:100) prior[i,i]=FALSE;


#----Bayesian lasso
set.seed(1);
Breast_BL=ODE_Bayesian_Lasso(TCG_series,pseudo_time,prior|prior_PPI);
#save(Breast_BL_MIomnipath_all_ver3,file="Breast_BL_MIomnipath_all_ver3_MI90.Rdata");
BL_to_csv(Breast_BL,0.95,"GSE7390_GRNs.csv",TCG_names);

#--Draw bubble plot
library(reshape2)
library(ggplot2)
BL3=Breast_BL;
am=BL3$Ajacent_Matrix;
rownames(am)=colnames(am)=TCG_names;
sd=BL3$Standard_Deviations;
row_id=(apply(abs(am),1,sum)>1e-5);
col_id=(apply(abs(am),2,sum)>1e-5);
am=am[row_id,col_id];
sd=sd[row_id,col_id];
data_melt=melt(am);
names(data_melt)=c('Gene1','Gene2','Value');
p=ggplot(data_melt,aes(x=Gene2,y=Gene1,size=abs(am)/sd,color=am))+geom_point()+theme(axis.text.x = element_text(angle=90,hjust=1))+
  ylab("Target genes")+xlab("Genes")+scale_colour_gradient(low="green",high="red");
p


#----Locate key genes
Eig_scores=Locate_Key_Genes(Breast_BL,TCG_names);
cut=5;
library(ggplot2)
trt=names(Eig_scores)[1:cut];
outcome=Eig_scores[1:cut];
df=data.frame(trt,outcome)
p=ggplot(df, aes(reorder(trt,-outcome), outcome)) +
  geom_bar(aes(fill=outcome),stat="identity")+xlab("")+ylab("Hub scores")+
  scale_fill_gradient(low = "#00FFFF", high = "#0000FF", na.value = NA)+
  theme_minimal(base_size=20)#+theme(axis.text.x = element_text(angle=90, hjust=1));
p
Time_course(Eig_scores,cut=5,TCG_series,pseudo_time);

#----Prepare survival time data
library(GEOquery)
Sys.setenv(VROOM_CONNECTION_SIZE=1e8);
my_id="GSE7390";
gset=getGEO(my_id,GSEMatrix =TRUE, getGPL=FALSE);
FUN7=function(strg){return(as.numeric(substring(strg,first=7)))};
FUN8=function(strg){return(as.numeric(substring(strg,first=8)))};
FUN9=function(strg){return(as.numeric(substring(strg,first=9)))};

os=as.array(gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.15"]]);
os=apply(os,1,FUN=FUN7)/365;
eos=as.array(gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.16"]]);
eos=apply(eos,1,FUN7);
rfs=as.array(gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.13"]]);
rfs=apply(rfs,1,FUN=FUN8)/365;
erfs=as.array(gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.14"]]);
erfs=apply(erfs,1,FUN8);
dmfs=as.array(gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.17"]]);
dmfs=apply(dmfs,1,FUN=FUN9)/365;
edmfs=as.array(gset[["GSE7390_series_matrix.txt.gz"]]@phenoData@data[["characteristics_ch1.18"]]);
edmfs=apply(edmfs,1,FUN=FUN9);
#os,rfs,dmfs: survival time data for samples with respect to overall survival, relapse free survival and distant metastasis-free survival
#eos,erfs,edmfs: 0/1 vector for samples identifying whether corresponding events are observed

#---KM analysis
cut=1;
top_id=names(Eig_scores);
par(mfcol=(c(cut,3)))
for(i in 1:cut){
  ng=top_id[i]; #ng is the gene that you want to apply KM analysis.
  KM_analysis(os,eos,ng);text(4,0.4,"OS",cex=1.5);
  KM_analysis(rfs,erfs,ng);text(4,0.4,"RFS",cex=1.5);
  KM_analysis(dmfs,edmfs,ng);text(4,0.4,"DMFS",cex=1.5);
  }

