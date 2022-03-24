KM_analysis=function(survtime,event,gene_name){
  gene_m=as.numeric(Gene_GSE7390[gene_name,]);
  med_m=median(gene_m);
  group_m=rep(1,length(survtime));
  group_m[gene_m>med_m]=2;
  library(survival)
  Lusurv=Surv(survtime,event);
  Lufit=survfit(Lusurv~group_m);
  sdf=survdiff(Lusurv~group_m);
  pval=1-pchisq(sdf$chisq,1);
  plot(Lufit,mark.time=T,conf.int='none',col=c("red","blue"),lwd=2,xlab="Time(Years)",ylab="Survival probability",main=gene_name,cex.axis=1.5,cex.lab=1.5,cex.main=1.5);
  #legend("bottomleft",c(paste(c("Low",gene_name),collapse=" "),paste(c("High",gene_name),collapse=" "))
     #    ,col=c("red","blue"),lwd=2,bty='n',cex=1.5);
  #sdf
  text(4,0.2,paste(c("p=",signif(pval,2)),collapse = " "),cex=1.5);
}
