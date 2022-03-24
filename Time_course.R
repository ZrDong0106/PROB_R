Time_course=function(eig_scores1,cut,TCG_series1,pseudo_time1){
#---Time course of top genes
library(RColorBrewer)
top_id=names(eig_scores1)[1:cut];
#top_id=c("KIF2C","UBE2C","MCM10","MELK");

#row.names(TCG_series)=c(TCG_names,"grade");
top_series=TCG_series1[top_id,];
top_series=t(apply(top_series,1,scale)); #standardize
cols=brewer.pal(cut,"Dark2");
plot(pseudo_time1,top_series[1,],ylim=c(min(top_series)-1,max(top_series)+1),type="l",xlab="Pseudotime",ylab="Standardized gene expression",col=cols[1],lwd=2,cex.lab=1.5,cex.axis=1.5);
for(i in 2:cut){
  lines(pseudo_time1,top_series[i,],col=cols[i],lwd=2);
}
legend(x='bottomright',legend=top_id,lty=1,col=cols,bty='n',lwd=2,cex=1);
}