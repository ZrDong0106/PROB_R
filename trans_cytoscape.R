trans_cytoscape=function(adj_mat,gene_names){
  a1=1;b1=1;c1=1;
  docu=data.frame(a1,b1,c1);
  colnames(docu)=c("gene1","gene2","pos/neg");
  flg=1;ng=length(gene_names);
  for(i in 1:ng){
    for(j in 1:ng)
    {
      if(adj_mat[i,j]!=0) {
        pos=1;if(adj_mat[i,j]<0) {pos=-1;}
        docu[flg,1]=gene_names[j];
        docu[flg,2]=gene_names[i];
        docu[flg,3]=pos;
        flg=flg+1;
      }
    }
  }
  return (docu);
}