Locate_Key_Genes=function(BLin,tcgname){
  BLAMs=abs(BLin$Ajacent_Matrix)/BLin$Standard_Deviations;
  BtB=BLAMs%*%t(BLAMs);
  BtB=BLAMs;
  
  principle_vec=eigen(BtB)$vectors[,1];
  names(principle_vec)=tcgname;
  principle_vec=abs(principle_vec);
  eig_scores=sort(principle_vec,decreasing=TRUE);
  
  return(eig_scores);
}