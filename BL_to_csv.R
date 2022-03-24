BL_to_csv=function(BLin,threshold,file_name,gene_names){
  docu=trans_cytoscape(BLin$Ajacent_Matrix*(BLin$Presence_Probability>threshold),gene_names);
  write.csv(docu,file=file_name);
}