Progression_Inference=function(DATA){
  #-----Progression_Inference-------
  ngene=nrow(DATA)-1;
  npatient=ncol(DATA);
  data=DATA[1:ngene,];
  grade=DATA[ngene+1,];
  
  #-----Cauculate Gaussian similarity function-----
  omega=matrix(rep(0,npatient*npatient),ncol=npatient,nrow=npatient);
  gamma=omega;S=omega;Ks=omega;
  for (i in 1:npatient)
    for(j in 1:npatient){
      omega[i,j]=(1+abs(grade[1,i]-grade[1,j]));
    }
  
  k=10; #number of neighbors considered in knn
  library(RANN)
  knnres=nn2(t(data),t(data),k+1); 
  idx_k=knnres$nn.idx[,k+1];
  dists=knnres$nn.dists;
  sigma=dists[,k+1];
  sigma2=sigma^2;
  for (i in 1:npatient)
    for(j in 1:npatient){
      gamma[i,j]=(sigma2[i]+sigma2[j])/omega[i,j];
    }
  
  for (i in 1:npatient)
    for(j in 1:npatient){
      S[i,j]=exp(-(sum((data[,i]-data[,j])^2))/gamma[i,j]); 
      #gaussian similarity function
      
    }
  
  #-----Cauculate transition matrix-----
  D=rep(0,npatient);
  H=S;E=D;P=S;
  for(i in 1:npatient)
    D[i]=sum(S[i,]);
  for(i in 1:npatient)
    for(j in 1:npatient)
      H[i,j]=S[i,j]/D[i]/D[j];
  for(i in 1:npatient)
  {H[i,i]=0;E[i]=sum(H[i,]);}
  for(i in 1:npatient)
    for(j in 1:npatient)
      P[i,j]=E[i]^(-0.5)*H[i,j]*E[j]^(-0.5); # transition matrix
  
  phi0=E/sqrt(sum(E^2)); # largest eigenvector of P
  phi0=as.matrix(phi0);
  library(MASS)
  Q=ginv((diag(npatient)-P+phi0%*%t(phi0)))-diag(npatient); # accumulated transition matrix
  
  #---Select root patient
  grade=as.numeric(grade);
  Ind_max=subset(1:npatient,grade==max(grade));
  b=length(Ind_max);
  #b=1;
  rn=ceiling(runif(1)*b);
  x_ref=Ind_max[rn];
  
  depth_to_root=function(Q,s){
    nQ=nrow(Q);
    dpt=rep(0,nQ);
    for(i in 1:nQ){
      dpt[i]=sqrt(sum((Q[i,]-Q[s,])^2));
    }
    return(dpt);
  }
  
  TPD_to_xref=depth_to_root(Q,x_ref);
  Ind_min=subset(1:npatient,grade==min(grade));
  Ind_sort=order(TPD_to_xref,decreasing=TRUE);
  
  for(i in 1:length(Ind_sort)){
    if(grade[Ind_sort[i]]==1) {
      root=Ind_sort[i];break;
    }
    #select the root with maximum TPD to x_ref in the minimum grade subgroup
  }
  #-----Transfer into a smoothed trajectory----
  #root=58;
  PPD=depth_to_root(Q,root);
  PPD_order=order(PPD);
  smoothL=10^(floor(log10(npatient))-1);
  
 # smthdata=data[1,PPD_order];
  #g_smooth=ksmooth(PPD_order,smthdata,kernel="normal",bandwidth=smoothL);
  
  
  point_selected=((smoothL/2+1):(npatient-smoothL/2))
  
  data_ordered=data;
  for(i in 1:ngene){
    smthdata=data[i,PPD_order];
    g_smooth=ksmooth(1:npatient,smthdata,kernel="normal",bandwidth=smoothL);
    data_ordered[i,]=g_smooth$y;
  }
  data_ordered[ngene+1,]=grade[PPD_order];
  data_ordered=data_ordered[,point_selected];
  time_sampled=point_selected/npatient;
  
  return(list(Accumulated_Transition_Matrix=Q,Temporal_Progression=PPD,Ordered_Data=data_ordered,Sampled_Time=time_sampled,Order=PPD_order));
  
}
