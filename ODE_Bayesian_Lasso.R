ODE_Bayesian_Lasso=function(Data_ordered,Time_Sampled,Prior_Network,penalty=1,max_regulator=NULL){
  ngene=nrow(Data_ordered)-1;
  npatient=ncol(Data_ordered);
  #-----Bayesian Lasso ODE---
  
  library(monomvn)
  library(Brq)
  d=Data_ordered[1:ngene,];
  t=Time_Sampled;
  dsize=ncol(d);
  AM_adjusted=diag(ngene);
  AM=diag(ngene);
  S=diag(ngene);
  Std=matrix(rep(1,ngene*ngene),nrow=ngene,ncol=ngene);
  for (i in 1:ngene){
    ngene_i=sum(Prior_Network[i,]);
    if(ngene_i==0) {AM[i,i]=S[i,i]=AM_adjusted[i,i]=0;}
    else{
      #--numerical differentiation---
      y_output=(d[i,2:dsize]-d[i,1:(dsize-1)])/(t[2:dsize]-t[1:(dsize-1)]);
      #s21=var(y_output-mean(y_output));
      #if(s21==0) {continue;}
      #--reform data---
      x_input=as.data.frame(d[i,]);
      for(j in 1:ngene){
        if(Prior_Network[i,j]==1)
          x_input=rbind(x_input,d[i,]*d[j,]);
      }
      ngene_i=sum(Prior_Network[i,]);
      
      x_input=t(x_input);
      x_input=x_input[1:(dsize-1),];
      x_input=cbind(x_input[,-1],x_input[,1]);
      colnames(x_input)=c(1:ngene_i,"X");
      y_output=t(y_output);
      
      
      if(is.null(max_regulator)) 
      {fit=blasso(x_input,y_output,icept=FALSE,normalize=FALSE,lambda2=penalty);}
      else if(max_regulator>ncol(x_input))
        {fit=blasso(x_input,y_output,icept=FALSE,lambda2=penalty,M=min(max_regulator,ncol(x_input)));}
      else if(max_regulator<=ncol(x_input))
        {fit=blasso(x_input,y_output,icept=FALSE,M=min(max_regulator,ncol(x_input)));}
      beta=fit$beta;
      for(j in 1:ngene){
        if(i==j | Prior_Network[i,j]==0){S[i,j]=0;AM[i,j]=0;}
        else{
          #densityj=density(beta[,j]);
          #betaj=densityj$x[densityj$y==max(densityj$y)];
          j0=sum(Prior_Network[i,1:j]);
          betaj=mean(beta[,j0]);
          AM[i,j]=betaj[1]; #AM stores the posterior mean of each para
          #---Cauculate presence probability of edges in the network
          S[i,j]=0;
          alpmax=0;
          stdj=sqrt(var(beta[,j0]));
          for(k in 1:99){
            alphak=0.01*k;
            #lbk=quantile(beta[,j],alpha/2,na.rm=TRUE);
            #rbk=quantile(beta[,j],1-alpha/2,na.rm=TRUE);
            
            lbk=AM[i,j]+stdj*qt(alphak/2,1e3-1);
            rbk=AM[i,j]-stdj*qt(alphak/2,1e3-1);
            if(lbk<0 & rbk>0){ alpmax=alphak;}
            #else {AM_adjusted[i,j]=AM[i,j];}
          }
          S[i,j]=1-alpmax;
          AM_adjusted[i,j]=S[i,j]*(S[i,j]>0.95);  
          Std[i,j]=stdj;
        }
        #fit2=BLqr(x_input,y_output);
        #plot(density(beta[,1]));
        
        
        
        
      }}}
  return(list(Ajacent_Matrix=AM,Adjusted_Ajacent_Matrix=AM_adjusted,Presence_Probability=S,Standard_Deviations=Std));
}