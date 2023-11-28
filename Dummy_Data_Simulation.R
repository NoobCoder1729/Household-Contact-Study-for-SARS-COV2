#rm(list=ls())
cor.bin=function(p.vec,rho,reps,k){
  n=length(p.vec)
  Sig.1=matrix(rho,n,n)
  diag(Sig.1)=rep(1,n)#the variance covariance matrix of the Gaussian as compound symmetry
  out=svd(Sig.1) #SVD
  
  z=matrix(rnorm(n*reps),n,reps) #Random std normal generation
  x=out$u%*%diag(out$d^{.5})%*%z #correlated normal generation
  u.vec=pnorm(x) #transformation to uniform
  ber.out=1*(u.vec<p.vec) #the output vector
  ber.out #storage
  return(ber.out[,1:k])
  
}

truncated.poisson.gen=function(n,lambda){
  k=1
  out=c()
  while (k<=n){
    x=rpois(1,lambda)
    if(x==0){
      k=k
    }else{
      out[k]=x
      k=k+1
    }
  }
  return(out)
}

#n=1368;rr.pl=.1;rr.trt=.05;rho=.465

data.gen.trial=function(cor.bin.fn,pois.gen.fn,n,rr.pl,rr.trt,rho.1){
  n.vec=c()
  tot.ct=0 #initializing
  l=1 #initializing
  while(tot.ct<n){ #loop till total number of people less than people in household size more than 1
    n.vec[l]=pois.gen.fn(1,4)
    tot.ct=sum(n.vec)
    l=l+1
  }
  
  n.vec[(l-1)]=n-(tot.ct-n.vec[(l-1)])#left over people are assigned to a single household
  
  trt.assgn=1*(runif(n)<.5)
  
  p.vec=c()
  p.vec[trt.assgn==1]=rr.trt
  p.vec[trt.assgn==0]=rr.pl
  
  ##Artificially created pointers for creating responses
  ##the pointers indicate different household groups
  cum.sum=c()#initialization
  cum.sum[1]=n.vec[1] #1st household
  #tagging the other households
  for(i in 2:n){
    cum.sum[i]=cum.sum[i-1]+n.vec[i]
  }
  
  
  ##The response
  #we create the response which are correlated binary values
  Visit=matrix(0,n,4) #initialization
  #the first line indicates that all responses cannot be 0 i.e. we need at least one infection
  #generating the correlated responses
  for(j in 1:4){
    while(sum(Visit[,j])==0){
      Visit[1:cum.sum[1],j]=cor.bin.fn(p.vec[1:cum.sum[1]],rho=rho.1,reps=1,k=1)  
      for (i in 1:(length(n.vec)-1)) {
        Visit[(cum.sum[i]+1):cum.sum[i+1],j]=cor.bin.fn(p.vec[(cum.sum[i]+1):cum.sum[i+1]],rho=rho.1,reps=1,k=1)  
      }
      
    }
  }
  
  subject.vec=c()
  for(i in 1:length(n.vec)){
    subject.vec=c(subject.vec,rep(i,n.vec[i]))
  }
  
  age=1*(runif(n)<.66)
  
  out.data=data.frame(cbind(trt.assgn,age,subject.vec,Visit))
  names(out.data)=c("Treatment_Assigned","Age","Household_Number","Visit_1","Visit_2","Visit_3","Visit_4")
  
  return(out.data)
  
}

df1=data.gen.trial(cor.bin.fn = cor.bin,pois.gen.fn = truncated.poisson.gen,n=1368,rr.pl=.1,rr.trt=.05,rho.1 =.465)

#install.packages("writexl")
library("writexl")
write_xlsx(df1,"C:\\Users\\riddh.bhattacharya\\Documents\\Simulation_Docs\\Study_Data_Frame.xlsx")
