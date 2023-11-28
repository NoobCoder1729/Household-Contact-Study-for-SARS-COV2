rm(list=ls())
#install.packages("geepack")
library("geepack")

###function to generate correlated bernoulli
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


###Truncated Poisson
#### where there is an upper bound for the realizations
new.truncated.poisson.1=function(n,lambda,u.b){
  k=1
  out=c()
  while(k<=n){
    y=rpois(1,lambda)
    if((y>1)&&(y<=u.b)){
      out=c(out,y)
      k=k+1
    }else{k=k}
  }
  return(out)
}

#new.truncated.poisson.1(10,4,4)




p.val.mse.fn.b1=function(ind.ct,p,prob.trt,prob.contr,one.percent,rho1,u.b1){
  
  #registerDoRNG(seed = seed)
  
  DataStor <- foreach(s = 1:(2*1e3), .combine = "rbind") %dopar% {
    #k=1
    cluster.ct.1=floor((one.percent*ind.ct)/100)
    ind.ct.n1=ind.ct-cluster.ct.1
    #k=1
    
    ###Individual level randomization###
    N=ind.ct
    
    n.vec=c()
    tot.ct=0
    l=1
    while(tot.ct<ind.ct.n1){
      n.vec[l]=new.truncated.poisson.1(1,4,u.b1)
      tot.ct=sum(n.vec)
      l=l+1
    }
    
    n.vec[(l-1)]=ind.ct.n1-(tot.ct-n.vec[(l-1)])
    
    
    X=matrix(0,N,p)
    X[,1]=rep(1,N)
    
    for(i in 1:N){
      j=sample(1:p,1)
      if(j != 1){
        X[i,j]=1
      }
      
    }
    
    ##Assigning Probabilities
    p.vec=rep(0,N)
    p.vec[X[,2]==1]=prob.trt
    p.vec[X[,2]==0]=prob.contr
    
    ##Artificial Pointer for Responses
    cum.sum=c()
    cum.sum[1]=n.vec[1]
    for(i in 2:length(n.vec)){
      cum.sum[i]=cum.sum[i-1]+n.vec[i]
    }
    
    ##The Response
    y=rep(0,N)
    while(sum(y)==0){
      y[1:cum.sum[1]]=cor.bin(p.vec[1:cum.sum[1]],rho=rho1,reps=1,k=1)  
      for (i in 1:(length(n.vec)-1)) {
        y[(cum.sum[i]+1):cum.sum[i+1]]=cor.bin(p.vec[(cum.sum[i]+1):cum.sum[i+1]],rho=rho1,reps=1,k=1)  
      }
      
    }
    
    for(i in (cum.sum[length(n.vec)]+1):N){
      y[i]=1*(runif(1)<p.vec[i])
    }
    
    ##subject ID
    subject.vec=c()
    for(i in 1:length(n.vec)){
      subject.vec=c(subject.vec,rep(i,n.vec[i]))
    }
    
    subject.vec.1=seq(from=(cum.sum[length(n.vec)]+1),to=N,by=1)
    
    subject.vec=c(subject.vec,subject.vec.1)
    
    ##Design Matrix
    X.design=cbind(X[,-1],subject.vec)
    data.frame.1=as.data.frame(cbind(y,X.design))
    #data.frame.1=data.frame.2[sample(nrow(data.frame.2)),]
    names(data.frame.1)=c("Resp","Treatment","Subject")
    
    ##GEE Model
    g1.mod=geeglm(Resp~Treatment,id=Subject,
                  family = binomial,data=data.frame.1,
                  corstr = "exchangeable",scale.fix = TRUE)
    s.g1=summary(g1.mod)
    
    a.out=anova(g1.mod)
    
    
    check.list=1*(a.out$`P(>|Chi|)`[1]<.05)
    
    mse.vec.gee=mean((predict(g1.mod,data.frame.1,type = "response")-prob.trt)^2)
    
    ##Logistic Model
    g2.mod=glm(Resp~Treatment,family = binomial,data=data.frame.1)
    
    check.vec.logit=1*(coef(summary(g2.mod))[2,4]<.05)
    
    mse.vec.logit=mean((predict(g2.mod,data.frame.1,type = "response")-prob.trt)^2)
    
    c(check.list, check.vec.logit, mse.vec.gee, mse.vec.logit)
  }
  
  outputMat <- apply(DataStor, 2, mean)
  
  return(outputMat)
}


library("xtable")
library(doParallel)
library(doRNG)

registerDoParallel(cores = 3)

seed = 1729

household.vec=c(5,10,20,30,40)
u.b1=4

start.time=Sys.time()

df1=data.frame(t(matrix(unlist(lapply(household.vec, p.val.mse.fn.b1,ind.ct=1368,p=2,prob.trt=.1,prob.contr=.1,rho1=.999,u.b1)),nrow=4,ncol=length(household.vec))))

print(xtable(df1,digits=6,type="html"))
end.time=Sys.time()
time=end.time-start.time
time


#############################################
#########Q-Q plots###########################

rm(list=ls())
#install.packages("geepack")
library("geepack")

###function to generate correlated bernoulli
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

###Truncated Poisson
#### where there is an upper bound for the realizations
new.truncated.poisson.1=function(n,lambda,u.b){
  k=1
  out=c()
  while(k<=n){
    y=rpois(1,lambda)
    if((y>1)&&(y<=u.b)){
      out=c(out,y)
      k=k+1
    }else{k=k}
  }
  return(out)
}



p.val.mse.fn.q1=function(ind.ct,p,prob.trt,prob.contr,one.percent,rho1,u.b1){
  
  #registerDoRNG(seed = seed)
  
  DataStor <- foreach(s = 1:(2*1e3), .combine = "rbind") %dopar% {
    #k=1
    cluster.ct.1=floor((one.percent*ind.ct)/100)
    ind.ct.n1=ind.ct-cluster.ct.1
    #k=1
    
    ###Individual level randomization###
    N=ind.ct
    
    n.vec=c()
    tot.ct=0
    l=1
    while(tot.ct<ind.ct.n1){
      n.vec[l]=new.truncated.poisson.1(1,4,u.b1)
      tot.ct=sum(n.vec)
      l=l+1
    }
    
    n.vec[(l-1)]=ind.ct.n1-(tot.ct-n.vec[(l-1)])
    
    
    X=matrix(0,N,p)
    X[,1]=rep(1,N)
    
    for(i in 1:N){
      j=sample(1:p,1)
      if(j != 1){
        X[i,j]=1
      }
      
    }
    
    ##Assigning probabilities
    p.vec=rep(0,N)
    p.vec[X[,2]==1]=prob.trt
    p.vec[X[,2]==0]=prob.contr
    
    ##Artificial Pointer for creating Responses
    cum.sum=c()
    cum.sum[1]=n.vec[1]
    for(i in 2:length(n.vec)){
      cum.sum[i]=cum.sum[i-1]+n.vec[i]
    }
    
    ##The Responses
    y=rep(0,N)
    while(sum(y)==0){
      y[1:cum.sum[1]]=cor.bin(p.vec[1:cum.sum[1]],rho=rho1,reps=1,k=1)  
      for (i in 1:(length(n.vec)-1)) {
        y[(cum.sum[i]+1):cum.sum[i+1]]=cor.bin(p.vec[(cum.sum[i]+1):cum.sum[i+1]],rho=rho1,reps=1,k=1)  
      }
      
    }
    
    for(i in (cum.sum[length(n.vec)]+1):N){
      y[i]=1*(runif(1)<p.vec[i])
    }
    
    ##Subject ID
    subject.vec=c()
    for(i in 1:length(n.vec)){
      subject.vec=c(subject.vec,rep(i,n.vec[i]))
    }
    
    subject.vec.1=seq(from=(cum.sum[length(n.vec)]+1),to=N,by=1)
    
    subject.vec=c(subject.vec,subject.vec.1)
    
    ##Design Matrix
    X.design=cbind(X[,-1],subject.vec)
    data.frame.1=as.data.frame(cbind(y,X.design))
    #data.frame.1=data.frame.2[sample(nrow(data.frame.2)),]
    names(data.frame.1)=c("Resp","Treatment","Subject")
    
    
    ####GEE Model
    g1.mod=geeglm(Resp~Treatment,id=Subject,
                  family = binomial,data=data.frame.1,
                  corstr = "exchangeable",scale.fix = TRUE)
    s.g1=summary(g1.mod)
    
    a.out=anova(g1.mod)
    
    
    check.list=a.out$`P(>|Chi|)`[1]
    
    
    ###Logistic Model
    g2.mod=glm(Resp~Treatment,family = binomial,data=data.frame.1)
    
    check.vec.logit=coef(summary(g2.mod))[2,4]
    
    
    
    c(check.list, check.vec.logit)
  }
  
  
  return(DataStor)
}

library(doParallel)
library(doRNG)

seed=1729

registerDoParallel(cores = 3)

df1=p.val.mse.fn.q1(1368,2,prob.trt=.1,prob.contr=.1,one.percent=40,rho1=.99,u.b1=2)

par(mar=c(5,5,2,2))
probs=ppoints(length(df1[,1]))
plot(quantile(runif(length(df1[,1])),probs)
     ,quantile(df1[,1],probs),
     main="QQ GEE",
     xlab = "Uniform Quantiles", ylab="Sample Quantiles",
     type="l",col="blue")
lines(quantile(runif(length(df1[,2])),probs)
      ,quantile(df1[,2],probs),type="l",col="red")
abline(0,1)

legend(0, 1, legend=c(" GEE", "Logistic"),
       col=c("blue","red"), lty=1:1, cex=0.6)












