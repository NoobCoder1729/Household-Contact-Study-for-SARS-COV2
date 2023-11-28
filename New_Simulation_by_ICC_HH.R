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


###This function generates truncated Poisson values, truncated at 1.
new.truncated.poisson=function(n,lambda){
  k=1 #initializing
  out=c() #initializing
  while(k<=n){                #loop till we get n samples
    y=rpois(1,lambda)         #generate poisson 
    if(y>1){                  #check if greater than 1
      out=c(out,y)
      k=k+1
    }else{k=k}
  }
  return(out)
}

#new.truncated.poisson.1(10,4,4)




p.val.mse.fn.new.1=function(ind.ct,p,prob.trt,prob.contr,one.percent,rho1){
  
  #registerDoRNG(seed = seed) #this is setting the seed for parallelization, may or may not set up to you
  
  DataStor <- foreach(s = 1:(2*1e3), .combine = "rbind") %dopar% {
    #the above line will create a matrix Datastor in which all output vectors shall be saved
    cluster.ct.1=floor((one.percent*ind.ct)/100) #the number of households with size 1.
    ind.ct.n1=ind.ct-cluster.ct.1 ##number of poeople in households with size more than 1.
    
    
    ###Individual level randomization###
    N=ind.ct ##total number of individuals
    
    ##generating random groups in individuals which forms households
    n.vec=c()#initializaing
    tot.ct=0 #initializing
    l=1 #initializing
    while(tot.ct<ind.ct.n1){ #loop till total number of people less than people in household size more than 1
      n.vec[l]=new.truncated.poisson(1,4)
      tot.ct=sum(n.vec)
      l=l+1
    }
    
    n.vec[(l-1)]=ind.ct.n1-(tot.ct-n.vec[(l-1)])#left over people are assigned to a single household
    
    
    ###Data Matrix
    X=matrix(0,N,p)#initialization
    X[,1]=rep(1,N)
    
    #random treatment assignment
    for(i in 1:N){
      j=sample(1:p,1)
      if(j != 1){
        X[i,j]=1
      }
      
    }
    
    ##Assigning Probability values
    p.vec=rep(0,N)
    p.vec[X[,2]==1]=prob.trt
    p.vec[X[,2]==0]=prob.contr
    
    ##Artificially created pointers for creating responses
    ##the pointers indicate different household groups
    cum.sum=c()#initialization
    cum.sum[1]=n.vec[1] #1st household
    #tagging the other households
    for(i in 2:length(n.vec)){
      cum.sum[i]=cum.sum[i-1]+n.vec[i]
    }
    
    
    ##The response
    #we create the response which are correlated binary values
    y=rep(0,N) #initialization
    #the first line indicates that all responses cannot be 0 i.e. we need at least one infection
    #generating the correlated responses
    while(sum(y)==0){
      y[1:cum.sum[1]]=cor.bin(p.vec[1:cum.sum[1]],rho=rho1,reps=1,k=1)  
      for (i in 1:(length(n.vec)-1)) {
        y[(cum.sum[i]+1):cum.sum[i+1]]=cor.bin(p.vec[(cum.sum[i]+1):cum.sum[i+1]],rho=rho1,reps=1,k=1)  
      }
      
    }
    
    #generating independent responses for single person houses
    for(i in (cum.sum[length(n.vec)]+1):N){
      y[i]=1*(runif(1)<p.vec[i])
    }
    
    ###Subject Index for GEE
    ###Each house is the same subject. 
    ###One can think of the data as Repeated measurements 
    ###from the same Household which is the subject
    
    #households with more than 1 person
    subject.vec=c()
    for(i in 1:length(n.vec)){
      subject.vec=c(subject.vec,rep(i,n.vec[i]))
    }
    
    #households with 1 person
    subject.vec.1=seq(from=(cum.sum[length(n.vec)]+1),to=N,by=1)
    
    subject.vec=c(subject.vec,subject.vec.1)
    
    ##design Matrix
    X.design=cbind(X[,-1],subject.vec)
    data.frame.1=as.data.frame(cbind(y,X.design))
    names(data.frame.1)=c("Resp","Treatment","Subject")
    
    
    ##GEE Model
    g1.mod=geeglm(Resp~Treatment,id=Subject,
                  family = binomial,data=data.frame.1,
                  corstr = "exchangeable",scale.fix = TRUE)
    s.g1=summary(g1.mod)
    
    a.out=anova(g1.mod) #anova
    
    
    check.list=1*(a.out$`P(>|Chi|)`[1]<.05) #checking if p-value more than cutoff
    
    #mse
    #mse.vec.gee=mean((predict(g1.mod,data.frame.1,type = "response")-prob.trt)^2)
    
    ###Logistic Model
    g2.mod=glm(Resp~Treatment,family = binomial,data=data.frame.1)
    
    #p-value more than cutoff or not
    check.vec.logit=1*(coef(summary(g2.mod))[2,4]<.05)
    
    #mse
    #mse.vec.logit=mean((predict(g2.mod,data.frame.1,type = "response")-prob.trt)^2)
    
    #output from loop
    c(check.list, check.vec.logit)
  }
  
  #the averages of the output from loop
  outputMat <- apply(DataStor, 2, mean)
  
  #value returned by function/final output
  return(outputMat)
}


library("xtable")
library(doParallel)
library(doRNG)

#number of cores
registerDoParallel(cores = 4)

#seed = 1729 #use this if you use RNG

#percentage of households with size 1
household.vec=c(10,20,40,70,80)

#p.val.mse.fn.new.1(1450,2,.05,.1,90,.999)

start.time=Sys.time()
rho.vec=c(0.340000,0.430000,0.500000,0.570000,0.690000,0.750000,0.800000,0.870000 )
gamma.vec=c(.15,.20,.25,.30,.40,.45,.50,.60)

df1=data.frame(t(matrix(unlist(lapply(rho.vec,p.val.mse.fn.new.1,ind.ct=1368,p=2,prob.trt=.05,prob.contr=.1,one.percent=30)),nrow=2,ncol=length(rho.vec))))

df2=data.frame(cbind(rho.vec,gamma.vec))

df=cbind(df2,df1)

print(xtable(df,digits=6,type="html"))
end.time=Sys.time()
time=end.time-start.time
time

######
###########QQ Plots##################
#############################################


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

inv.logit=function(x){
  1/(1+exp(-x)) #inverse logit function
}

new.truncated.poisson=function(n,lambda){
  k=1
  out=c()
  while(k<=n){
    y=rpois(1,lambda)
    if(y>1){
      out=c(out,y)
      k=k+1
    }else{k=k}
  }
  return(out)
}


#new.truncated.poisson.1(10,4,4)




p.val.mse.fn.q2=function(ind.ct,p,prob.trt,prob.contr,one.percent,rho1){
  
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
      n.vec[l]=new.truncated.poisson(1,4)
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
    
    ###Artificial pointer for creating responses
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
    
    ###Subject ID
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

df1=p.val.mse.fn.q2(1450,2,prob.trt=.1,prob.contr=.1,one.percent=80,rho1=.999)

par(mfrow=c(1,2))
probs=ppoints(length(df1[,1]))
plot(quantile(runif(length(df1[,1])),probs)
     ,quantile(df1[,1],probs),
     main="QQ GEE",
     xlab = "Uniform Quantiles", ylab="Sample Quantiles")
abline(0,1)

probs=ppoints(length(df1[,2]))
plot(quantile(runif(length(df1[,2])),probs)
     ,quantile(df1[,2],probs),
     main="Q-Q Logistic",
     xlab = "Uniform Quantiles", ylab="Sample Quantiles")
abline(0,1)




















