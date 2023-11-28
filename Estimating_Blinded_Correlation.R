##Function to estimate Correlation
##Pearson Method
pearson.cor.fn=function(y,hh.tag){
  if(length(y)==length(hh.tag)){
    Z=aggregate(y,by=list(Category=hh.tag),FUN=sum)[,2]
    one=rep(1,length(y))
    n=aggregate(one,by=list(Category=hh.tag),FUN=sum)[,2]
    mu=sum(Z*(n-1))/sum(n*(n-1))
    rho=(1/(mu*(1-mu)))*((sum(Z*(Z-1))/sum(n*(n-1)))-mu^2)
    return(rho)
  }else{
    return("Error")
  }
}

y.1=1*(runif(30)<.5)
hh.tag.1=c()
for(i in 1:10){
  hh.tag.1=c(rep(i,3),hh.tag.1)
}

pearson.cor.fn(y.1,hh.tag.1)
#length(y)==length(hh.tag)














