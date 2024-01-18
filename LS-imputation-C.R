##total number of batches
m=
##specify the df of each batch
df=

var<-list()
beta<-list()
t_stat<-list()
p<-list()

## load var, beta of all the batches
for(j in 1:m){
  var[[j]]<-as.numeric(unlist(read.table(file=paste0("path/var",j,".txt"))))
  beta[[j]]<-as.numeric(unlist(read.table(file=paste0("path/beta",j,".txt"))))
  t_stat[[j]]<-beta[[j]]/sqrt(var[[j]])
  df_final=df[j]
  p[[j]]<-2*pt(-abs(t_stat[[j]]), df=df_final)
}




##IVW of beta
beta1<-do.call(cbind,beta)
var1<-do.call(cbind,var)
betafinal<-c()
for(i in 1:nrow(beta1)){
  betafinal[i]<-sum(beta1[i,]/var1[i,])/sum(1/var1[i,])
}
betaivw<-betafinal


##combine all the p-value
p1<-do.call(cbind,p)

##cauchy method
tc=c()
for(i in 1:nrow(p1)){
  tc[i]=sum(1/ncol(p1)*tan((0.5-p1[i,])*pi))
}

pcauchy=0.5-atan(tc)/pi


 ##averaging method
##r=B-1
r=(ncol(p1)-1)
prb_1=c()
for(i in 1:nrow(p1)){
  prb_1[i]=ncol(p1)^{1/r}*{sum(p1[i,]^r)/ncol(p1)}^{1/r}
}

##r=1
pr1=2*apply(p1,1,mean)

##r=inf
prinf=apply(p1,1,max)



##r=-inf
pr_inf=ncol(p1)*apply(p1,1,min)





##r=1/B-1
r=1/(ncol(p1)-1)
pr1b_1=c()
for(i in 1:nrow(p1)){
  pr1b_1[i]=(r+1)^{1/r}*{sum(p1[i,]^r)/ncol(p1)}^{1/r}
}




##r=B+1
r=(ncol(p1)+1)
prbadd1=c()
for(i in 1:nrow(p1)){
  prbadd1[i]=ncol(p1)^{1/r}*{sum(p1[i,]^r)/ncol(p1)}^{1/r}
}


##r=B+5
r=(ncol(p1)+5)
prbadd5=c()
for(i in 1:nrow(p1)){
  prbadd5[i]=ncol(p1)^{1/r}*{sum(p1[i,]^r)/ncol(p1)}^{1/r}
}



##r=6/B-1
r=6/(ncol(p1)-1)
pr6b_1=c()
for(i in 1:nrow(p1)){
  pr6b_1[i]=(r+1)^{1/r}*{sum(p1[i,]^r)/ncol(p1)}^{1/r}
}




##r=B-2
r=(ncol(p1)-2)
prb_2=c()
for(i in 1:nrow(p1)){
  prb_2[i]=(r+1)^{1/r}*{sum(p1[i,]^r)/ncol(p1)}^{1/r}
}


##Use proper function to save the result
save(betaivw,pcauchy,prb_1,pr1,prinf,pr_inf,pr1b_1,prbadd1,prbadd5,pr6b_1,prb_2,"path")



