#Use proper function to load SNP matrix
snp="path"

##Use proper function to load trait values
trait="path"

beta<-c()
beta_sd<-c()
n_snp<-c()
p_value<-c()
for(i in 1:ncol(snp)){
  n_snp[i]<-sum(is.na(snp[,i]))
  if(n_snp[i]==0){
    m1<-lm(trait~snp[,i])
    beta[i]<-m1$coef[2]
    beta_sd[i]<-summary(m1)$coef[2,2]
    p_value[i]<-summary(m1)$coef[2,4]
  }else{
    na_id<-which(is.na(snp[,i]==T))
    m1<-lm(trait[-na_id]~snp[-na_id,i])
    beta[i]<-m1$coef[2]
    beta_sd[i]<-summary(m1)$coef[2,2]
    p_value[i]<-summary(m1)$coef[2,4]
  }
}
##Use proper function to save the result
save(n_snp,beta,beta_sd,p_value,file="path")






