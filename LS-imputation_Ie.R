##Use proper function to load var and beta
##beta1, var1 is a p*m matrix, m is the number of batches 
beta1="path"
var1="path"
m=ncol(beta1)

beta<-apply(beta1,1,mean)
sd<-sqrt(apply(var1,1,sum)/m^2)
z_stat=beta/sd
pie<-2*pnorm(-abs(z_stat))


##Use proper function to save the result
save(pie,"path")
