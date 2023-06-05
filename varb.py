import math
import scipy
import scipy.linalg 
import numpy as np
import os
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

##batch number
i=
##Use proper function to load SNP matrix
snp="path"
snp=np.array(snp)
##Use proper function to load elements of \Sigma^*
sd="path"
sdfinal=np.diag(sd)
lam=1e-6

p=np.shape(snp)[1]
res=[]
res2=[]
for j in range(p):
  c=snp[:,j]
  cp=np.matmul(c.T,c)
  cp1=1/cp
  cp2=math.sqrt(cp1)
  res.append(cp1)
  res2.append(cp2)

res3=np.diag(res2)
res1=np.diag(res)
cxt=np.matmul(res1,snp.T)
cxt3=np.matmul(res3,snp.T)

##approximate LD matrix
LD=np.matmul(cxt3,cxt3.T)

#xtinv=np.linalg.pinv(cxt,rcond=math.sqrt(np.finfo(float).eps))

xxt=np.matmul(cxt.T,cxt)
a1=np.diag(xxt)
a2=a1+lam
np.fill_diagonal(xxt,a2)
xxtinv=np.linalg.inv(xxt)
xtinv=np.matmul(xxtinv,cxt.T)



v1=np.matmul(xtinv,sdfinal)
v2=np.matmul(v1,LD)
v3=np.matmul(v2,sdfinal)
v=np.matmul(v3,xtinv.T)

##Use proper function to save the result
np.save("path".format(i), v)
