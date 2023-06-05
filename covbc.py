import math
import scipy
import scipy.linalg 
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


## Use proper function to load SNP matrix for batch b and batch c
snpb="path"
snpc="path"
##Use proper function to load elements of \sigma^*
sd="path"
sdfinal=np.diag(sd)

pb=np.shape(snpb)[1]
resb=[]
resb2=[]
for j in range(pb):
  b=snpb[:,j]
  bp=np.matmul(b.T,b)
  bp1=1/bp
  bp2=math.sqrt(bp1)
  resb.append(bp1)
  resb2.append(bp2)

resb3=np.diag(resb2)
resb1=np.diag(resb)

bxt=np.matmul(resb1,snpb.T)
bxt3=np.matmul(resb3,snpb.T)


LD=np.matmul(bxt3,bxt3.T)

xxt=np.matmul(bxt.T,bxt)
a1=np.diag(xxt)
lambda1=1e-6
a2=a1+lambda1
np.fill_diagonal(xxt,a2)


##n*n, first term
xtinv1=np.linalg.inv(xxt)
xtinv=np.matmul(xtinv1,bxt.T)
print(xtinv)

##first term times sigma
v1=np.matmul(xtinv,sdfinal1)
##first term times sigma * LD
v2=np.matmul(v1,LD)
## above*sigma
v3=np.matmul(v2,sdfinal1)



pc=np.shape(snpc)[1]
resc=[]
for k in range(p):
  c=snpc[:,k]
  cp=np.matmul(c.T,c)
  cp1=1/cp
  resc.append(cp1)
resc1=np.diag(resc)

cxt=np.matmul(resc1,snpc.T)

xxt=np.matmul(cxt.T,cxt)
a1=np.diag(xxt)
lambda1=1e-6
a2=a1+lambda1
np.fill_diagonal(xxt,a2)


xtinv1=np.linalg.inv(xxt)
xtinv=np.matmul(cxt,xtinv1)

v=np.matmul(v3,xtinv)

##Use proper function to save the result
np.save("path", v)
