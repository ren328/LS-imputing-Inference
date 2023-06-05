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
##Use proper function to load elements of \Sigma^*
sd="path"
sdfinal=np.diag(sd)
lambda1=1e-6

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
xxtb=np.matmul(bxt.T,bxt)
b1=np.diag(xxtb)

b2=b1+lambda1
np.fill_diagonal(xxtb,b2)



xtinvb1=np.linalg.inv(xxtb)
xtinvb=np.matmul(xtinvb1,bxt.T)



v1=np.matmul(xtinvb,sdfinal1)

v2=np.matmul(v1,LD)

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

xxtc=np.matmul(cxt.T,cxt)
c1=np.diag(xxtc)
c2=c1+lambda1
np.fill_diagonal(xxtc,c2)


xtinvc1=np.linalg.inv(xxtc)
xtinvc=np.matmul(cxt,xtinvc1)

v=np.matmul(v3,xtinvc)

##Use proper function to save the result
np.save("path", v)
