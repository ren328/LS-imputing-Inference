import math
import scipy
import scipy.linalg 
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

##batch number
batch=1 
##Use proper function to load $Var(Y_(b)|X,X^*)$
varb=np.load("path")
##Use proper function to load $\hat_sigma(b)^2$
sigmab="path"

vary=np.diag(varb)
vary1=np.diag(1/np.sqrt(vary))

cp=np.matmul(vary1,varb)
R=np.matmul(cp,vary1)






##centered SNP matrix 
r = robjects.r
r['source']('/panfs/jay/groups/20/panwei/ren00075/0515newvar/20000py/batch1/f1.R')
rfunction = robjects.globalenv['f1']
with localconverter(robjects.default_converter + pandas2ri.converter):
     snp= rfunction(chr)

r = robjects.r
r['source']('/panfs/jay/groups/20/panwei/ren00075/0515newvar/20000py/batch1/f1.R')
rfunction = robjects.globalenv['f2']
with localconverter(robjects.default_converter + pandas2ri.converter):
     sigma= rfunction(chr)

snp=np.array(snp)
sigma=np.array(sigma)
print(sigma.shape)
print(snp.shape)
print(type(snp))

p=snp.shape[1]
##D,p*p
res=[]
for j in range(p):
  print(j)
  c=snp[:,j]
  cp=np.matmul(c.T,c)
  cp1=1/cp
  cp2=cp1*c.T*sigma[j]
  cp3=np.matmul(cp2,rb)
  cp4=np.matmul(cp3,c)
  cp5=cp4*cp1
  res.append(cp5)

np.savetxt("/panfs/jay/groups/20/panwei/ren00075/0515newvar/20000py/batch1/var{}.txt".format(chr),res)

