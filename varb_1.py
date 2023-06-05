import math
import scipy
import scipy.linalg 
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

##batch number
i=
##Use proper function to load $Var(Y_(b)|X,X^*)$
varb=np.load("path")
##Use proper function to load $\hat_sigma(b)^2$
sigmab="path"
sigmab=np.array(sigmab)
##Use proper functin to load SNP matrix
snp="path"
snp=np.array(snp)

vary=np.diag(varb)
vary1=np.diag(1/np.sqrt(vary))

bp=np.matmul(vary1,varb)
rb=np.matmul(bp,vary1)


p=snp.shape[1]
##D,p*p
res=[]
for j in range(p):
  print(j)
  c=snp[:,j]
  cp=np.matmul(c.T,c)
  cp1=1/cp
  cp2=cp1*c.T*sigmab[j]
  cp3=np.matmul(cp2,rb)
  cp4=np.matmul(cp3,c)
  cp5=cp4*cp1
  res.append(cp5)

##Use proper function to save the result
np.savetxt("path",res)

