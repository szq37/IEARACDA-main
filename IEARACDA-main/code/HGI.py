import copy
import numpy as np
def HG_Inference(Y,circRNA_si,disease_si,alpha):
# the normalization of circRNA
  circRNA_si_1=copy.deepcopy(circRNA_si)
  for mm1 in range(m):
      for mm2 in range(m):
          circRNA_si[mm1,mm2]=circRNA_si[mm1,mm2]/(np.sqrt(np.sum(circRNA_si_1[mm1,:]))*np.sqrt(np.sum(circRNA_si_1[mm2,:])))
#the normalization of disease
  disease_si_1=copy.deepcopy(disease_si)
  for nn1 in range(n):
      for nn2 in range(n):
          disease_si[nn1,nn2]=disease_si[nn1,nn2]/(np.sqrt(np.sum(disease_si_1[nn1,:]))*np.sqrt(np.sum(disease_si_1[nn2,:])))
#calculate the score matrix Y
  S=Y
  Si=alpha*np.dot(np.dot(circRNA_si, S), disease_si)+(1-alpha)*Y
  while np.linalg.norm(Si-S,1)>10**-6:
      S=Si
      Si=alpha*np.dot(np.dot(circRNA_si, S), disease_si)+(1-alpha)*Y
  Y_recovery = Si
  return Y_recovery
