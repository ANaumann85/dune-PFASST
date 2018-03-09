#!/usr/bin/env python
# vim: ai:ts=2:sw=2:et

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../')

#from diff import getDiff

def getDiff(ref, res):
  refD=np.loadtxt(ref)
  resD=np.loadtxt(res)
  resD -= refD
  return [np.linalg.norm(resD),np.amax(resD)]

nSteps=[2,4,8,16,32,64]
nIters=[1,2,4,8]
nCores=[2,4]

te=20.0
degree=2
resFormat="heat2dMoving_result_pfasst_nCoarse_1_nFine_1_nIter_%(nIter)d_nSteps_%(nSteps)d_nCore_%(nCore)d.csv"
refFile="../../../../movingHeatSDC/movingHeat/movingHeat/heat2d_for_pfasst/degree_%(deg)d/heat_ros2_8192.csv" % {'deg':degree}
for nCore in nCores:
  diffData=np.zeros((len(nSteps), len(nIters)))
  for nIterP in range(len(nIters)):
    nIter = nIters[nIterP]
    for nStepP in range(len(nSteps)):
      nStep=nSteps[nStepP]
      data={'nIter' : nIter, 'nSteps' : nStep*nCore, 'nCore':nCore}
      res=resFormat % data
      diff=[0,0]
      if os.path.isfile(res):
        diff=getDiff(refFile, res)
      diffData[nStepP,nIterP]=diff[1]
  print diffData
  lines=plt.loglog(te/np.array(nSteps), diffData)
  plt.title('nCore=%(nCore)d' % data)
  plt.legend(lines,['nIter=%d' %k for k in nIters])
  plt.show()

