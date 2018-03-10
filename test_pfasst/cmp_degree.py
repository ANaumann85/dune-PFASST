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
  return [np.linalg.norm(resD),np.amax(np.abs(resD))]

if len(sys.argv) != 2:
  print("usage %s <degree>\n" % argv[0])
  sys.exit(1)
nSteps=[2,4,8,16,32,64, 128, 256]
nIters=[1,2,3,4,8]
nCores=[2,4]

te=20.0
degree=int(sys.argv[1])
resFormat="heat2dMoving_result_pfasst_nCoarse_%(nCoarse)d_nFine_%(nFine)d_nIter_%(nIter)d_nSteps_%(nSteps)d_nCore_%(nCore)d.csv"
#refFile="../../../../movingHeatSDC/movingHeat/movingHeat/heat2d_for_pfasst/degree_%(deg)d/heat_ros2_81920.csv" % {'deg':degree}
refFile="../../../../movingHeatSDC/movingHeat/movingHeat/heat2d_for_pfasst/degree_%(deg)d/heat_ros2_8192.csv" % {'deg':degree}

for nCore in nCores:
  diffData=np.zeros((len(nSteps), len(nIters)+1))
  for nIterP in range(len(nIters)):
    nIter = nIters[nIterP]
    for nStepP in range(len(nSteps)):
      nStep=nSteps[nStepP]
      data={'nIter' : nIter, 'nSteps' : nStep*nCore, 'nCore':nCore, 'nCoarse':1 , 'nFine':1}
      res=resFormat % data
      diff=[0,0]
      if os.path.isfile(res):
        diff=getDiff(refFile, res)
      diffData[nStepP,0]=nStep*nCore
      diffData[nStepP,nIterP+1]=diff[1]
  print np.array2string(diffData,precision=3)
  lines=plt.loglog(te/np.array(nSteps), diffData[:,1:])
  plt.title('nCore=%(nCore)d' % data)
  plt.legend(lines,['nIter=%d' %k for k in nIters])
  plt.show()

