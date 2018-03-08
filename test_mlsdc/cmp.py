#!/usr/bin/env python
# vim: ai:ts=2:sw=2:et

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../')

from diff import getDiff

nSteps=[2,4,8,16,32,64]
nIters=[1,2,4,8]
nFines=[1,2,4]
nCoarses=[1,2,4]

resFormat="heat2dMoving_result_mlsdc_nCoarse_%(nCoarse)d_nFine_%(nFine)d_nIter_%(nIter)d_nSteps_%(nSteps)d.vtu"
refFile="../../../movingHeat_from_taurus/movingHeat_pfasst_ros2/movingHeat/movingHeat/heat2d_for_pfasst/heat_ros2_8192-00001.vtu"
for nFine in nFines:
  for nCoarse in nCoarses:
    diffData=np.zeros((len(nSteps), len(nIters)))
    for nIterP in range(len(nIters)):
      nIter = nIters[nIterP]
      for nStepP in range(len(nSteps)):
        nStep=nSteps[nStepP]
        data={'nFine':nFine, 'nCoarse' : nCoarse, 'nIter' : nIter, 'nSteps' : nStep}
        res=resFormat % data
        diff=getDiff(refFile, res)
        diffData[nStepP,nIterP]=diff[1]
    #print diffData
    lines=plt.loglog(nSteps, diffData)
    plt.title('nCoarse=%(nCoarse)d, nFine=%(nFine)d' % data)
    plt.legend(lines,['nIter=%d' %k for k in nIters])
    plt.show()

