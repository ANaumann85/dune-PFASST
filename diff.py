#!/usr/bin/env python
from paraview.simple import *
import math, sys

def getDiff(coarseName, fineName, comp=0):
    #print "current fine:", fineName
    #print "current coarse:", coarseName
    #read fine
    finePVD=XMLUnstructuredGridReader(FileName=fineName) #PVDReader(FileName=fineName)
    #read coarse
    coarsePVD=XMLUnstructuredGridReader(FileName=coarseName) #PVDReader(FileName=coarseName)
    #request number of timesteps
    #rename values
    #request original value names
    fineArName=finePVD.GetPointDataInformation().GetArray(comp).GetName()
    coarseArName=coarsePVD.GetPointDataInformation().GetArray(comp).GetName()
    fineCalcName="fine_"+fineArName
    coarseCalcName="coarse_"+coarseArName

    #create calculators
    fineCalc=Calculator()
    fineCalc.Input=finePVD
    fineCalc.Function=fineArName
    fineCalc.ResultArrayName=fineCalcName

    coarseCalc=Calculator()
    coarseCalc.Input=coarsePVD
    coarseCalc.Function=coarseArName
    coarseCalc.ResultArrayName=coarseCalcName

    #resample contains coarse data only! append them
    appAtt=AppendAttributes()
    appAtt.Input=[coarseCalc, fineCalc]

    #compute difference
    diffCalc=Calculator()
    diffCalc.Input=appAtt
    diffCalc.ResultArrayName="diff"

    diffCalc.Function="("+fineCalcName+"-"+coarseCalcName+")"
    diffCalc.Function=diffCalc.Function+"*"+diffCalc.Function

    #integrate the differences
    integ=IntegrateVariables()
    #request difference in last timestep

    resAr=integ.GetPointDataInformation()
    integral=resAr.GetArray("diff").GetRange()[0]
    l2Norm=0
    if integral > 0:
            l2Norm =math.sqrt(integral)

    maxCalc=Calculator()
    maxCalc.Input=appAtt;
    maxCalc.ResultArrayName="diff"
    maxCalc.Function="abs("+fineCalcName+"-"+coarseCalcName+")"
    resAr=maxCalc.GetPointDataInformation()
    maxNorm=resAr.GetArray("diff").GetRange()[1]
    #writer=XMLUnstructuredGridWriter(FileName="diff.vtu")
    #writer.UpdatePipeline()
    return [l2Norm, maxNorm]

if __name__ == "__main__":
    if len(sys.argv) < 3:
            print "usage: ", sys.argv[0], "<vtu1> <vtu2> [comp]"
            exit(1)

    coarseName=sys.argv[1]
    fineName=sys.argv[2]

    comp=0
    if len(sys.argv) >= 4:
        comp = int(sys.argv[3])
    res=getDiff(coarseName, fineName, comp)
    print "l2 norm:", res[0] #l2Norm
    print "maxNorm:", res[1] #maxNorm
