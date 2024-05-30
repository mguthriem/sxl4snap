# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import sys

import SNAPSXLTools as snp
import importlib
importlib.reload(snp)

### experiment data info ###
run = 48028
isLite = False
### -------------------- ###

### sample info ###
a = 18.2928
b = 18.6383
c = 6.5848
alpha = 90
beta = 90
gamma = 90
centering = 'F'
### ----------- ###

#get SNAP config info
instConfig = snp.instConfig("/SNS/SNAP/shared/Calibration_dynamic/SNAPInstPrm.json")
stateID,stateDict,errorState=snp.stateFromRunFunction(run,instConfig)
print(f"stateID: {stateID}")
print(f"stateDict: {stateDict}")

### incident wavelength ###
wavelength_band = [stateDict["wav"]-instConfig.neutronBandwidth/2,
    stateDict["wav"]+instConfig.neutronBandwidth/2]


ipts = GetIPTS(runNumber=run,
    instrument=instConfig.name)
    
if run<=40000:   #TODO: figure out exact run where naming changed
    Fname = f"{ipts}data/SNAP_{run}_event.nxs"
else:
    Fname = f"{ipts}nexus/SNAP_{run}.nxs.h5"
        
if isLite:
    liteName = f"{ipts}shared/lite/SNAP_{run}.lite.nxs.h5"
    snp.makeLite(Fname,liteName)
    LoadEventNexus(Filename=liteName,
    OutputWorkspace=f"TOF_lite_raw_{run}")
    CloneWorkspace(InputWorkspace=f"TOF_lite_raw_{run}",OutputWorkspace="data")
else:
    LoadEventNexus(Filename=Fname,
    OutputWorkspace=f"TOF_raw_{run}")
    CloneWorkspace(InputWorkspace=f"TOF_raw_{run}",OutputWorkspace="data")


SetGoniometer(Workspace='data',
              Axis0='BL3:Mot:omega,0,1,0,1', #because it points up,
              Axis1='BL3:Mot:phi,0.707,0.707,0,1',
              Average=True)

PreprocessDetectorsToMD(InputWorkspace='data',
                        OutputWorkspace='detectors')

#Convert to MD:
two_theta = mtd['detectors'].column('TwoTheta')
Q_max = 4*np.pi/min(wavelength_band)*np.sin(0.5*max(two_theta))
ConvertToMD(InputWorkspace='data',
            QDimensions='Q3D',
            dEAnalysisMode='Elastic',
            Q3DFrames='Q_sample',
            LorentzCorrection=True,
            MinValues=[-Q_max,-Q_max,-Q_max],
            MaxValues=[+Q_max,+Q_max,+Q_max],
            PreprocDetectorsWS='detectors',
            OutputWorkspace='md')