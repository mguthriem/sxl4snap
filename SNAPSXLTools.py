#This is an assembly of tools to manage SNAP's reconfigurability, modelled on 
#SNAPRed work.

import json
from mantid.simpleapi import *


class instConfig:
  
  def __init__(self,instJsonPath):
    
    try:
        with open(instJsonPath, "r") as json_file:
            inst = json.load(json_file)
    except:
        print(f"error opening inst config file: {instJsonPath}")
        return
    
    #copy a subset of parameters from input json file to matching class attributes
    self.version = inst["version"]
    self.facility = inst["facility"]
    self.name = inst["name"]
    self.nexusFileExtension = inst["nexusFileExtension"]
    self.calibrationDirectory = inst["calibrationDirectory"]
    self.calibrationFilePrefix = inst["calibrationFilePrefix"]
    self.pixelGroupingDirectory = inst["pixelGroupingDirectory"]
    self.simpleContainerDirectory = inst["simpleContainerDirectory"]
    self.sharedDirectory = inst["sharedDirectory"]
    self.nexusDirectory = inst["nexusDirectory"]
    self.neutronBandwidth = inst["neutronBandwidth"]
    self.extendedNeutronBandwidth = inst["extendedNeutronBandwidth"]
    self.L1 = inst["L1"]
    self.L2 = inst["L2"]
       

def stateFromRunFunction(runNum,inst):
  #returns stateID and info given a run numner
  #
  #20220809 - modified error handling. Instead of requesting input from within
  #function, added an additional variable errorState, with the following meaning
  #
  #20230322 - modified to also generate SNAPRed style state ID using hash
  #
  #errorState = 0 : all is well in the World
  #errorState = 1 : GetIPTS failed. IPTS doesn't exist
  #errorState = 2 : Insufficient log information in nexus file to determine instrument state 

  import h5py
  from os.path import exists
  import sys


  #print(iPrm.inst)
  #initialise errorState dictionary
  errorState=dict()
  errorState['value']=0
  errorState['function']='StateFromRunFunction'
  errorState['message']=''
  errorState['parameters']=[]

#   name= inst.name
  nexusLoc = inst.nexusDirectory
  nexusExt = inst.nexusFileExtension

#   print(f"calling GetIPTS with ")
  try:
    IPTSLoc = GetIPTS(RunNumber=runNum,Instrument=inst.name)
  except:
    errorState['value']=1
    errorState['message']='mantid GetIPTS algorithm failed'
    errorState['parameters']=[runNum]
    return '00000-00',dict(),errorState

  fName = IPTSLoc + nexusLoc + '/SNAP_' + str(runNum) + nexusExt
  print(fName)
  if exists(fName):
    f = h5py.File(fName, 'r')
  else:
    errorState['value']=2
    errorState['message']='error opening run nexus file'
    errorState['parameters']=[fName]
    return '00000-00',dict(),errorState

  fail = False

  stateDict = dict() #dictionary to store state variable values

  missingLogVal = []

  try:
    det_arc1 = f.get('entry/DASlogs/det_arc1/value')[0]
    stateDict['det_arc1']=det_arc1
  except:
    missingLogVal.append('det_arc1')
    fail = True
  try:    
    det_arc2 = f.get('entry/DASlogs/det_arc2/value')[0]
    stateDict['det_arc2']=det_arc2
  except:
    missingLogVal.append('det_arc2')
    fail = True
  try:
    wav = f.get('entry/DASlogs/BL3:Chop:Skf1:WavelengthUserReq/value')[0]
    stateDict['wav']=wav
  except:
    missingLogVal.append('wav')
    fail = True
  try:
    freq = f.get('entry/DASlogs/BL3:Det:TH:BL:Frequency/value')[0]
    stateDict['freq']=freq
  except:
    missingLogVal.append('freq')
    fail = True
  try:
    GuideIn = f.get('entry/DASlogs/BL3:Mot:OpticsPos:Pos/value')[0]
    stateDict['GuideStat']=GuideIn
  except:
    missingLogVal.append('GuideStat')
    fail = True

  if not fail:
    # stateID,errorState = checkSNAPState([det_arc1,det_arc2,wav,freq,0.0],[GuideIn,0])
    stateID = genSNAPState([det_arc1,det_arc2,wav,freq,0.0],[GuideIn,0])
  else:
    errorState['value']=3
    errorState['message']='Insufficient log data, can\'t determine state'
    errorState['parameters']=missingLogVal
    f.close()
    return '00000-00',dict(),errorState
  f.close()
  return stateID,stateDict,errorState

def genSNAPState(floatVars,IntVars):

  import hashlib
  import json

  det_arc1 = floatVars[0]
  det_arc2 = floatVars[1]
  wav = floatVars[2]
  freq = floatVars[3]
  GuideIn = IntVars[0]

  from dataclasses import dataclass

  # https://docs.python.org/3/library/dataclasses.html
  @dataclass
  class StateId:
      vdet_arc1: float
      vdet_arc2: float
      WavelengthUserReq: float
      Frequency: int
      Pos: int

      # Round inputs to reduced number of possible states
      def __init__(self, vdet_arc1: float, vdet_arc2: float, WavelengthUserReq: float, Frequency: float, Pos: int):
          self.vdet_arc1 = float(round(vdet_arc1 * 2) / 2)
          self.vdet_arc2 = float(round(vdet_arc2 * 2) / 2)
          self.WavelengthUserReq = float(round(WavelengthUserReq, 1))
          self.Frequency = int(round(Frequency))
          self.Pos = int(Pos)


  stateID = StateId(vdet_arc1=det_arc1, vdet_arc2=det_arc2, WavelengthUserReq=wav, Frequency =freq, Pos=GuideIn)
  hasher = hashlib.shake_256()

  decodedKey = json.dumps(stateID.__dict__).encode('utf-8')

  hasher.update(decodedKey)

  hashedKey = hasher.digest(8).hex()
  # print(hashedKey)
  return hashedKey

def makeLite(inFileName,outFileName):    

    #makeLite accepts the full path of an input nxs.h5 file and outputs a lite version in the 
    #shared directory of the original files IPTS

    import numpy as np
    import h5py
    import shutil
    import time
    import os
    print('h5py version:',h5py.__version__)

    time_0 = time.time()

    sumNeigh = [8,8]

    #step 0) check if lite file already exists
    if os.path.exists(outFileName):
      print('requested output file already exists!')
      return 

    #
    # Step 1A) create output directory if it doesn't exist
    #
    if not os.path.exists(os.path.dirname(outFileName)):
      os.mkdir(os.path.dirname(outFileName))
    #
    # Step 1B) make a copy of the original file
    #
    print('Making Lite version of file:')
    fSize = os.path.getsize(inFileName)
    print(f'  Copying original .nxs.h5 file: {inFileName} size: {fSize/1e9:.4f} Gb')
    shutil.copyfile(inFileName,outFileName)

    stepTime = time.time() - time_0
    totTime = stepTime
    print(f'    Time to complete step: {stepTime:3.1f} sec. Total time to execute: {totTime:3.1f}')

    print(f'    Relabelling pixel IDs')
    h5obj = h5py.File(outFileName,'r+')

    #
    # Step 2) Relabel pixel IDs
    #
    time_0 = time.time()

    #Create list of SNAP panel IDs
    detpanel = []
    for i in range(1,7):
        for j in range(1,4):
            detpanel.append(str(i)+str(j))
            
    #Relabel pixel IDs and write to copied file
    for panel in detpanel:
        h5eventIDs = h5obj[f'entry/bank{panel}_events/event_id']
        eventIDs = np.array(h5eventIDs[:])
        superEventIDs = superID(eventIDs,sumNeigh[0],sumNeigh[0])
        h5eventIDs[...]=superEventIDs

    stepTime = time.time() - time_0
    totTime += stepTime
    print(f'    Time to complete step: {stepTime:3.1f}s')
    print(f'    Total time to generate Lite file: {totTime:3.1f}')

    #
    # Step 3: Update instrument definition in nxs file
    #
    time_0 = time.time()

    print('Updating instrument definition')
    h5IDF = h5obj['entry/instrument/instrument_xml/data']
    stringIDF = str(h5IDF[0],encoding='ascii')#[:] #a string containing the IDF
    lines = stringIDF.split('\n')
    newLines = []
    for line in lines:
        if '<component type="panel"' in line:
            splitLine = line.split("\"")
    #        print('Original pixel numbering:',int(splitLine[3]),'step by row:',int(splitLine[7]))
            idstart = int(splitLine[3])
            idstepbyrow = int(splitLine[7])
            newidstart=str(int((idstart/65536)*32**2))
            newidstepbyrow = str(int(idstepbyrow/sumNeigh[0]))
            splitLine[3]=newidstart
            splitLine[7]=newidstepbyrow
            newLines.append("\"".join(splitLine))
    #        print('New line:\n',"\"".join(splitLine))
        elif 'xpixels=' in line:
            splitLine = line.split("\"")
            splitLine[1]="32"
            splitLine[3]="-0.076632"
            splitLine[5]="+0.004944"
            newLines.append("\"".join(splitLine))
        elif 'ypixels=' in line:
            splitLine = line.split("\"")
            splitLine[1]="32"
            splitLine[3]="-0.076632"
            splitLine[5]="+0.004944"
            newLines.append("\"".join(splitLine))
        elif 'left-front-bottom-point' in line:
            splitLine = line.split("\"")
            splitLine[1]='-0.002472'
            splitLine[3]='-0.002472'
            newLines.append("\"".join(splitLine))
        elif 'left-front-top-point' in line:
            splitLine = line.split("\"")
            splitLine[1]='0.002472'
            splitLine[3]='-0.002472'
            newLines.append("\"".join(splitLine))
        elif 'left-back-bottom-point' in line:
            splitLine = line.split("\"")
            splitLine[1]='-0.002472'
            splitLine[3]='-0.002472'
            newLines.append("\"".join(splitLine))
        elif 'right-front-bottom-point' in line:
            splitLine = line.split("\"")
            splitLine[1]='0.002472'
            splitLine[3]='-0.002472'
            newLines.append("\"".join(splitLine))        
        else:
            newLines.append(line)
            
    newXMLString = "\n".join(newLines)
    h5IDF[...]=newXMLString#.encode(encoding='ascii')
    h5obj.close()
    #
    # Finish up
    #
    stepTime = time.time() - time_0
    totTime += stepTime
    print(f'    Time to complete step: {stepTime:.4f} sec. Total time to execute: {totTime:.4f}')
    return outFileName

def superID(nativeID,xdim,ydim):
    import numpy as np
    #accepts a numpy array of native ID from standard SNAP nexus file and returns a numpy array with 
    # super pixel ID according to 
    #provided dimensions xdim and ydim of the super pixel. xdim and ydim shall be multiples of 2
    #

    Nx = 256 #native number of horizontal pixels 
    Ny = 256 #native number of vertical pixels 
    NNat = Nx*Ny #native number of pixels per panel
    
    firstPix = (nativeID // NNat)*NNat
    redID = nativeID % NNat #reduced ID beginning at zero in each panel
    
    (i,j) = divmod(redID,Ny) #native (reduced) coordinates on pixel face
    superi = divmod(i,xdim)[0]
    superj = divmod(j,ydim)[0]

    #some basics of the super panel   
    superNx = Nx/xdim #32 running from 0 to 31
    superNy = Ny/ydim
    superN = superNx*superNy

    superFirstPix = (firstPix/NNat)*superN
    
    super = superi*superNy+superj+superFirstPix
    super = super.astype('int')

#     print('native ids: ',nativeID)
#     print('first pixels: ',firstPix)
#     print('super first pixels: ',superFirstPix)
    
    
    return super