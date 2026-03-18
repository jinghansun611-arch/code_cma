# -*- coding: mbcs -*
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from textRepr import *
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import numpy as np
import scipy.io as sio
import sys
import os
#Open log file
#Open log file
fmeg=open('pylog.txt','w')
fmeg.truncate(0)
#read data
ptdata=sio.loadmat('mat-odb.mat')
PHE= str(ptdata['Phe'][0][0])
T  = str(ptdata['T'][0][0])
matname='mat-'+T+'-'+PHE+'.mat'
#CIR_PATH POINT 1: (X,Y,Z)=[0] [1] [2]  
#CIR_PATH POINT 2: (X,Y,Z)=[3] [4] [5]
matdata = sio.loadmat(matname)['sketch'][0][0]
pathpoint = matdata['path']
odbname=str('Job-'+T+'-'+PHE+'.odb')
#Odb name

#Open odb
o1 = session.openOdb(
    name=odbname)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#Create Path_P1
for kl in range(0,np.size(pathpoint)) :
    #Clear history
    file = open("TXT-RESULT-X-" + str(kl+1) +".txt", 'w').close()
    file = open("TXT-RESULT-Z-" + str(kl+1) +".txt", 'w').close()
    file = open("TXT-RESULT-Y-" + str(kl+1) +".txt", 'w').close()
    p = []
    PP = pathpoint[kl][0]
    for i in range(0, np.size(PP,0)):
        p.append((PP[i][0],PP[i][1],PP[i][2]))
    session.Path(name='Path-P1', type=POINT_LIST, expression=(p))

# Set background
session.graphicsOptions.setValues(
    backgroundStyle=SOLID, 
    backgroundColor='#FFFFFF')

session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    renderShellThickness=ON)

session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=FREE)
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    renderBeamProfiles=ON)
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    beamScaleFactor=0.25)
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    numSweepSegmentsArs=200)
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#F1F1F1')
cmap = session.viewports['Viewport: 1'].colorMappings['Default']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()