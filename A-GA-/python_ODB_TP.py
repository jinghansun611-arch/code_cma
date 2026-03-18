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

matname='mat-0.mat'
matdata = sio.loadmat(matname)['sketch'][0][0]
pathpoint = matdata['path']
odbname =str('Job-0.odb')
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
    #Create XY_data : U1 (X)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 'COOR1'), )
    pth = session.paths['Path-P1']
    session.XYDataFromPath(name='XYData-X', path=pth, includeIntersections=True, 
        projectOntoMesh=True, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE, 
        removeDuplicateXYPairs=True, includeAllElements=False)
    #Create XY_data : U2 (Y)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 'COOR2'), )
    pth = session.paths['Path-P1']
    session.XYDataFromPath(name='XYData-Y', path=pth, includeIntersections=True, 
        projectOntoMesh=True, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE, 
        removeDuplicateXYPairs=True, includeAllElements=False)
    #Create XY_data : U3 (Z)
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
        variableLabel='COORD', outputPosition=NODAL, refinement=(COMPONENT, 'COOR3'), )
    pth = session.paths['Path-P1']
    session.XYDataFromPath(name='XYData-Z', path=pth, includeIntersections=True, 
        projectOntoMesh=True, pathStyle=PATH_POINTS, numIntervals=10, 
        projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE, 
        removeDuplicateXYPairs=True, includeAllElements=False)
    #Output data
    x0 = session.xyDataObjects['XYData-X']
    session.writeXYReport(fileName="TXT-RESULT-X-" + str(kl+1) +".txt", xyData=(x0, ))
    x0 = session.xyDataObjects['XYData-Y']
    session.writeXYReport(fileName="TXT-RESULT-Y-" + str(kl+1) +".txt", xyData=(x0, ))
    x0 = session.xyDataObjects['XYData-Z']
    session.writeXYReport(fileName="TXT-RESULT-Z-" + str(kl+1) +".txt", xyData=(x0, ))
#Close odb
session.odbs[odbname].close()

#Write log information
meg='get data SUCCESSFULLY\n'
fmeg.write('------------%s:  %s'%(odbname,meg))
print(meg)
fmeg.close()