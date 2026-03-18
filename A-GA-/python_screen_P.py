from abaqus import *
from abaqusConstants import *

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