# python_export_cae.py
# 作用：根据 mat-cae.mat 指定的 (T,Phe)，读取 mat-T-Phe.mat 建模并保存 CAE，然后退出（不求解）

import sys
import time
import numpy as np
import scipy.io as sio

from abaqus import mdb
from abaqusConstants import *
from regionToolset import Region


# ===================== 读取导出指令 =====================
mat_cae = sio.loadmat('mat-cae.mat')['FEA'][0][0]

def get_int_field(field_name, default=0):
    if field_name in mat_cae.dtype.names:
        return int(mat_cae[field_name][0][0])
    return default

export_mode = get_int_field('EXPORT_MODE', 0)
if export_mode != 1:
    print('EXPORT_MODE != 1，本脚本仅用于导出 CAE。直接退出。')
    sys.exit(0)

T_target   = get_int_field('EXPORT_T', 0)
Phe_target = get_int_field('EXPORT_PHE', 0)

if T_target <= 0 or Phe_target <= 0:
    print('EXPORT_T / EXPORT_PHE 无效：', T_target, Phe_target)
    sys.exit(1)

General = str(T_target)
name_mat = 'mat-' + General + '-' + str(Phe_target) + '.mat'
print('>>> Export target:', name_mat)

# ===================== 读取目标个体参数 =====================
try:
    mat_model = sio.loadmat(name_mat)
except Exception as e:
    print('ERROR: cannot load', name_mat, e)
    sys.exit(1)

# ============================================================
# ✅ 关键：把你 python_CAE_Press_M.py 里“建模部分”复制到这里
# 要求：
#   1) 只做建模：创建 model/part/assembly/材料/截面/网格/BC/Step/Load 等
#   2) 绝对不要出现：submit(), waitForCompletion()
#   3) 不需要读多个人（这里只导出一个 mat_model）
#
# 通常你原脚本里是：
#   for i in range(1, num_Phe+1):
#       mat_model = sio.loadmat(name_mat)
#       ... build ...
#
# 你现在只需要把循环体里的“build代码”拿出来，直接对 mat_model 执行一次
# ============================================================

# ----------------- 在这里粘贴/调用你的建模代码 -----------------
# 例如（示意）：
# build_from_mat(mat_model, General)    # 如果你愿意以后我也可以帮你封装成函数
#
# 你如果暂时不想封装，直接把原来的建模段落复制进来即可。
# ----------------- 这里结束 -----------------


Phe        = str( mat_model['Phe'][0][0] )
NAME_JOB   = str( 'Job-'+ General + '-' + Phe )
NAME_MODEL = str( 'M-' + Phe )

sketchdata = mat_model['sketch'][0][0]
feadata = mat_model['FEA'][0][0]

P1_AZ = float(sketchdata['Assemble_P1Z'])
P2_AZ = float(sketchdata['Assemble_P2Z'])
TAR_LINE = sketchdata['tarline']
P1_DIV = sketchdata['div']
P1_CTR = sketchdata['center']
P1_FREE = sketchdata['find_free']
P1_BOND = sketchdata['find_bond']
#NUM_LINE = sketchdata['num_line']
NUM_LINE = int(np.array(sketchdata['num_line']).squeeze())
RING_P = sketchdata['ring']
BOND_U = sketchdata['U']

BASE_R = float(mat_model['tar'][0][0]['r'])
P1_T   = float(feadata['t'])
LOAD_PRESS_ALL = float(feadata['load_pressall'])
LOAD_PRESS_B1   = float(feadata['load_pressbond'])
LOAD_DISTURB   = float(feadata['load_disturb'])
P1_MESH  = float(feadata['mesh'][0])
STEP_TIME = feadata['step_time']
STEP_STABILIZE = feadata['step_stabilize']
PENALTY_1 = float(feadata['penalty_1'])

# ABAQUS START

## Create Model: M-i
mdb.Model(modelType=STANDARD_EXPLICIT, name=(NAME_MODEL))
M = mdb.models[NAME_MODEL]
#--------------------------------------------------------Part
## Create PART 1 : Target 2D
ske1 = 'P1-Line-1'
M.ConstrainedSketch(name = ske1, sheetSize = 200.0)
for j in range(1,5):
    PP = sketchdata['line'+str(j)][0][0]
    for k in range(np.size( PP, 0)-1):
            M.sketches[ske1].Line(
                point1=(PP[k  ][0], PP[k  ][1] ), 
                point2=(PP[k+1][0], PP[k+1][1] ) )
M.Part(dimensionality=THREE_D, name='PART-1', type=DEFORMABLE_BODY)
M.parts['PART-1'].BaseShell(sketch=M.sketches[ske1])
# IF Multilines :
if NUM_LINE>1:
    for kl in range(1,NUM_LINE) :
        ske1 = 'P1-Line-' + str(kl+1)
        M.ConstrainedSketch(name = ske1, sheetSize = 200.0)
        for j in range(1,5):
            PP = sketchdata['line'+str(j)][kl][0]
            for k in range(np.size(PP,0)-1):
                    M.sketches[ske1].Line(
                        point1=(PP[k  ][0], PP[k  ][1] ), 
                        point2=(PP[k+1][0], PP[k+1][1] ) )
        M.parts['PART-1'].Shell(sketch=M.sketches[ske1])

## Create BASE : Cylinder
ske2 = 'P2'
M.ConstrainedSketch(name= ske2, sheetSize=400.0)
M.sketches[ske2].ConstructionLine(
    point1=(0.0, -200.0), point2=(0.0, 200.0))
M.sketches[ske2].Line(
    point1=(BASE_R, 250.0), point2=(BASE_R, -25.0))
M.Part(dimensionality=THREE_D, name='BASE', type=ANALYTIC_RIGID_SURFACE)
M.parts['BASE'].AnalyticRigidSurfRevolve(
    sketch=M.sketches[ske2])
RP_P2 = M.parts['BASE'].ReferencePoint(point=(0.0, 200.0, 0.0))

## Create TAR : CenterLine
M.Part(dimensionality=THREE_D, name='TAR', type=DEFORMABLE_BODY)
M.parts['TAR'].ReferencePoint(point=(0.0, 0.0, 0.0))
for kl in range( 0, NUM_LINE):
    p = []
    PP = TAR_LINE[kl][0]
    for j in range(0, np.size(PP,0)):
        p.append(( PP[j][0], PP[j][1], PP[j][2]))
    M.parts['TAR'].WireSpline(
        mergeType=IMPRINT, meshable=ON, points=(p), smoothClosedSpline=ON)

## Create Ring : For Bonding  ( Total : Num of Bonding )
for kl in range(0,NUM_LINE) :
    PP = RING_P[kl][0]
    for j in range(0,np.size(PP,0)):
        Ring_name = 'Ring-' + str(kl+1) + '-' + str(j+1)
        ske3 = 'Ring-' + str(kl+1) + '-' + str(j)
        M.ConstrainedSketch(name=(ske3), sheetSize=400.0)
        M.sketches[(ske3)].ConstructionLine(
            point1=(0.0, -200.0), point2=(0.0, 200.0))
        M.sketches[(ske3)].Line(
            point1=( BASE_R, PP[j][0]+10), point2=(BASE_R, PP[j][0]-10))
        M.Part(dimensionality=THREE_D, name=Ring_name, type=ANALYTIC_RIGID_SURFACE)
        M.parts[Ring_name].AnalyticRigidSurfRevolve(sketch=M.sketches[ske3] ) 
        locals()['RP-R'+ str(kl+1) + '-' + str(j)] = M.parts[Ring_name].ReferencePoint(
            point=(0.0, float(PP[j][0]), 0.0))

## Divide : P1
sked1 = 'P1_DIV'
M.ConstrainedSketch(name=sked1, sheetSize=200)
M.parts['PART-1'].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES, 
    sketch=M.sketches[sked1])
for kl in range(0,NUM_LINE):
    PP = P1_CTR[kl][0]
    for j in range( np.size(PP,0) - 1 ):
        M.sketches[sked1].Line(
            point1=(PP[j]  [0], PP[j]  [1]), 
            point2=(PP[j+1][0], PP[j+1][1]) )
    PP = P1_DIV[kl][0]
    for j in range(np.size(PP,0)):
        M.sketches[sked1].Line(
            point1=(PP[j][0], PP[j][1]), 
            point2=(PP[j][2], PP[j][3]) )

M.parts['PART-1'].PartitionFaceBySketch(
    faces=M.parts['PART-1'].faces, 
    sketch=M.sketches[sked1], 
    sketchOrientation=LEFT)

## ----------------------------------Material
## Create Material : PI
M.Material(name='PI')
M.materials['PI'].Elastic(table=((2500.0, 0.3), ))
## Create Section : PI-PART-1
M.HomogeneousShellSection(
        idealization=NO_IDEALIZATION, integrationRule=SIMPSON, 
        material='PI', name='PI-PART-1', 
        nodalThicknessField='', numIntPts=5, 
        poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
        thickness=P1_T, 
        thicknessField='', thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
#Designate Section : P1
M.parts['PART-1'].SectionAssignment(
        offset=0.0, offsetField='', 
        offsetType=BOTTOM_SURFACE, 
        region=Region(faces=M.parts['PART-1'].faces), 
        sectionName='PI-PART-1', 
        thicknessAssignment=FROM_SECTION)
#Designate Section : Line  ( TAR_LINE )
M.CircularProfile(name='Profile-1', r=1.0)
M.BeamSection(
    consistentMassMatrix=False, 
    integration=DURING_ANALYSIS, 
    material='PI', 
    name='Line', 
    poissonRatio=0.0, 
    profile='Profile-1', 
    temperatureVar=LINEAR)
M.parts['TAR'].SectionAssignment(
    offset=0.0, 
    offsetField='', 
    offsetType=MIDDLE_SURFACE, 
    region=Region(edges=M.parts['TAR'].edges),
    sectionName='Line', 
    thicknessAssignment=FROM_SECTION ) 
M.parts['TAR'].assignBeamSectionOrientation(
    method=N1_COSINES, 
    n1=(0.0, 0.0, -1.0), 
    region=Region(edges=M.parts['TAR'].edges))

#--------------------------------------------Assemble
#Create Instance : P1
#Translate : Z + P1_AZ
M.rootAssembly.Instance(
    dependent=ON, 
    name='PART-1-1', 
    part=M.parts['PART-1'])
M.rootAssembly.translate(
        instanceList=('PART-1-1', ), 
        vector=(0.0, 0.0, P1_AZ))

#Create Instance : P2
#Rotation : URZ = -90.0
M.rootAssembly.Instance(
    dependent=ON, 
    name='BASE-1', 
    part=M.parts['BASE'])
M.rootAssembly.rotate(
    angle = -90.0, 
    axisDirection=(0.0, 0.0, 10.0), 
    axisPoint    =(0.0, 0.0, 0.0), 
    instanceList=('BASE-1', ))
M.rootAssembly.translate(instanceList=('BASE-1', ), vector=(
0.0, 0.0, -P2_AZ))

#Create Instance : P3 ( TAR_LINE )
M.rootAssembly.Instance(
    dependent=ON, 
    name='TAR-1', 
    part=M.parts['TAR'])

#Create Instance : Ring
#Rotation : URZ = -90.0
for kl in range(0,NUM_LINE) :
    PP = RING_P[kl][0]
    for j in range(0,np.size(PP,0)):
        Ring_name = 'Ring-' + str(kl+1) + '-' + str(j+1)
        Ring_A    = Ring_name + '-1'
        M.rootAssembly.Instance(
            dependent=ON, 
            name = Ring_A, 
            part=M.parts[Ring_name])
        M.rootAssembly.rotate(
            angle = -90.0, 
            axisDirection=(0.0, 0.0, 10.0), 
            axisPoint    =(0.0, 0.0, 0.0), 
            instanceList=( Ring_A, ))   

#------------------------------------------------Step
#Step-1 : Load   : Disturbance
#Step-2 : Release: Deformation
#Step-3 : Unload : Disturbance
M.StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, 
                                    initialInc=0.01, 
                                    maxNumInc=400, 
                                    minInc=1e-7, 
                                    name='Step-1', nlgeom=ON, previous='Initial',
                                    stabilizationMagnitude= STEP_STABILIZE[0][0], 
                                    stabilizationMethod=DAMPING_FACTOR, timePeriod=float(STEP_TIME[0][0]))
M.StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, 
                                    initialInc=0.001, 
                                    maxNumInc=400, 
                                    minInc = 1e-7, 
                                    name='Step-2', previous='Step-1', 
                                    stabilizationMagnitude= STEP_STABILIZE[1][0],
                                    stabilizationMethod=DAMPING_FACTOR, timePeriod=float(STEP_TIME[1][0]))
M.StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, 
                                    initialInc=0.001, 
                                    maxNumInc=400, 
                                    minInc = 1e-7, 
                                    name='Step-3', previous='Step-2', 
                                    stabilizationMagnitude= STEP_STABILIZE[2][0],
                                    stabilizationMethod=DAMPING_FACTOR, timePeriod=float(STEP_TIME[2][0]))
M.StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, 
                                    initialInc=0.001, 
                                    maxNumInc=400, 
                                    minInc = 1e-7, 
                                    name='Step-4', previous='Step-3', 
                                    stabilizationMagnitude= STEP_STABILIZE[3][0],
                                    stabilizationMethod=DAMPING_FACTOR, timePeriod=float(STEP_TIME[3][0]))
if 1==0:
    M.StaticStep(adaptiveDampingRatio=None, continueDampingFactors=False, 
                                    initialInc=0.5, 
                                    maxInc=5.0, 
                                    maxNumInc=400, 
                                    minInc=1e-7, 
                                    name='Step-5', previous='Step-4', 
                                    stabilizationMagnitude= STEP_STABILIZE[4][0], 
                                    stabilizationMethod=DAMPING_FACTOR, timePeriod=float(STEP_TIME[4][0]))
#Adjust : OutputFile
M.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S','PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'CFORCE','CSTATUS', 'COORD'),
        numIntervals=5, timeMarks=OFF,)   #Only 5 result per step
        #frequency=LAST_INCREMENT)        #Only Last increment
        #)                                #None: All increment

# Adjust : Restart
M.steps['Step-1'].Restart(frequency=0, numberIntervals=1, 
    overlay=OFF, timeMarks=OFF)
M.steps['Step-2'].Restart(frequency=0, numberIntervals=1, 
    overlay=OFF, timeMarks=OFF)
M.steps['Step-3'].Restart(frequency=0, numberIntervals=1, 
    overlay=OFF, timeMarks=OFF)
M.steps['Step-4'].Restart(frequency=0, numberIntervals=1, 
    overlay=OFF, timeMarks=OFF)

#Adjust : IA
M.steps['Step-1'].control.setValues(
    allowPropagation=OFF,
    resetDefaultValues=OFF, 
    timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0, 12.0, 10.0, 6.0, 3.0, 50.0))

#-------------------------Surface and Set
#Surface
##############  side1Face --> + vector   side2Face --> - vector 
    
# P1 - BOTTOM
M.rootAssembly.Surface(
    name='P1-BOTTOM',
    side2Faces=M.rootAssembly.instances['PART-1-1'].faces)

# P1 - FREE
sur = M.rootAssembly.instances['PART-1-1'].faces
freeface = sur.findAt(((P1_FREE[0][0][0],P1_FREE[0][1][0], P1_AZ),))
for kl in range(0,NUM_LINE):
    newface = sur.findAt(((P1_FREE[0][0][kl],P1_FREE[0][1][kl], P1_AZ),),
                        ((P1_FREE[1][0][kl],P1_FREE[1][1][kl], P1_AZ),),
                        ((P1_FREE[2][0][kl],P1_FREE[2][1][kl], P1_AZ),),
                        ((P1_FREE[3][0][kl],P1_FREE[3][1][kl], P1_AZ),) )
    freeface = freeface + newface
M.rootAssembly.Surface(
    name='FREE', 
    side2Faces=freeface)


# BASE - OUT 
M.rootAssembly.Surface(
    name='BASE-OUTER', 
    side1Faces=M.rootAssembly.instances['BASE-1'].faces.findAt(
        (( 0, 0, BASE_R - P2_AZ ),) ))

# P1 - Bond
for kl in range(0,NUM_LINE):
    PP = P1_BOND[kl][0]
    for j in range(0,np.size(PP,0)):
        M.rootAssembly.Surface(
            name='P1-Bond-' + str(kl+1) + '-' + str(j+1),
            side2Faces=M.rootAssembly.instances['PART-1-1'].faces.findAt( 
                ((PP[j][0], PP[j][1], P1_AZ),) ) )
        
# Ring - Surface
for kl in range(0,NUM_LINE):
    PP = RING_P[kl][0]
    for j in range(0,np.size(PP,0)):
        Ring_name = 'Ring-' + str(kl+1) + '-' + str(j+1)
        Ring_A    = Ring_name + '-1'
        Sur_name = 'R-' + str(kl+1) + '-' + str(j+1) + '-Bond'
        M.rootAssembly.Surface(
            name= Sur_name, 
            side1Faces=M.rootAssembly.instances[Ring_A].faces.findAt((( float(PP[j][0]), 0,BASE_R),) ))    


#-----------------------Interaction
#Create : Attribute : Nof
M.ContactProperty('NoF')
M.interactionProperties['NoF'].TangentialBehavior(formulation=FRICTIONLESS)
M.interactionProperties['NoF'].NormalBehavior(
        allowSeparation=ON, constraintEnforcementMethod=DEFAULT, pressureOverclosure=HARD)
#Create : Attribute : Rough - FIX
M.ContactProperty('Rough')
M.interactionProperties['Rough'].TangentialBehavior(
    formulation=ROUGH)
M.interactionProperties['Rough'].NormalBehavior(
    allowSeparation=OFF, constraintEnforcementMethod=DEFAULT, 
    pressureOverclosure=HARD)
#Create : Attribute : Rough - FREE
M.ContactProperty('Rough-Penalty')
M.interactionProperties['Rough-Penalty'].TangentialBehavior(
    dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
    formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
    pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
    table=((0.5, ), ), temperatureDependency=OFF)
M.interactionProperties['Rough-Penalty'].NormalBehavior(
    allowSeparation=OFF, constraintEnforcementMethod=DEFAULT, 
    pressureOverclosure=HARD)
#Create : Attribute : Penalty-0.1
M.ContactProperty('Penalty-1')
M.interactionProperties['Penalty-1'].TangentialBehavior(
    dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
    formulation=PENALTY, 
    fraction=0.005, maximumElasticSlip=FRACTION, pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
    table=((0.1, ), ), 
    temperatureDependency=OFF)
M.interactionProperties['Penalty-1'].NormalBehavior(
        allowSeparation=ON, constraintEnforcementMethod=DEFAULT, pressureOverclosure=HARD)

#Create:Interaction
M.SurfaceToSurfaceContactStd(
        adjustMethod=NONE, 
        clearanceRegion=None, 
        createStepName='Initial', 
        datumAxis=None, 
        initialClearance=OMIT, 
        enforcement=NODE_TO_SURFACE, 
        #interactionProperty='NoF', 
        interactionProperty='NoF', 
        main=M.rootAssembly.surfaces['BASE-OUTER'], 
        name='P1-P2', 
        secondary =M.rootAssembly.surfaces['FREE'], 
        sliding=FINITE, 
        thickness=ON)

for kl in range(0,NUM_LINE):
    PP = RING_P[kl][0]
    for j in range(0,np.size(PP,0)):
        Bond_name = 'Bond-' + str(kl+1) + str(j+1)
        Ring_surname = 'R-' + str(kl+1) + '-' + str(j+1) + '-Bond'
        Part_surname = 'P1-Bond-' + str(kl+1) + '-' + str(j+1) 
        M.SurfaceToSurfaceContactStd(
            adjustMethod=NONE, 
            clearanceRegion=None, 
            createStepName='Initial', 
            datumAxis=None, 
            initialClearance=OMIT, 
            enforcement=NODE_TO_SURFACE, 
            interactionProperty='Rough-Penalty', 
            main=M.rootAssembly.surfaces[Ring_surname], #这边把master修改成了main 下面的slave修改成了secondary
            name=Bond_name, 
            secondary =M.rootAssembly.surfaces[Part_surname], 
            sliding=FINITE, 
            thickness=ON)       
        M.interactions[Bond_name].setValuesInStep(interactionProperty='Rough', stepName='Step-2')


#---------------------Load
# Create : AMP
M.SmoothStepAmplitude(
    data=((0.0, 0.0), (0.5, 0.0), (1, 1.0)), 
    name='Amp-Bond', 
    timeSpan=STEP)
M.SmoothStepAmplitude(
    data=((0.0, 0.0), (0.5, 1), (1,1)), 
    name='Amp-All', 
    timeSpan=STEP)
#Load : Disturbance
M.Pressure(
    amplitude='Amp-All', createStepName='Step-2', distributionType=UNIFORM, field='', 
    magnitude = LOAD_DISTURB, 
    name='Load-Disturbance', 
    region=M.rootAssembly.surfaces['FREE'])
M.loads['Load-Disturbance'].deactivate('Step-3')

#Load : Press-ALL
M.Pressure(
        amplitude='Amp-All', createStepName='Step-1', distributionType=UNIFORM, field='', 
        magnitude = -LOAD_PRESS_ALL, 
        name='Load-Press-All', 
        region=M.rootAssembly.surfaces['P1-BOTTOM'])
M.loads['Load-Press-All'].deactivate('Step-2')
#Load : Press-Bond
for kl in range(0,NUM_LINE):
    PP = RING_P[kl][0]
    for j in range(0,np.size(PP,0)):
        Sur_name = 'P1-Bond-' + str(kl+1) + '-' + str(j+1)
        Load_name = 'Load-Press-Bond-'+ str(kl+1) + '-' + str(j+1)
        M.Pressure(
            amplitude='Amp-Bond', createStepName='Step-1', distributionType=UNIFORM, field='', 
            magnitude = -LOAD_PRESS_B1, 
            name=Load_name, 
            region=M.rootAssembly.surfaces[Sur_name])
        M.loads[Load_name].deactivate('Step-2')
    
#Bondary Conditions : Rigid Surface
M.DisplacementBC(
    amplitude=UNSET, 
    createStepName='Initial', 
    distributionType=UNIFORM, 
    fieldName='', localCsys=None, 
    name='BC-Surface', 
    region=Region(referencePoints=(
    M.rootAssembly.instances['BASE-1'].referencePoints[RP_P2.id], )), 
    u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
M.boundaryConditions['BC-Surface'].setValuesInStep(stepName='Step-1', 
    u3= P2_AZ )

#Bondary Conditions : Bondings
for kl in range(0,NUM_LINE):
    PP = RING_P[kl][0]
    for j in range(0,np.size(PP,0)):
        Ring_A   = 'Ring-' + str(kl+1) + '-' + str(j+1) + '-1'
        BC_name  = 'BC-Bond-' + str(kl+1) + '-' + str(j+1)
        M.DisplacementBC(
            amplitude=UNSET, 
            createStepName='Step-1', 
            distributionType=UNIFORM, 
            fieldName='', 
            fixed=OFF, 
            localCsys=None, 
            name=BC_name, 
            region=Region(referencePoints=(M.rootAssembly.instances[Ring_A].referencePoints[(locals()['RP-R'+ str(kl+1) + '-' + str(j)]).id], )), 
            u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
        M.boundaryConditions[BC_name].setValuesInStep(
            stepName='Step-3', 
            u1  = BOND_U[kl][0][j][0], 
            u2  = BOND_U[kl][0][j][1], 
            u3  = BOND_U[kl][0][j][2],
            ur1 = BOND_U[kl][0][j][3], 
            ur2 = BOND_U[kl][0][j][4], 
            ur3 = BOND_U[kl][0][j][5] )
        M.boundaryConditions[BC_name].setValuesInStep(
            stepName='Step-4',  
            u1  = BOND_U[kl][1][j][0], 
            u2  = BOND_U[kl][1][j][1], 
            u3  = BOND_U[kl][1][j][2],
            ur1 = BOND_U[kl][1][j][3], 
            ur2 = BOND_U[kl][1][j][4], 
            ur3 = BOND_U[kl][1][j][5] )


#------------------------------------Mesh
#Mesh : P1
M.parts['PART-1'].setMeshControls(
        elemShape=QUAD, 
        regions=M.parts['PART-1'].faces, 
        technique=FREE )
    #technique=STRUCTURED) 
M.parts['PART-1'].seedPart(
        deviationFactor=0.1, 
        minSizeFactor=0.1, 
        size=P1_MESH)
#M.parts['PART-1'].setElementType(
#    elemTypes=(
#        ElemType(elemCode=S4, elemLibrary=STANDARD, secondOrderAccuracy=OFF), 
#        ElemType(elemCode=S3, elemLibrary=STANDARD)), 
#        regions=(M.parts['PART-1'].faces, ))

M.parts['PART-1'].generateMesh()

#Mesh : P3
M.parts['TAR'].seedPart(
    deviationFactor=0.1, 
    minSizeFactor=0.1, 
    size=0.5)
M.parts['TAR'].generateMesh()

M.rootAssembly.regenerate()

# ===================== 保存 CAE 并退出 =====================
cae_name = 'CAE-T' + str(T_target) + '-Phe' + str(Phe_target)
try:
    mdb.saveAs(cae_name)
    print('>>> Saved:', cae_name + '.cae')
except Exception as e:
    print('ERROR: saveAs failed:', e)
    sys.exit(1)

time.sleep(1)
sys.exit(0)
