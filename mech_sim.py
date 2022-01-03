import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import csv

E11_b=True
E22_b=True
E33_b=False
G12_b=True
G23_b=False
G13_b=False
# properties:
#set 3
# E_f=340
# E_m=2.9
# Poisson_f=0.23
# Poisson_m=0.34
# rho_f=3.9E-9
# rho_m=1.2E-9
# Cond_f=30
# Cond_m=0.3
#set 2
# E_f=66.9
# E_m=2.9
# Poisson_f=0.29
# Poisson_m=0.34
# rho_f=2.52E-9
# rho_m=1.2E-9
# Cond_f=1.4
# Cond_m=0.3
#set 1
E_f=110
E_m=2.5
Poisson_f=0.31
Poisson_m=0.34
rho_f=4.45E-9
rho_m=1.42E-9
Cond_f=7.1
Cond_m=0.12
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.ModelFromInputFile(name='Mesh_model_2', 
    inputFileName='C:/temp/Mesh_model_2.inp')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Mesh_model_2'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
from abaqus import *
from abaqusConstants import *
import __main__
mdb.Job(name='Job-1', model='Mesh_model_2', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1, numGPUs=0)
mdb.jobs['Job-1'].writeInput(consistencyChecking=OFF)
# change inp file and import again

#%%
import numpy as np
import os


def search_string_in_file(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append((line_number, line.rstrip()))
    return list_of_results

def delete_multiple_lines(original_file, line_numbers):
    """In a file, delete the lines at line number in given list"""
    is_skipped = False
    counter = 0
    dummy_file = original_file + '.bak'
    with open(original_file, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        for line in read_obj:
            if counter not in line_numbers:
                write_obj.write(line)
            else:
                is_skipped = True
            counter += 1
    if is_skipped:
        os.remove(original_file)
        os.rename(dummy_file, original_file)
    else:
        os.remove(dummy_file)

def copy_multiple_lines(new_file,original_file, line_numbers):
    
    # print('line numbers',line_numbers)
    counter = 0
  
    with open(original_file, 'r') as read_obj, open(new_file, 'w+') as write_obj:
        for line in read_obj:
            
            if counter in line_numbers:
                # print('counter',counter)
                write_obj.write(line)
            
            counter+=1

#  copy paste 

file_name='C:\\temp\Job-1.inp'
string_to_search='VOLUME'
datatocpy=search_string_in_file(file_name, string_to_search)
# print(datatocpy[-1][0]+1)

new_file= "C:\\temp\part_set.inp"
copy_multiple_lines(new_file,file_name, np.arange(datatocpy[0][0]-1,datatocpy[-1][0]+1))

#clean assy data from file
file_name= "C:\\temp\part_set.inp"
outfile = "C:\\temp\cleaned_file.inp"

delete_list = [", instance=PART-1-1"]
fin = open(file_name)
fout = open(outfile, "w+")
for line in fin:
    for word in delete_list:
       line = line.replace(word, "")
    fout.write(line)
fin.close()
fout.close()


file_name='C:\\temp\Job-1.inp'
in_file= 'C:\\temp\cleaned_file.inp'
string_to_search='*End Part'
row_to_insert=search_string_in_file(file_name, string_to_search)

exp = row_to_insert[0][0]-1
print('inser in line',exp)
with open (file_name,'r') as b:
    lines = b.readlines()
with open(file_name,'w') as b:
    for i,line in enumerate(lines):
        if i == exp:
            a=open(in_file,'r')
            for l in a:
                b.write(l)
            
        b.write(line)


with open('C:\\temp\section_data.txt','r') as s:
    data=s.readlines()
    fibs_s=int(data[0])
    fibs=int(data[1])
    nvol=int(data[2])
    model_num=data[6]
    model_num=model_num.replace('\n','')
    dir='C:/temp/model'+str(model_num)
    with open(dir+'/sim_data.txt','w') as r:
        
        for i in data:
            r.write(i+'\n')
#%%
print('started simulation in model num ',model_num)

mdb.ModelFromInputFile(name='Job-1', 
    inputFileName='C:/temp/Job-1.inp')


import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Job-1'].parts['PART-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['Job-1'].Material(name='Material-1')
mdb.models['Job-1'].materials['Material-1'].Density(table=((rho_f, ), ))
mdb.models['Job-1'].materials['Material-1'].Elastic(table=((E_f, Poisson_f), ))
mdb.models['Job-1'].materials['Material-1'].Conductivity(table=((Cond_f, ), ))
mdb.models['Job-1'].Material(name='Material-2')
mdb.models['Job-1'].materials['Material-2'].Density(table=((rho_m, ), ))
mdb.models['Job-1'].materials['Material-2'].Elastic(table=((E_m, Poisson_m), ))
mdb.models['Job-1'].materials['Material-2'].Conductivity(table=((Cond_m, ), ))
mdb.models['Job-1'].HomogeneousSolidSection(name='Section-1', 
    material='Material-1', thickness=None)
mdb.models['Job-1'].HomogeneousSolidSection(name='Section-2', 
    material='Material-2', thickness=None)

p = mdb.models['Job-1'].parts['PART-1']

if(fibs_s!=1):
    for i in range (1,fibs_s):
        region = p.sets['VOLUME'+str(i)]
        p = mdb.models['Job-1'].parts['PART-1']
        p.SectionAssignment(region=region, sectionName='Section-2')    
for i in range (fibs_s,fibs):

    region = p.sets['VOLUME'+str(i)]
    p = mdb.models['Job-1'].parts['PART-1']
    p.SectionAssignment(region=region, sectionName='Section-1')

for i in range (fibs,nvol+1):
    region = p.sets['VOLUME'+str(i)]
    p = mdb.models['Job-1'].parts['PART-1']
    p.SectionAssignment(region=region, sectionName='Section-2')


#copy mesh data

#geometry snap
dir='C:/temp/model'+str(model_num)
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Material']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Material']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.printToFile(fileName=dir+'/snap', format=PNG, 
    canvasObjects=(session.viewports['Viewport: 1'], ))
#easypbc initial

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import sys


sys.path.insert(15, r'c:/temp/abaqus_plugins/EasyPBC V.1.3')
import easypbc
easypbc.feasypbc(part='Job-1', inst='PART-1-1', meshsens=1E-05, CPU=4, 
        E11=E11_b, E22=E22_b, E33=E33_b, G12=G12_b, G13=G13_b, G23=G23_b, 
        onlyPBC=False)


with open('C:\\temp\Job-1_elastic_properties.txt','r') as s:
    data=s.readlines()

    dir='C:/temp/model'+str(model_num)
    with open(dir+'/sim_data.txt','a') as r:
        r.write('MECHANICAL\n')
        for i in data:
            r.write(i+'\n')

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import csv

from abaqus import *
from abaqusConstants import *
import __main__
mdb.Job(name='Job-1', model='Mesh_model_2', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=1, numGPUs=0)
mdb.jobs['Job-1'].writeInput(consistencyChecking=OFF)
# change inp file and import again

#%%
import numpy as np
import os


def search_string_in_file(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append((line_number, line.rstrip()))
    return list_of_results

def delete_multiple_lines(original_file, line_numbers):
    """In a file, delete the lines at line number in given list"""
    is_skipped = False
    counter = 0
    dummy_file = original_file + '.bak'
    with open(original_file, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        for line in read_obj:
            if counter not in line_numbers:
                write_obj.write(line)
            else:
                is_skipped = True
            counter += 1
    if is_skipped:
        os.remove(original_file)
        os.rename(dummy_file, original_file)
    else:
        os.remove(dummy_file)

def copy_multiple_lines(new_file,original_file, line_numbers):
    
    # print('line numbers',line_numbers)
    counter = 0
  
    with open(original_file, 'r') as read_obj, open(new_file, 'w+') as write_obj:
        for line in read_obj:
            
            if counter in line_numbers:
                # print('counter',counter)
                write_obj.write(line)
            
            counter+=1

#  copy paste 

file_name='C:\\temp\Job-1.inp'
string_to_search='VOLUME'
datatocpy=search_string_in_file(file_name, string_to_search)
# print(datatocpy[-1][0]+1)

new_file= "C:\\temp\part_set.inp"
copy_multiple_lines(new_file,file_name, np.arange(datatocpy[0][0]-1,datatocpy[-1][0]+1))

#clean assy data from file
file_name= "C:\\temp\part_set.inp"
outfile = "C:\\temp\cleaned_file.inp"

delete_list = [", instance=PART-1-1"]
fin = open(file_name)
fout = open(outfile, "w+")
for line in fin:
    for word in delete_list:
       line = line.replace(word, "")
    fout.write(line)
fin.close()
fout.close()


file_name='C:\\temp\Job-1.inp'
in_file= 'C:\\temp\cleaned_file.inp'
string_to_search='*End Part'
row_to_insert=search_string_in_file(file_name, string_to_search)

exp = row_to_insert[0][0]-1
print('inser in line',exp)
with open (file_name,'r') as b:
    lines = b.readlines()
with open(file_name,'w') as b:
    for i,line in enumerate(lines):
        if i == exp:
            a=open(in_file,'r')
            for l in a:
                b.write(l)
            
        b.write(line)


with open('C:\\temp\section_data.txt','r') as s:
    data=s.readlines()
    fibs_s=int(data[0])
    fibs=int(data[1])
    nvol=int(data[2])
    dx=np.double(data[3])
    dy=np.double(data[4])
    dz=np.double(data[5])
    model_num=data[6]
    print(fibs,nvol)
#%%
mdb.ModelFromInputFile(name='Job-1', 
    inputFileName='C:/temp/Job-1.inp')


import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Job-1'].parts['PART-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['Job-1'].Material(name='Material-1')
mdb.models['Job-1'].materials['Material-1'].Density(table=((rho_f, ), ))
mdb.models['Job-1'].materials['Material-1'].Elastic(table=((E_f, Poisson_f), ))
mdb.models['Job-1'].materials['Material-1'].Conductivity(table=((Cond_f, ), ))
mdb.models['Job-1'].Material(name='Material-2')
mdb.models['Job-1'].materials['Material-2'].Density(table=((rho_m, ), ))
mdb.models['Job-1'].materials['Material-2'].Elastic(table=((E_m, Poisson_m), ))
mdb.models['Job-1'].materials['Material-2'].Conductivity(table=((Cond_m, ), ))
mdb.models['Job-1'].HomogeneousSolidSection(name='Section-1', 
    material='Material-1', thickness=None)
mdb.models['Job-1'].HomogeneousSolidSection(name='Section-2', 
    material='Material-2', thickness=None)

p = mdb.models['Job-1'].parts['PART-1']

if(fibs_s!=1):
    for i in range (1,fibs_s):
        region = p.sets['VOLUME'+str(i)]
        p = mdb.models['Job-1'].parts['PART-1']
        p.SectionAssignment(region=region, sectionName='Section-2')    
for i in range (fibs_s,fibs):

    region = p.sets['VOLUME'+str(i)]
    p = mdb.models['Job-1'].parts['PART-1']
    p.SectionAssignment(region=region, sectionName='Section-1')

for i in range (fibs,nvol+1):
    region = p.sets['VOLUME'+str(i)]
    p = mdb.models['Job-1'].parts['PART-1']
    p.SectionAssignment(region=region, sectionName='Section-2')

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
p1 = mdb.models['Mesh_model_2'].parts['PART-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
del mdb.models['Mesh_model_2']
session.viewports['Viewport: 1'].setValues(displayedObject=None)
del mdb.models['Model-1']
p = mdb.models['Job-1'].parts['PART-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Job-1'].rootAssembly
a.regenerate()
session.viewports['Viewport: 1'].setValues(displayedObject=a)
del mdb.jobs['Job-1']




# %% heat transfer
print('started theraml simulation in model num ',model_num)
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

elemType1 = mesh.ElemType(elemCode=DC3D4, elemLibrary=STANDARD)
# elemType1 = mesh.ElemType(elemCode=DC3D8, elemLibrary=STANDARD)
p = mdb.models['Job-1'].parts['PART-1']
z1 = p.elements
n_ele=len(z1)
elems1 = z1[0:int(n_ele)]
pickedRegions =(elems1, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
a = mdb.models['Job-1'].rootAssembly
a.regenerate()
# session.viewports['Viewport: 1'].setValues(displayedObject=a)
mdb.Job(name='Job-2', model='Job-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
    numDomains=2, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
    numCpus=2, numGPUs=0)

import sys
sys.path.insert(15, r'c:/temp/abaqus_plugins/MicroMechanics')
import microMechanics
from microMechanics.mmpBackend import Interface
from microMechanics.mmpBackend.mmpInterface.mmpRVEConstants import *
from microMechanics.mmpBackend.mmpKernel.mmpLibrary import mmpRVEGenerationUtils as RVEGenerationUtils
Interface.Loading.ThermalModelMaker(constraintType='PERIODIC', 
    modelName='Job-1', jobName='Job-2', doNotSubmit=False, 
    homogenizeProperties=(True, False, False))

mdb.jobs['Job-2'].waitForCompletion()

import sys
import itertools
from odbAccess import *
import csv

def thermal_calculator(steptoread,stepnum,dL):
    frametoread=steptoread.frames[-1]
    t=frametoread.fieldOutputs ['TEMP']
    v=frametoread.fieldOutputs ['IVOL']
    fl=frametoread.fieldOutputs ['HFL']
    Temp=[]
    elements=[]
    
    #xfl=fl.values[i].data[0] xflux value

    sum=0
    Vtot=0
    k=0
    for (i,j,m) in zip(fl.values,v.values,t.values) :
        elements.append(m.elementLabel)
        Temp.append(m.data)
        sum+=i.data[stepnum]*j.data
        Vtot+=j.data
        k+=1
    nele=max(elements)

    tdiff=max(Temp)-min(Temp)
    fbar=(sum/Vtot)
    K=fbar*(dL/tdiff)
    return nele,K

a = mdb.models['Job-1'].rootAssembly
# weight=a.getMassProperties()['mass']

odb =openOdb('job-2.odb')
K=[]
steptoread=odb.steps['HomogK_X-Initial']
dL=dx
nele,ki=thermal_calculator(steptoread,0,dL)
K.append(ki)
###############################################################
steptoread=odb.steps['HomogK_Y-Initial']
dL=dy
nele,ki=thermal_calculator(steptoread,1,dL)
K.append(ki)
###############################
steptoread=odb.steps['HomogK_Z-Initial']
dL=dz
nele,ki=thermal_calculator(steptoread,2,dL)
K.append(ki)
# print('V_f=',V_f)
print('Kx=',K)
# print('E=',E)
model_num=model_num.replace('\n','')
dir='C:/temp/model'+str(model_num)
with open(dir+'/sim_data.txt','a') as f:
    f.write('THERMAL'+'\n')
    f.write('n_ele='+str(n_ele)+'\n')
    f.write('Kx='+str(K[0])+'\n')
    f.write('Ky='+str(K[1])+'\n')
    f.write('Kz='+str(K[2])+'\n')





dir='C:\\temp\simulated_model.txt'
with open(dir,'w') as c:
    c.write(model_num)

# %%
