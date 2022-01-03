##      EasyPBC Ver. 1.3   (08/10/2018) updated on (03/11/2018) for error No. 4 and 2D check.
##      EasyPBC is an ABAQUS CAE plugin developed to estimate the homogenised effective elastic properties of user created periodic(RVE)
##      Copyright (C) 2018  Sadik Lafta Omairey
##
##      This program is distributed in the hope that it will be useful,
##      but WITHOUT ANY WARRANTY; without even the implied warranty of
##      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##      GNU General Public License for more details.
##
##      You should have received a copy of the GNU General Public License
##      along with this program.  If not, see <https://www.gnu.org/licenses/>.
##      Kindly do not re-distribute  
##      Citation: Omairey S, Dunning P, Sriramula S (2018) Development of an ABAQUS plugin tool for periodic RVE homogenisation.
##      Engineering with Computers. https://doi.org/10.1007/s00366-018-0616-4
##      Email sadik.omairey@gmail.com to obtain the latest version of the software.




## Importing ABAQUS Data and Python modules ##

from abaqus import *
from abaqusConstants import *
import __main__
import math
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
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import time
import os
import sys
import ctypes
import multiprocessing
path = os.getcwd()


## Plugin main GUI function ##

def feasypbc(part,inst,meshsens,E11,E22,E33,G12,G13,G23,CPU,onlyPBC):
        for T in (range(1)):
                start = time.time()
                modelName = part
                instanceName = inst
                upperName= inst.upper()
                
                fail = []
                keycheck2 =[inst]
                
                if part not in (mdb.models.keys()):
                        Er2=0
                        messageBox2 = ctypes.windll.user32.MessageBoxA
                        returnValue = messageBox2(Er2,'Model name is incorrect, please input the correct Model name.','EasyPBC Start-up error 02',0x30 | 0x0)
                        print('Start-up error 02. Refer EasyPBC user guide')
                        continue
                
                a = mdb.models[modelName].rootAssembly
                errorcheck1 = mdb.models[modelName].rootAssembly.instances.keys()
                if errorcheck1 == fail:
                        Er1=0
                        messageBox1 = ctypes.windll.user32.MessageBoxA
                        returnValue = messageBox1(Er1,'Model part is not created!\nPlease create part and try again','EasyPBC Start-up error 01',0x30 | 0x0)
                        print('Start-up error 01. Refer EasyPBC user guide')
                        continue
                
                if (mdb.models[modelName].rootAssembly.instances.keys()) != keycheck2:                       
                        Er3=0
                        messageBox3 = ctypes.windll.user32.MessageBoxA
                        returnValue = messageBox3(Er3,'Instance name is incorrect, please input the correct instance name.','EasyPBC Start-up error 03',0x30 | 0x0)
                        print('Start-up error 03. Refer EasyPBC user guide')
                        continue

                if CPU <= 0:
                        Er5=0
                        messageBox5 = ctypes.windll.user32.MessageBoxA
                        returnValue = messageBox5(Er5,'Specified number of CPUs is <= zero, please set it to a value larger than zero.','EasyPBC Start-up error 05',0x30 | 0x0)
                        print('Start-up error 05. Refer EasyPBC user guide')
                        continue

                
                CPUs = int(round(CPU))
                if CPUs > multiprocessing.cpu_count():
                        CPUs = multiprocessing.cpu_count()
                        print ('Warning: Specified number of CPUs is greater than the available. The maximum available number of CPUs is used (%s CPU(s)).' % CPUs)


                Nodeset = mdb.models[modelName].rootAssembly.instances[instanceName].nodes

                ## Start of sets creation ##                
                j = 0
                x=[]
                y=[]
                z=[]
                c1=[]
                c2=[]
                c3=[]
                c4=[]
                c5=[]
                c6=[]
                c7=[]
                c8=[]
                Max=[]
                ftedgexyz={}
                btedgexyz={}
                fbedgexyz={}
                bbedgexyz={}
                fledgexyz={}
                bledgexyz={}
                fredgexyz={}
                bredgexyz={}
                ltedgexyz={}
                rtedgexyz={}
                lbedgexyz={}
                rbedgexyz={}
                frontsxyz={}
                backsxyz={}
                topsxyz={}
                botsxyz={}
                leftsxyz={}
                rightsxyz={}
                frontbcxyz={}
                backbcxyz={}
                topbcxyz={}
                botbcxyz={}
                leftbcxyz={}
                rightbcxyz={}
                ftedge=[]
                btedge=[]
                fbedge=[]
                bbedge=[]
                fledge=[]
                fredge=[]
                bledge=[]
                bredge=[]
                ltedge=[]
                lbedge=[]
                rtedge=[]
                rbedge=[]
                fronts=[]
                backs=[]
                lefts=[]
                rights=[]
                tops=[]
                bots=[]
                backs=[]
                frontbc=[]
                backbc=[]
                leftbc=[]
                rightbc=[]
                topbc=[]
                botbc=[]
                backbc=[]
                errorset=[]
                coc1={}
                coc2={}
                coc3={}
                coc4={}
                coc5={}
                coc6={}
                coc7={}
                coc8={}

                error=False


                print ('----------------------------------')
                print ('-------- Start of EasyPBC --------')
                print ('----------------------------------')



                ## Identifying RVE size ##    
                for i in Nodeset:
                    x.insert(j,i.coordinates[0])
                    y.insert(j,i.coordinates[1])
                    z.insert(j,i.coordinates[2])
                    j=j+1



                errorcheck4 = x                   
                if not errorcheck4:
                        Er4=0
                        messageBox4 = ctypes.windll.user32.MessageBoxA
                        returnValue = messageBox4(Er4,'Instance not detected! Make sure:\n1- Instance is created;\n2- Double click on instance to refresh it before running EasyPBC;\n3- Part/instnace is meshed.','EasyPBC Start-up error 04',0x30 | 0x0)
                        print('Start-up error 04. Refer EasyPBC user guide')
                        continue


                Max = max(x)
                May = max(y)
                Maz = max(z)
                Mnx = min(x)
                Mny = min(y)
                Mnz = min(z)


                
                if (Maz - Mnz)<=meshsens:  ## 2D Model Check

                ############################################ 2D Section #####################################

                        L=abs(Max-Mnx)
                        H=abs(May-Mny)
                        
                        Dispx = L*0.2
                        Dispy = H*0.2

                        #### This part is commented out to reduce complexity for users whom merge parts into a single instance. ##############
                        
##                        partName = mdb.models[modelName].rootAssembly.instances[inst].partName                ## Extracts part name used in the instance
##                        AssSec = mdb.models[modelName].parts[partName].sectionAssignments[0].sectionName
##                        if AssSec == None:
##                                print 'Warning: No shell sections detected, please create a section and try again!!'
##                                print '         Refer to error 09 troubleshooting in easyPBC user guide.'
##                                error=True
##                                continue

                        SecName = mdb.models['Model-1'].sections.keys()[0]                                      ## Finding the name of first section            
                        Thikness = mdb.models[modelName].sections[SecName].thickness
                        skipAtt = False
                        
                        if Thikness == None or Thikness == 0:
                                Thikness = 1
                                print 'Attention: EasyPBC did not detected a shell thickness. Thus, it assumes thickness is equal to 1.0 unit length'
                                print '           If another thickness value is desired, divide the elastic property(ies) by the actual thickness value.'
                                skipAtt = True
                        if skipAtt == False:
                                print 'Attention: EasyPBC detected a shell thickness of %s unit length from the first section in ABAQUS property module.' % Thikness
                                print '           If this value is incorrect (this section is not used), multiply the elastic property(ies) %s and divide it by the actual value.' % Thikness

                        
                        ## Creating Ref. Points (Same R.Ps as in the 3D case) ##
                        for i in a.features.keys():
                            if i.startswith('RP'):
                                del a.features['%s' % (i)]
                        a.ReferencePoint(point=(Max+0.8*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP6: G23
                        a.ReferencePoint(point=(Max+0.6*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP5: G13
                        a.ReferencePoint(point=(Max+0.4*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP4: G12
                        a.ReferencePoint(point=(Max+0.2*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP3: Rigid body movement X-axis
                        a.ReferencePoint(point=(Max-0.5*(Max-Mnx), May-0.5*(May-Mny), Maz+0.2*abs(Maz-Mnz)))  ## RP2: Rigid body movement Z-axis
                        a.ReferencePoint(point=(Max-0.5*(Max-Mnx), May+0.2*abs(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP1: Rigid body movement Y-axis

                        r1 = a.referencePoints

                        ## Naming Ref. Points ##
                        d=1
                        for i in r1.keys():
                            refPoints1=(r1[i], )
                            a.Set(referencePoints=refPoints1, name='RP%s' % (d))
                            d=d+1
                          
                        ## Identifying boundary nodes ##
                        for i in Nodeset:
                            if (Mnx+meshsens) < i.coordinates[0] < (Max-meshsens) and (Mny+meshsens) < i.coordinates[1] < (May-meshsens):
                                continue

                            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens:
                                c1.insert(0,i.label)
                                coc1[i.label]=[i.coordinates[0], i.coordinates[1]]
                            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens:
                                c2.insert(0,i.label)
                                coc2[i.label]=[i.coordinates[0], i.coordinates[1]]
                            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens:
                                c5.insert(0,i.label)
                                coc5[i.label]=[i.coordinates[0], i.coordinates[1]]
                            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens:
                                c6.insert(0,i.label)
                                coc6[i.label]=[i.coordinates[0], i.coordinates[1]]
                            if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                                frontsxyz[i.label]=[i.coordinates[0], i.coordinates[1]]
                            if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                                backsxyz[i.label]=[i.coordinates[0], i.coordinates[1]] 
                            if abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                                topsxyz[i.label]=[i.coordinates[0], i.coordinates[1]]
                            if abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                                botsxyz[i.label]=[i.coordinates[0], i.coordinates[1]]

                        ## Checking number of nodes of opposite/associated sets ##
                        if len(frontsxyz) != len(backsxyz):
                         print 'Warning: Number of Nodes in Front surface (fronts) not equal to number of nodes in Back surface (backs). These sets will not be created!!'
                         print '         Refer to error 06 troubleshooting in easyPBC user guide.'
                         frontsxyz={}
                         error=True
                        if len(topsxyz) != len(botsxyz):
                         print 'Warning: Number of Nodes in Top surface (tops) not equal to number of nodes in Bottom surface (bots). These sets will not be created!!'
                         print '         Refer to error 06 in easyPBC user guide.'
                         topsxyz={}
                         error=True
       
                        ## Sorting and appending sets ##
                        for i in frontsxyz.keys():
                                for k in backsxyz.keys():
                                        if abs(frontsxyz[i][1] - backsxyz[k][1])<=meshsens:
                                                fronts.append(i)
                                                backs.append(k)

                        if len(frontsxyz)!= len(fronts) or len(backsxyz)!= len(backs):
                            print 'Warning: Node(s) in Front and/or Back surface (fronts and/or backs) was not imported. effected sets will not be created!!'
                            print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                            for i, k in zip(frontsxyz.keys(),backsxyz.keys()):
                                    if i not in fronts:
                                            errorset.append(i)
                                    if k not in backs:
                                            errorset.append(k)
                            fronts=[]
                            backs=[]
                            error=True                    
                        if len(fronts)!=len(set(fronts)) or len(backs)!=len(set(backs)):
                            print 'Warning: Node(s) in either Front or Back surface (fronts or backs) being linked with more than one opposite node. effected sets will not be created!!'
                            print '         Refer to error 08 in easyPBC user guide.'
                            fronts=[]
                            backs=[]
                            error=True

                        for i in topsxyz.keys():
                            for k in botsxyz.keys():
                                if abs(topsxyz[i][0] - botsxyz[k][0]) <=meshsens:
                                    tops.append(i)
                                    bots.append(k)
                        if len(topsxyz)!= len(tops) or len(botsxyz)!= len(bots):
                            print 'Warning: Node(s) in Top and/or Bottom surface (tops and/or bots) was not imported. effected sets will not be created!!'
                            print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                            for i, k in zip(topsxyz.keys(),botsxyz.keys()):
                                    if i not in tops:
                                            errorset.append(i)
                                    if k not in bots:
                                            errorset.append(k)
                            tops=[]
                            bots=[]
                            error=True
                        if len(tops)!=len(set(tops)) or len(bots)!=len(set(bots)):
                            print 'Warning: Node(s) in either Top or Bottom surface (tops or bots) being linked with more than one opposite node. effected sets will not be created!!'
                            print '         Refer to error 08 in easyPBC user guide.'
                            tops=[]
                            bots=[]
                            error=True



                        ## Creating ABAQUS sets ##
                        a.SetFromNodeLabels(name='c1', nodeLabels=((instanceName,c1),))
                        a.SetFromNodeLabels(name='c2', nodeLabels=((instanceName,c2),))
                        a.SetFromNodeLabels(name='c5', nodeLabels=((instanceName,c5),))
                        a.SetFromNodeLabels(name='c6', nodeLabels=((instanceName,c6),))
                        a.SetFromNodeLabels(name='fronts', nodeLabels=((instanceName,fronts),))
                        a.SetFromNodeLabels(name='backs', nodeLabels=((instanceName,backs),))
                        a.SetFromNodeLabels(name='tops', nodeLabels=((instanceName,tops),))
                        a.SetFromNodeLabels(name='bots', nodeLabels=((instanceName,bots),))
                        print ('------ End of Sets Creation ------')

                        ## Extracting model mass ##
                        prop=mdb.models[modelName].rootAssembly.getMassProperties()
                        mass=prop['mass']

                        a = mdb.models[modelName].rootAssembly
                        Nodeset = mdb.models[modelName].rootAssembly.instances[instanceName].nodes
                        mdb.models[modelName].StaticStep(name='Step-1', previous='Initial')


                        ## Creating single-node ABAQUS sets ##                
                        if error==False:
                                
                                if E11==True or E22==True or E33==True or G12==True or G13==True or G23==True:
                                        for i,k in zip(tops,bots):
                                            a.SetFromNodeLabels(name='tops%s' % (i), nodeLabels=((instanceName,[i]),))
                                            a.SetFromNodeLabels(name='bots%s' % (k), nodeLabels=((instanceName,[k]),))

                                        for i,k in zip(fronts,backs):
                                            a.SetFromNodeLabels(name='fronts%s' % (i), nodeLabels=((instanceName,[i]),))
                                            a.SetFromNodeLabels(name='backs%s' % (k), nodeLabels=((instanceName,[k]),))


                                ## Creating constraints for elastic moduli ##
                                if E11==True or E22==True:
                                        for i in mdb.models[modelName].constraints.keys():
                                                del mdb.models[modelName].constraints[i]
                                                
                                        for i,k in zip(tops,bots):
                                            mdb.models[modelName].Equation(name='E-1-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 1), (-1.0, 'bots%s'%k, 1)))
                                        for i,k in zip(tops,bots):
                                            mdb.models[modelName].Equation(name='E-2-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 2), (-1.0, 'bots%s'%k, 2),(-1.0, 'RP5', 2)))


                                        for i,k in zip(fronts,backs):
                                            mdb.models[modelName].Equation(name='E-1-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 1), (-1.0, 'backs%s'%k, 1),(-1.0, 'RP4', 1)))
                                        for i,k in zip(fronts,backs):
                                            mdb.models[modelName].Equation(name='E-2-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 2), (-1.0, 'backs%s'%k, 2)))


                                        mdb.models[modelName].Equation(name='E-1-c62', terms=((1.0, 'c6', 1), (-1.0, 'c2', 1)))
                                        mdb.models[modelName].Equation(name='E-1-c21', terms=((1.0, 'c2', 1), (-1.0, 'c1', 1),(1.0, 'RP4', 1)))
                                        mdb.models[modelName].Equation(name='E-1-c15', terms=((1.0, 'c1', 1), (-1.0, 'c5', 1)))
                                       
                                        mdb.models[modelName].Equation(name='E-2-c62', terms=((1.0, 'c6', 2), (-1.0, 'c2', 2),(1.0, 'RP5', 2)))
                                        mdb.models[modelName].Equation(name='E-2-c21', terms=((1.0, 'c2', 2), (-1.0, 'c1', 2)))
                                        mdb.models[modelName].Equation(name='E-2-c15', terms=((1.0, 'c1', 2), (-1.0, 'c5', 2),(-1.0, 'RP5', 2)))


                                ## Elastic modulus E11 ##
                                if E11==True and onlyPBC == False:
                                        for i in mdb.models[modelName].loads.keys():
                                                del mdb.models[modelName].loads[i]
                                        for i in mdb.models[modelName].boundaryConditions.keys():
                                                del mdb.models[modelName].boundaryConditions[i]


                                        region = a.sets['RP4']
                                        mdb.models[modelName].DisplacementBC(name='E11-1', createStepName='Step-1', 
                                            region=region, u1=Dispx, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                            localCsys=None)

                                        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                            createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                            sectionPoints=DEFAULT, rebar=EXCLUDE)

                                        import os, glob

                                        mdb.Job(name='job-E11', model= modelName, description='', type=ANALYSIS, 
                                            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                        mdb.jobs['job-E11'].submit(consistencyChecking=OFF)
                                        mdb.jobs['job-E11'].waitForCompletion()
                                        o3 = session.openOdb(name='%s' % (path+'\job-E11.odb'))
                                     
                                        odb = session.odbs['%s' % (path+'\job-E11.odb')]

                                        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                        odbName=session.viewports[session.currentViewportName].odbDisplay.name



                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]

                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                            NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP4', ))

                                        forceE11 = 0
                                        for i in session.xyDataObjects.keys():
                                            forceE11=forceE11+(session.xyDataObjects[i][0][1])

                                        stressE11 = abs(forceE11/(H*Thikness))


                                        E11 = stressE11/(Dispx/L)                               


                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]
                                        
                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                            NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
                                        
                                        C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
                                        
                                        C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
                                        Dis = abs(C1U1new - C2U1new)

                                        E11U1= abs(L - Dis)


                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]
                                            

                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                            NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
                                        
                                        C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
                                        
                                        C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
                                        Dis = abs(C1U2new - C5U2new)

                                        E11U2= abs(H - Dis)

                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]



                                        V12=(E11U2/H)/(E11U1/L)


                                ## Elastic modulus E22 ##                                
                                if E11==False or onlyPBC == True:
                                        E11='N/A'
                                        V12='N/A'

                                if E22==True and onlyPBC == False:

                                        for i in mdb.models[modelName].loads.keys():
                                                del mdb.models[modelName].loads[i]
                                        for i in mdb.models[modelName].boundaryConditions.keys():
                                                del mdb.models[modelName].boundaryConditions[i]

                                        region = a.sets['RP5']
                                        mdb.models[modelName].DisplacementBC(name='E22-1', createStepName='Step-1', 
                                            region=region, u1=UNSET, u2=Dispy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                            localCsys=None)


                                        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                            createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                            sectionPoints=DEFAULT, rebar=EXCLUDE)

                                        import os, glob

                                        mdb.Job(name='job-E22', model= modelName, description='', type=ANALYSIS, 
                                            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                        mdb.jobs['job-E22'].submit(consistencyChecking=OFF)
                                        mdb.jobs['job-E22'].waitForCompletion()
                                        o3 = session.openOdb(name='%s' % (path+'\job-E22.odb'))
                                     
                                        odb = session.odbs['%s' % (path+'\job-E22.odb')]

                                        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                        odbName=session.viewports[session.currentViewportName].odbDisplay.name


                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]

                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                            NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('RP5', ))

                                        forceE22 = 0
                                        for i in session.xyDataObjects.keys():
                                            forceE22=forceE22+(session.xyDataObjects[i][0][1])

                                        stressE22 = abs(forceE22/(L*Thikness))


                                        E22 = stressE22/(Dispy/H)                                




                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]


                                        
                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                            NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
                                        
                                        C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
                                        
                                        C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
                                        Dis = abs(C1U1new - C2U1new)

                                        E22U1= abs(L - Dis)


                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]
                                            

                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                            NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
                                        
                                        C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
                                        
                                        C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
                                        Dis = abs(C1U2new - C5U2new)

                                        E22U2= abs(H - Dis)

                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]




                                        V21=(E22U1/L)/(E22U2/H)


                                ## Elastic modulus E33 ##
                                if E22==False or onlyPBC == True:
                                        E22='N/A'
                                        V21='N/A'



                                ## Creating constraints for shear moduli ##
                                if G12==True:
                                        if onlyPBC == False:
                                                for i in mdb.models[modelName].constraints.keys():
                                                        del mdb.models[modelName].constraints[i]

                                        for i,k in zip(tops,bots):
                                            mdb.models[modelName].Equation(name='G-1-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 1), (-1.0, 'bots%s'%k, 1),(-1.0, 'RP4', 1)))
                                        for i,k in zip(tops,bots):
                                            mdb.models[modelName].Equation(name='G-2-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 2), (-1.0, 'bots%s'%k, 2),(-1.0, 'RP1', 2)))
                                       

                                        for i,k in zip(fronts,backs):
                                            mdb.models[modelName].Equation(name='G-1-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 1), (-1.0, 'backs%s'%k, 1),(-1.0, 'RP3', 1)))
                                        for i,k in zip(fronts,backs):
                                            mdb.models[modelName].Equation(name='G-2-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 2), (-1.0, 'backs%s'%k, 2),(-1.0, 'RP4', 2)))


                                        mdb.models[modelName].Equation(name='G-1-c62', terms=((1.0, 'c6', 1), (-1.0, 'c2', 1),(1.0, 'RP4', 1)))
                                        mdb.models[modelName].Equation(name='G-1-c21', terms=((1.0, 'c2', 1), (-1.0, 'c1', 1),(1.0, 'RP3', 1)))
                                        mdb.models[modelName].Equation(name='G-1-c15', terms=((1.0, 'c1', 1), (-1.0, 'c5', 1),(-1.0, 'RP4', 1)))

                                            
                                        mdb.models[modelName].Equation(name='G-2-c62', terms=((1.0, 'c6', 2), (-1.0, 'c2', 2),(1.0, 'RP1', 2)))
                                        mdb.models[modelName].Equation(name='G-2-c21', terms=((1.0, 'c2', 2), (-1.0, 'c1', 2),(1.0, 'RP4', 2)))
                                        mdb.models[modelName].Equation(name='G-2-c15', terms=((1.0, 'c1', 2), (-1.0, 'c5', 2),(-1.0, 'RP1', 2)))


                                ## Shear modulus G12 ##
                                if G12==True and onlyPBC == False:
                                        for i in mdb.models[modelName].loads.keys():
                                                del mdb.models[modelName].loads[i]
                                        for i in mdb.models[modelName].boundaryConditions.keys():
                                                del mdb.models[modelName].boundaryConditions[i]

                                        region = a.sets['RP4']
                                        mdb.models[modelName].DisplacementBC(name='G12-1', createStepName='Step-1', 
                                            region=region, u1=Dispx, u2=Dispy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                            localCsys=None)


                                        regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                        mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                            createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                            sectionPoints=DEFAULT, rebar=EXCLUDE)

                                        import os, glob

                                        mdb.Job(name='job-G12', model= modelName, description='', type=ANALYSIS, 
                                            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                            scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                        mdb.jobs['job-G12'].submit(consistencyChecking=OFF)

                                        mdb.jobs['job-G12'].waitForCompletion()


                                        o3 = session.openOdb(name='%s' % (path+'\job-G12.odb'))


                                        odb = session.odbs['%s' % (path+'\job-G12.odb')]

                                        session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                        odbName=session.viewports[session.currentViewportName].odbDisplay.name


                                        for i in session.xyDataObjects.keys():
                                            del session.xyDataObjects['%s' % (i)]

                                        session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                            NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP4', ))

                                        forceG12 = 0
                                        for i in session.xyDataObjects.keys():
                                            forceG12=forceG12+(session.xyDataObjects[i][0][1])

                                        stressG12 = abs(forceG12/(L*Thikness))

                                        G12 = stressG12/((Dispx/H)+(Dispy/L))

                                ## Shear modulus G13 ##                                
                                if G12==False or onlyPBC == True:
                                        G12='N/A'


                                E33 ='N/A'
                                G13 ='N/A'
                                G23 ='N/A'
                                
                                density = 0
                                if mass != None:
                                        density = mass/(L*W*H)
                                
                                print ('----------------------------------------------------')
                                print ('----------------------------------------------------')
                                print ('The homogenised elastic properties:')
                                print ('E11=%s Stress units' % (E11))
                                print ('V12=%s ratio' % (V12))
                                print ('E22=%s Stress units' % (E22))
                                print ('V21=%s ratio' % (V21))
                                print ('G12=%s Stress units' % (G12))
                                print ('----------------------------------------------------')
                                print ('Total mass=%s Mass units' % (mass))
                                print ('Homogenised density=%s Density units' % (density))
                                print ('----------------------------------------------------')
                                print ('Processing duration %s seconds' % (time.time()-start))
                                print ('----------------------------------------------------')
                                
                                filename = ('%s_elastic_properties.txt' % part)
                                print ('The homogenised elastic properties are saved in ABAQUS Work Directory under %s' % filename)
                                f = open(filename,'w')
                                f.write('{0:^10}{1:^20}{2:^20}\n'.format('Property','Value','Unit'))
                                f.write('{0:^10}{1:^20}{2:^20}\n'.format('E11',E11,'Stress units'))
                                f.write('{0:^10}{1:^20}{2:^20}\n'.format('V12',V12,'ratio'))
                                f.write('{0:^10}{1:^20}{2:^20}\n'.format('E22',E22,'Stress units'))
                                f.write('{0:^10}{1:^20}{2:^20}\n'.format('V21',V21,'ratio'))
                                f.write('{0:^10}{1:^20}{2:^20}\n'.format('G12',G12,'Stress units'))

                                f.write ('Total mass=%s Mass units \n' % (mass))
                                f.write ('Homogenised density=%s Density units \n' % (density))
                           
                                f.write ('processing duration %s Seconds' % (time.time()-start))

                                f.close()

                                print ('Citation: Omairey S, Dunning P, Sriramula S (2018) Development of an ABAQUS plugin tool for periodic RVE homogenisation.')
                                print ('Engineering with Computers. https://doi.org/10.1007/s00366-018-0616-4')


                                filename = ('%s_elastic_properties(easycopy).txt' % part)
                                f = open(filename,'w')
                                f.write('{0:^10}\n'.format(E11))
                                f.write('{0:^10}\n'.format(E22))
                                f.write('{0:^10}\n'.format(G12))
                                f.write('{0:^10}\n'.format(V12))
                                f.write('{0:^10}\n'.format(V21))
                                f.write ('{0:^10}\n' .format(mass))
                                f.write ('{0:^10}\n' .format(density))                   
                                f.write ('{0:^10}\n' .format((time.time()-start)))

                                f.close()

                                print ('----------------------------------------------------')
                                if onlyPBC == True:
                                        print ('EasyPBC created Period Boundary Conditions only. For further investigation, used relevant Reference Points to apply loads/displacements based on your needs. Details on the use of Reference Points are illustrated in Table 1 of the referred paper.')

                                        
                                
                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]
                                print ('---------------------------------------')
                                print ('--------- End of EasyPBC (2D) ---------')
                                print ('---------------------------------------')

                                if len(session.odbData.keys()) >= 1:                              
                                        odb.close(odb, write=TRUE)
                                        a = mdb.models[modelName].rootAssembly
                                        session.viewports['Viewport: 1'].setValues(displayedObject=a)

                        continue


                ## 3D model ##########################################################
                
                L=abs(Max-Mnx)
                H=abs(May-Mny)
                W=abs(Maz-Mnz)
                
                Dispx = L*0.02
                Dispy = H*0.02
                Dispz = W*0.02
                
                ## Creating Ref. Points ##
                for i in a.features.keys():
                    if i.startswith('RP'):
                        del a.features['%s' % (i)]
                a.ReferencePoint(point=(Max+0.8*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP6: G23
                a.ReferencePoint(point=(Max+0.6*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP5: G13
                a.ReferencePoint(point=(Max+0.4*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP4: G12
                a.ReferencePoint(point=(Max+0.2*abs(Max-Mnx), May-0.5*(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP3: Rigid body movement X-axis
                a.ReferencePoint(point=(Max-0.5*(Max-Mnx), May-0.5*(May-Mny), Maz+0.2*abs(Maz-Mnz)))  ## RP2: Rigid body movement Z-axis
                a.ReferencePoint(point=(Max-0.5*(Max-Mnx), May+0.2*abs(May-Mny), Maz-0.5*(Maz-Mnz)))  ## RP1: Rigid body movement Y-axis

                r1 = a.referencePoints

                ## Naming Ref. Points ##
                d=1
                for i in r1.keys():
                    refPoints1=(r1[i], )
                    a.Set(referencePoints=refPoints1, name='RP%s' % (d))
                    d=d+1
                  
                ## Identifying boundary nodes ##
                for i in Nodeset:
                    if (Mnx+meshsens) < i.coordinates[0] < (Max-meshsens) and (Mny+meshsens) < i.coordinates[1] < (May-meshsens) and (Mnz+meshsens) < i.coordinates[2] < (Maz-meshsens):
                        continue
                    if abs(i.coordinates[0]-Max)<=meshsens:
                        frontbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[0]-Mnx)<=meshsens:
                        backbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[2]-Maz)<=meshsens:
                        leftbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[2]-Mnz)<=meshsens:
                        rightbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[1]-May)<=meshsens:
                        topbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[1]-Mny)<=meshsens:
                        botbcxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                        c1.insert(0,i.label)
                        coc1[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                        c2.insert(0,i.label)
                        coc2[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                        c3.insert(0,i.label)
                        coc3[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                        c4.insert(0,i.label)
                        coc4[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                        c5.insert(0,i.label)
                        coc5[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens:
                        c6.insert(0,i.label)
                        coc6[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                        c7.insert(0,i.label)
                        coc7[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens:
                        c8.insert(0,i.label)
                        coc8[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        ftedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        fbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        btedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        bbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                        fledgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                        fredgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                        bledgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens:
                        bredgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                        ltedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                        lbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                        rtedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                        rbedgexyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[0]-Max)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        frontsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[0]-Mnx)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        backsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]] 
                    if abs(i.coordinates[2]-Maz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                        leftsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]	
                    if abs(i.coordinates[2]-Mnz)<=meshsens and abs(i.coordinates[1]-May)>meshsens and abs(i.coordinates[1]-Mny)>meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens:
                        rightsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]   
                    if abs(i.coordinates[1]-May)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        topsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]
                    if abs(i.coordinates[1]-Mny)<=meshsens and abs(i.coordinates[0]-Max)>meshsens and abs(i.coordinates[0]-Mnx)>meshsens and abs(i.coordinates[2]-Maz)>meshsens and abs(i.coordinates[2]-Mnz)>meshsens:
                        botsxyz[i.label]=[i.coordinates[0], i.coordinates[1], i.coordinates[2]]

                ## Checking number of nodes of opposite/associated sets ##
                if len(frontsxyz) != len(backsxyz):
                 print 'Warning: Number of Nodes in Front surface (fronts) not equal to number of nodes in Back surface (backs). These sets will not be created!!'
                 print '         Refer to error 06 troubleshooting in easyPBC user guide.'
                 frontsxyz={}
                 error=True
                if len(topsxyz) != len(botsxyz):
                 print 'Warning: Number of Nodes in Top surface (tops) not equal to number of nodes in Bottom surface (bots). These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 topsxyz={}
                 error=True
                if len(leftsxyz) != len(rightsxyz):
                 print 'Warning: Number of Nodes in Left surface (lefts) not equal to number of nodes in Right surface (rights). These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 leftsxyz={}
                 error=True
                if len(ftedgexyz) != len(btedgexyz) or len(btedgexyz) != len(bbedgexyz) or len(bbedgexyz) != len(ftedgexyz):
                 print 'Warning: Number of nodes in front-top ,back-top, back-bottom, front-bottom (ftedge, btedge, bbedge and fbedge) are not equal. These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 ftedgexyz={}
                 error=True
                if len(fledgexyz) != len(bledgexyz) or len(bledgexyz) != len(bredgexyz) or len(bredgexyz) != len(fredgexyz):
                 print 'Warning: Number of nodes in front-left, back-left, back-right, front-right edge (fledge, bledge, bredge and fredge) are not equal. These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 fledgexyz={}
                 error=True
                if len(ltedgexyz) != len(rtedgexyz) or len(rtedgexyz) != len(rbedgexyz) or len(rbedgexyz) != len(lbedgexyz):
                 print 'Warning: Number of nodes in left-top, right-top, right-bottom, front-bottom edge (ltedge, rtedge, rbedge and fbedge). are not equal. These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 ltedgexyz={}
                 error=True
                if len(frontbcxyz) != len(backbcxyz):
                 print 'Warning: Number of Nodes in Front BC surface (frontbc) not equal to number of nodes in Back BC surface (backbc). These sets will not be created!!'
                 print '         Refer to error 06 troubleshooting in easyPBC user guide.'
                 frontbcxyz={}
                 error=True
                if len(topbcxyz) != len(botbcxyz):
                 print 'Warning: Number of Nodes in Top BC surface (topbc) not equal to number of nodes in Bottom BC surface (botbc). These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 topbcxyz={}
                 error=True
                if len(leftbcxyz) != len(rightbcxyz):
                 print 'Warning: Number of Nodes in Left BC surface (leftbc) not equal to number of nodes in Right BC surface (rightbc). These sets will not be created!!'
                 print '         Refer to error 06 in easyPBC user guide.'
                 leftbcxyz={}
                 error=True

                ## Sorting and appending sets ##
                for i in frontsxyz.keys():
                        for k in backsxyz.keys():
                                if abs(frontsxyz[i][1] - backsxyz[k][1])<=meshsens and abs(frontsxyz[i][2] - backsxyz[k][2])<=meshsens:
                                        fronts.append(i)
                                        backs.append(k)

                if len(frontsxyz)!= len(fronts) or len(backsxyz)!= len(backs):
                    print 'Warning: Node(s) in Front and/or Back surface (fronts and/or backs) was not imported. effected sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    for i, k in zip(frontsxyz.keys(),backsxyz.keys()):
                            if i not in fronts:
                                    errorset.append(i)
                            if k not in backs:
                                    errorset.append(k)
                    fronts=[]
                    backs=[]
                    error=True                    
                if len(fronts)!=len(set(fronts)) or len(backs)!=len(set(backs)):
                    print 'Warning: Node(s) in either Front or Back surface (fronts or backs) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    fronts=[]
                    backs=[]
                    error=True

                for i in topsxyz.keys():
                    for k in botsxyz.keys():
                        if abs(topsxyz[i][0] - botsxyz[k][0]) <=meshsens and abs(topsxyz[i][2] - botsxyz[k][2]) <=meshsens:
                            tops.append(i)
                            bots.append(k)
                if len(topsxyz)!= len(tops) or len(botsxyz)!= len(bots):
                    print 'Warning: Node(s) in Top and/or Bottom surface (tops and/or bots) was not imported. effected sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    for i, k in zip(topsxyz.keys(),botsxyz.keys()):
                            if i not in tops:
                                    errorset.append(i)
                            if k not in bots:
                                    errorset.append(k)
                    tops=[]
                    bots=[]
                    error=True
                if len(tops)!=len(set(tops)) or len(bots)!=len(set(bots)):
                    print 'Warning: Node(s) in either Top or Bottom surface (tops or bots) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    tops=[]
                    bots=[]
                    error=True


                for i in leftsxyz.keys():
                    for k in rightsxyz.keys():
                        if abs(leftsxyz[i][0] - rightsxyz[k][0])<=meshsens and abs(leftsxyz[i][1] - rightsxyz[k][1]) <=meshsens:
                            lefts.append(i)
                            rights.append(k)
                if len(leftsxyz)!= len(lefts) or len(rightsxyz)!= len(rights):
                    print 'Warning: Node(s) in Left and/or Right surface (lefts and/or rights) was not imported. effected sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    for i, k in zip(leftsxyz.keys(),rightsxyz.keys()):
                            if i not in lefts:
                                    errorset.append(i)
                            if k not in rights:
                                    errorset.append(k)                    
                    lefts=[]
                    rights=[]
                    error=True                    
                if len(lefts)!=len(set(lefts)) or len(rights)!=len(set(rights)):
                    print 'Warning: Node(s) in either Left or Right surface (lefts or rights) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    lefts=[]
                    rights=[]
                    error=True

                for i in frontbcxyz.keys():
                    for k in backbcxyz.keys():
                        if abs(frontbcxyz[i][1] - backbcxyz[k][1])<=meshsens and abs(frontbcxyz[i][2] - backbcxyz[k][2])<=meshsens:
                            frontbc.append(i)
                            backbc.append(k)
                if len(frontbcxyz)!= len(frontbc) or len(backbcxyz)!= len(backbc):
                    print 'Warning: Node(s) in Front BC and/or Back BC surface (frontbc and/or backbc) was not imported. effected sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    for i, k in zip(frontbcxyz.keys(),backbcxyz.keys()):
                            if i not in frontbc:
                                    errorset.append(i)
                            if k not in backbc:
                                    errorset.append(k)
                    frontbc=[]
                    backbc=[]
                    error=True
                if len(frontbc)!=len(set(frontbc)) or len(backbc)!=len(set(backbc)):
                    print 'Warning: Node(s) in either Front BC or Back BC surface (frontbc or backbc) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    frontbc=[]
                    backbc=[]
                    error=True


                for i in topbcxyz.keys():
                    for k in botbcxyz.keys():
                        if abs(topbcxyz[i][0] - botbcxyz[k][0]) <=meshsens and abs(topbcxyz[i][2] - botbcxyz[k][2]) <=meshsens:
                            topbc.append(i)
                            botbc.append(k)
                if len(topbcxyz)!= len(topbc) or len(botbcxyz)!= len(botbc):
                    print 'Warning: Node(s) in Top BC and/or Bottom BC surface (topbc and/or botbc) was not imported. effected sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    for i, k in zip(topbcxyz.keys(),botbcxyz.keys()):
                            if i not in topbc:
                                    errorset.append(i)
                            if k not in botbc:
                                    errorset.append(k)
                    topbc=[]
                    botbc=[]
                    error=True
                if len(topbc)!=len(set(topbc)) or len(botbc)!=len(set(botbc)):
                    print 'Warning: Node(s) in either Top BC or Bottom BC surface (topbc or botbc) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    topbc=[]
                    botbc=[]
                    error=True


                for i in leftbcxyz.keys():
                    for k in rightbcxyz.keys():
                        if abs(leftbcxyz[i][0] - rightbcxyz[k][0])<=meshsens and abs(leftbcxyz[i][1] - rightbcxyz[k][1]) <=meshsens:
                            leftbc.append(i)
                            rightbc.append(k)
                if len(leftbcxyz)!= len(leftbc) or len(rightbcxyz)!= len(rightbc):
                    print 'Warning: Node(s) in Left BC and/or Right BC surface (lefts and/or rights) was not imported. effected sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    for i, k in zip(leftbcxyz.keys(),rightbcxyz.keys()):
                            if i not in leftbc:
                                    errorset.append(i)
                            if k not in rightbc:
                                    errorset.append(k)            
                    leftbc=[]
                    rightbc=[]
                    error=True
                if len(leftbc)!=len(set(leftbc)) or len(rightbc)!=len(set(rightbc)):
                    print 'Warning: Node(s) in either Left BC or Right BC surface (leftbc or rightbc) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    leftbc=[]
                    rightbc=[]
                    error=True


                for i in ftedgexyz.keys():
                    for k in btedgexyz.keys():
                        if abs(ftedgexyz[i][1] - btedgexyz[k][1])<=meshsens and abs(ftedgexyz[i][2] - btedgexyz[k][2])<=meshsens:
                            ftedge.append(i)
                            btedge.append(k)
                for i in btedge:
                    for k in bbedgexyz.keys():
                        if abs(btedgexyz[i][0] - bbedgexyz[k][0]) <=meshsens and abs(btedgexyz[i][2] - bbedgexyz[k][2]) <=meshsens:
                            bbedge.append(k)    
                for i in bbedge:
                    for k in fbedgexyz.keys():
                        if abs(bbedgexyz[i][1] - fbedgexyz[k][1]) <=meshsens and abs(bbedgexyz[i][2] - fbedgexyz[k][2]) <=meshsens:
                            fbedge.append(k) 
                if len(ftedge)!=len(set(ftedge)) or len(btedge)!=len(set(btedge)) or len(bbedge)!=len(set(bbedge)) or len(fbedge)!=len(set(fbedge)):
                    print 'Warning: Node(s) in either front-top, back-top, back-bottom and front-bottom edge(ftedge, btedge, bbedge and fbedge) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    ftedge=[]
                    btedge=[]
                    bbedg=[]
                    fbedge=[]
                    error==True
                if len(ftedgexyz)!= len(ftedge) or len(btedgexyz)!= len(btedge) or len(bbedgexyz)!= len(bbedge) or len(fbedgexyz)!= len(fbedge):
                    print 'Warning: Node(s) in front-top, back-top, back-bottom and front-bottom edge(ftedge, btedge, bbedge and fbedge) were not imported. these sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    ftedge=[]
                    btedge=[]
                    bbedg=[]
                    fbedge=[]
                    error=True

                for i in ltedgexyz.keys():
                    for k in rtedgexyz.keys():
                        if abs(ltedgexyz[i][0] - rtedgexyz[k][0])<=meshsens and abs(ltedgexyz[i][1] - rtedgexyz[k][1])<=meshsens:
                            ltedge.append(i)
                            rtedge.append(k)
                for i in rtedge:
                    for k in rbedgexyz.keys():
                        if abs(rtedgexyz[i][0] - rbedgexyz[k][0])<=meshsens and abs(rtedgexyz[i][2] - rbedgexyz[k][2])<=meshsens:
                            rbedge.append(k)    
                for i in rbedge:
                    for k in lbedgexyz.keys():
                        if abs(rbedgexyz[i][0] - lbedgexyz[k][0])<=meshsens and abs(rbedgexyz[i][1] - lbedgexyz[k][1])<=meshsens:
                            lbedge.append(k) 

                if len(ltedge)!=len(set(ltedge)) or len(rtedge)!=len(set(rtedge)) or len(rbedge)!=len(set(rbedge)) or len(lbedge)!=len(set(lbedge)):
                    print 'Warning: Node(s) in either front-top, back-bottom and front-bottom edge(ltedge, rtedge, rbedge and lbedge) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    ltedge=[]
                    rtedge=[]
                    rbedg=[]
                    lbedge=[]
                    error=True

                if len(ltedgexyz)!= len(ltedge) or len(rtedgexyz)!= len(rtedge) or len(rbedgexyz)!= len(rbedge) or len(lbedgexyz)!= len(lbedge):
                    print 'Warning: Node(s) in left-top, right-top, right-bottom, left-bottom edge (ltedge, rtedge, rbedge and lbedge) were not imported. these sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user guide and created Error set (if applicable).'
                    ltedge=[]
                    rtedge=[]
                    rbedg=[]
                    lbedge=[]
                    error=True

                for i in fledgexyz.keys():
                    for k in bledgexyz.keys():
                        if abs(fledgexyz[i][1] - bledgexyz[k][1])<=meshsens and abs(fledgexyz[i][2] - bledgexyz[k][2])<=meshsens:
                            fledge.append(i)
                            bledge.append(k)
                for i in bledge:
                    for k in bredgexyz.keys():
                        if abs(bledgexyz[i][0] - bredgexyz[k][0])<=meshsens and abs(bledgexyz[i][1] - bredgexyz[k][1])<=meshsens:
                            bredge.append(k)    
                for i in bredge:
                    for k in fredgexyz.keys():
                        if abs(bredgexyz[i][1] - fredgexyz[k][1])<=meshsens and abs(bredgexyz[i][2] - fredgexyz[k][2])<=meshsens:
                            fredge.append(k) 

                if len(fledge)!=len(set(fledge)) or len(bledge)!=len(set(bledge)) or len(bredge)!=len(set(bredge)) or len(fredge)!=len(set(fredge)):
                    print 'Warning: Node(s) in either front-left, back-left, back-right and front-right edge(fledge, bledge, bredge and fredge) being linked with more than one opposite node. effected sets will not be created!!'
                    print '         Refer to error 08 in easyPBC user guide.'
                    fledge=[]
                    bledge=[]
                    bredg=[]
                    fredge=[]
                    error=True
                if len(fledgexyz)!= len(fledge) or len(bledgexyz)!= len(bledge) or len(bredgexyz)!= len(bredge) or len(fredgexyz)!= len(fredge):
                    print 'Warning: Node(s) in front-left, back-left, back-right and front-right edge (fledge, bledge, bredge and fredge) were not imported. these sets will not be created!!'
                    print '         Refer to error 07 in easyPBC user user guide.'
                    fledge=[]
                    bledge=[]
                    bredg=[]
                    fredge=[]
                    error=True

                ## Creating ABAQUS sets ##
                a.SetFromNodeLabels(name='c1', nodeLabels=((instanceName,c1),))
                a.SetFromNodeLabels(name='c2', nodeLabels=((instanceName,c2),))
                a.SetFromNodeLabels(name='c3', nodeLabels=((instanceName,c3),))
                a.SetFromNodeLabels(name='c4', nodeLabels=((instanceName,c4),))
                a.SetFromNodeLabels(name='c5', nodeLabels=((instanceName,c5),))
                a.SetFromNodeLabels(name='c6', nodeLabels=((instanceName,c6),))
                a.SetFromNodeLabels(name='c7', nodeLabels=((instanceName,c7),))
                a.SetFromNodeLabels(name='c8', nodeLabels=((instanceName,c8),))
                a.SetFromNodeLabels(name='ftedge', nodeLabels=((instanceName,ftedge),))
                a.SetFromNodeLabels(name='fbedge', nodeLabels=((instanceName,fbedge),))
                a.SetFromNodeLabels(name='btedge', nodeLabels=((instanceName,btedge),))
                a.SetFromNodeLabels(name='bbedge', nodeLabels=((instanceName,bbedge),))
                a.SetFromNodeLabels(name='fledge', nodeLabels=((instanceName,fledge),))
                a.SetFromNodeLabels(name='fredge', nodeLabels=((instanceName,fredge),))
                a.SetFromNodeLabels(name='bledge', nodeLabels=((instanceName,bledge),))
                a.SetFromNodeLabels(name='bredge', nodeLabels=((instanceName,bredge),))
                a.SetFromNodeLabels(name='ltedge', nodeLabels=((instanceName,ltedge),))
                a.SetFromNodeLabels(name='lbedge', nodeLabels=((instanceName,lbedge),))
                a.SetFromNodeLabels(name='rtedge', nodeLabels=((instanceName,rtedge),))
                a.SetFromNodeLabels(name='rbedge', nodeLabels=((instanceName,rbedge),))
                a.SetFromNodeLabels(name='fronts', nodeLabels=((instanceName,fronts),))
                a.SetFromNodeLabels(name='backs', nodeLabels=((instanceName,backs),))
                a.SetFromNodeLabels(name='lefts', nodeLabels=((instanceName,lefts),))
                a.SetFromNodeLabels(name='rights', nodeLabels=((instanceName,rights),))
                a.SetFromNodeLabels(name='tops', nodeLabels=((instanceName,tops),))
                a.SetFromNodeLabels(name='bots', nodeLabels=((instanceName,bots),))
                a.SetFromNodeLabels(name='frontbc', nodeLabels=((instanceName,frontbc),))
                a.SetFromNodeLabels(name='backbc', nodeLabels=((instanceName,backbc),))
                a.SetFromNodeLabels(name='leftbc', nodeLabels=((instanceName,leftbc),))
                a.SetFromNodeLabels(name='rightbc', nodeLabels=((instanceName,rightbc),))
                a.SetFromNodeLabels(name='topbc', nodeLabels=((instanceName,topbc),))
                a.SetFromNodeLabels(name='botbc', nodeLabels=((instanceName,botbc),))
                print ('------ End of Sets Creation ------')

                ## Extracting model mass ##
                prop=mdb.models[modelName].rootAssembly.getMassProperties()
                mass=prop['mass']

                a = mdb.models[modelName].rootAssembly
                Nodeset = mdb.models[modelName].rootAssembly.instances[instanceName].nodes
                mdb.models[modelName].StaticStep(name='Step-1', previous='Initial')


                ## Creating single-node ABAQUS sets ##                
                if error==False:
                        
                        if E11==True or E22==True or E33==True or G12==True or G13==True or G23==True:
                                for i,k in zip(tops,bots):
                                    a.SetFromNodeLabels(name='tops%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='bots%s' % (k), nodeLabels=((instanceName,[k]),))

                                for i,k in zip(fronts,backs):
                                    a.SetFromNodeLabels(name='fronts%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='backs%s' % (k), nodeLabels=((instanceName,[k]),))

                                for i,k in zip(lefts,rights):
                                    a.SetFromNodeLabels(name='lefts%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='rights%s' % (k), nodeLabels=((instanceName,[k]),))

                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    a.SetFromNodeLabels(name='ftedge%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='btedge%s' % (k), nodeLabels=((instanceName,[k]),))
                                    a.SetFromNodeLabels(name='bbedge%s' % (j), nodeLabels=((instanceName,[j]),))
                                    a.SetFromNodeLabels(name='fbedge%s' % (l), nodeLabels=((instanceName,[l]),))

                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    a.SetFromNodeLabels(name='fledge%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='bledge%s' % (k), nodeLabels=((instanceName,[k]),))
                                    a.SetFromNodeLabels(name='bredge%s' % (j), nodeLabels=((instanceName,[j]),))
                                    a.SetFromNodeLabels(name='fredge%s' % (l), nodeLabels=((instanceName,[l]),))

                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    a.SetFromNodeLabels(name='ltedge%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='lbedge%s' % (k), nodeLabels=((instanceName,[k]),))
                                    a.SetFromNodeLabels(name='rbedge%s' % (j), nodeLabels=((instanceName,[j]),))
                                    a.SetFromNodeLabels(name='rtedge%s' % (l), nodeLabels=((instanceName,[l]),))

                                for i,k in zip(topbc,botbc):
                                    a.SetFromNodeLabels(name='topbc%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='botbc%s' % (k), nodeLabels=((instanceName,[k]),))

                                for i,k in zip(frontbc,backbc):
                                    a.SetFromNodeLabels(name='frontbc%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='backbc%s' % (k), nodeLabels=((instanceName,[k]),))

                                for i,k in zip(leftbc,rightbc):
                                    a.SetFromNodeLabels(name='leftbc%s' % (i), nodeLabels=((instanceName,[i]),))
                                    a.SetFromNodeLabels(name='rightbc%s' % (k), nodeLabels=((instanceName,[k]),))

                        ## Creating constraints for elastic moduli ##
                        if E11==True or E22==True or E33==True:
                                for i in mdb.models[modelName].constraints.keys():
                                        del mdb.models[modelName].constraints[i]
                                        
                                for i,k in zip(tops,bots):
                                    mdb.models[modelName].Equation(name='E-1-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 1), (-1.0, 'bots%s'%k, 1)))
                                for i,k in zip(tops,bots):
                                    mdb.models[modelName].Equation(name='E-2-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 2), (-1.0, 'bots%s'%k, 2),(-1.0, 'RP5', 2)))
                                for i,k in zip(tops,bots):
                                    mdb.models[modelName].Equation(name='E-3-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 3), (-1.0, 'bots%s'%k, 3)))

                                for i,k in zip(lefts,rights):
                                    mdb.models[modelName].Equation(name='E-1-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 1), (-1.0, 'rights%s'%k, 1)))
                                for i,k in zip(lefts,rights):
                                    mdb.models[modelName].Equation(name='E-2-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 2), (-1.0, 'rights%s'%k, 2)))
                                for i,k in zip(lefts,rights):
                                    mdb.models[modelName].Equation(name='E-3-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 3), (-1.0, 'rights%s'%k, 3),(-1.0, 'RP6', 3)))

                                for i,k in zip(fronts,backs):
                                    mdb.models[modelName].Equation(name='E-1-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 1), (-1.0, 'backs%s'%k, 1),(-1.0, 'RP4', 1)))
                                for i,k in zip(fronts,backs):
                                    mdb.models[modelName].Equation(name='E-2-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 2), (-1.0, 'backs%s'%k, 2)))
                                for i,k in zip(fronts,backs):
                                    mdb.models[modelName].Equation(name='E-3-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 3), (-1.0, 'backs%s'%k, 3)))


                                mdb.models[modelName].Equation(name='E-1-c12', terms=((1.0, 'c6', 1), (-1.0, 'c2', 1)))
                                mdb.models[modelName].Equation(name='E-1-c23', terms=((1.0, 'c2', 1), (-1.0, 'c3', 1)))
                                mdb.models[modelName].Equation(name='E-1-c34', terms=((1.0, 'c3', 1), (-1.0, 'c4', 1),(1.0, 'RP4', 1)))
                                mdb.models[modelName].Equation(name='E-1-c45', terms=((1.0, 'c4', 1), (-1.0, 'c8', 1)))
                                mdb.models[modelName].Equation(name='E-1-c56', terms=((1.0, 'c8', 1), (-1.0, 'c5', 1)))
                                mdb.models[modelName].Equation(name='E-1-c67', terms=((1.0, 'c5', 1), (-1.0, 'c1', 1)))
                                mdb.models[modelName].Equation(name='E-1-c78', terms=((1.0, 'c1', 1), (-1.0, 'c7', 1),(-1.0, 'RP4', 1)))
                               
                                mdb.models[modelName].Equation(name='E-2-c12', terms=((1.0, 'c6', 2), (-1.0, 'c2', 2),(1.0, 'RP5', 2)))
                                mdb.models[modelName].Equation(name='E-2-c23', terms=((1.0, 'c2', 2), (-1.0, 'c3', 2)))
                                mdb.models[modelName].Equation(name='E-2-c34', terms=((1.0, 'c3', 2), (-1.0, 'c4', 2)))
                                mdb.models[modelName].Equation(name='E-2-c45', terms=((1.0, 'c4', 2), (-1.0, 'c8', 2),(-1.0, 'RP5', 2)))
                                mdb.models[modelName].Equation(name='E-2-c56', terms=((1.0, 'c8', 2), (-1.0, 'c5', 2)))
                                mdb.models[modelName].Equation(name='E-2-c67', terms=((1.0, 'c5', 2), (-1.0, 'c1', 2),(1.0, 'RP5', 2)))
                                mdb.models[modelName].Equation(name='E-2-c78', terms=((1.0, 'c1', 2), (-1.0, 'c7', 2),(-1.0, 'RP5', 2)))

                                mdb.models[modelName].Equation(name='E-3-c12', terms=((1.0, 'c6', 3), (-1.0, 'c2', 3)))
                                mdb.models[modelName].Equation(name='E-3-c23', terms=((1.0, 'c2', 3), (-1.0, 'c3', 3),(-1.0, 'RP6', 3)))
                                mdb.models[modelName].Equation(name='E-3-c34', terms=((1.0, 'c3', 3), (-1.0, 'c4', 3)))
                                mdb.models[modelName].Equation(name='E-3-c45', terms=((1.0, 'c4', 3), (-1.0, 'c8', 3)))
                                mdb.models[modelName].Equation(name='E-3-c56', terms=((1.0, 'c8', 3), (-1.0, 'c5', 3),(1.0, 'RP6', 3)))
                                mdb.models[modelName].Equation(name='E-3-c67', terms=((1.0, 'c5', 3), (-1.0, 'c1', 3)))
                                mdb.models[modelName].Equation(name='E-3-c78', terms=((1.0, 'c1', 3), (-1.0, 'c7', 3),(-1.0, 'RP6', 3)))
                                       

                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    mdb.models[modelName].Equation(name='E-1-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 1), (-1.0, 'btedge%s'%k, 1),(-1.0, 'RP4', 1)))
                                    mdb.models[modelName].Equation(name='E-1-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 1), (-1.0, 'bbedge%s'%j, 1)))
                                    mdb.models[modelName].Equation(name='E-1-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 1), (-1.0, 'fbedge%s'%l, 1),(1.0, 'RP4', 1)))
                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    mdb.models[modelName].Equation(name='E-2-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 2), (-1.0, 'btedge%s'%k, 2)))
                                    mdb.models[modelName].Equation(name='E-2-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 2), (-1.0, 'bbedge%s'%j, 2),(-1.0, 'RP5', 2)))
                                    mdb.models[modelName].Equation(name='E-2-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 2), (-1.0, 'fbedge%s'%l, 2)))
                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    mdb.models[modelName].Equation(name='E-3-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 3), (-1.0, 'btedge%s'%k, 3)))
                                    mdb.models[modelName].Equation(name='E-3-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 3), (-1.0, 'bbedge%s'%j, 3)))
                                    mdb.models[modelName].Equation(name='E-3-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 3), (-1.0, 'fbedge%s'%l, 3)))

                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    mdb.models[modelName].Equation(name='E-1-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 1), (-1.0, 'bledge%s'%k, 1),(-1.0, 'RP4', 1)))
                                    mdb.models[modelName].Equation(name='E-1-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 1), (-1.0, 'bredge%s'%j, 1)))
                                    mdb.models[modelName].Equation(name='E-1-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 1), (-1.0, 'fredge%s'%l, 1),(1.0, 'RP4', 1)))
                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    mdb.models[modelName].Equation(name='E-2-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 2), (-1.0, 'bledge%s'%k, 2)))
                                    mdb.models[modelName].Equation(name='E-2-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 2), (-1.0, 'bredge%s'%j, 2)))
                                    mdb.models[modelName].Equation(name='E-2-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 2), (-1.0, 'fredge%s'%l, 2)))
                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    mdb.models[modelName].Equation(name='E-3-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 3), (-1.0, 'bledge%s'%k, 3)))
                                    mdb.models[modelName].Equation(name='E-3-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 3), (-1.0, 'bredge%s'%j, 3),(-1.0, 'RP6', 3)))
                                    mdb.models[modelName].Equation(name='E-3-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 3), (-1.0, 'fredge%s'%l, 3)))

                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    mdb.models[modelName].Equation(name='E-1-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 1), (-1.0, 'lbedge%s'%k, 1)))
                                    mdb.models[modelName].Equation(name='E-1-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 1), (-1.0, 'rbedge%s'%j, 1)))
                                    mdb.models[modelName].Equation(name='E-1-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 1), (-1.0, 'rtedge%s'%l, 1)))                                    
                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    mdb.models[modelName].Equation(name='E-2-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 2), (-1.0, 'lbedge%s'%k, 2),(-1.0, 'RP5', 2)))
                                    mdb.models[modelName].Equation(name='E-2-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 2), (-1.0, 'rbedge%s'%j, 2)))
                                    mdb.models[modelName].Equation(name='E-2-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 2), (-1.0, 'rtedge%s'%l, 2),(1.0, 'RP5', 2)))
                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    mdb.models[modelName].Equation(name='E-3-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 3), (-1.0, 'lbedge%s'%k, 3)))
                                    mdb.models[modelName].Equation(name='E-3-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 3), (-1.0, 'rbedge%s'%j, 3),(-1.0, 'RP6', 3)))
                                    mdb.models[modelName].Equation(name='E-3-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 3), (-1.0, 'rtedge%s'%l, 3)))

                        ## Elastic modulus E11 ##
                        if E11==True and onlyPBC == False:
                                for i in mdb.models[modelName].loads.keys():
                                        del mdb.models[modelName].loads[i]
                                for i in mdb.models[modelName].boundaryConditions.keys():
                                        del mdb.models[modelName].boundaryConditions[i]


                                region = a.sets['RP4']
                                mdb.models[modelName].DisplacementBC(name='E11-1', createStepName='Step-1', 
                                    region=region, u1=Dispx, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                    createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                    sectionPoints=DEFAULT, rebar=EXCLUDE)

                                import os, glob

                                mdb.Job(name='job-E11', model= modelName, description='', type=ANALYSIS, 
                                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                    scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                mdb.jobs['job-E11'].submit(consistencyChecking=OFF)
                                mdb.jobs['job-E11'].waitForCompletion()
                                o3 = session.openOdb(name='%s' % (path+'\job-E11.odb'))
                             
                                odb = session.odbs['%s' % (path+'\job-E11.odb')]

                                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                odbName=session.viewports[session.currentViewportName].odbDisplay.name



                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                    NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP4', ))

                                forceE11 = 0
                                for i in session.xyDataObjects.keys():
                                    forceE11=forceE11+(session.xyDataObjects[i][0][1])

                                stressE11 = abs(forceE11/(H*W))


                                E11 = stressE11/(Dispx/L)                                


                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]
                                
                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
                                
                                C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
                                
                                C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
                                Dis = abs(C1U1new - C2U1new)

                                E11U1= abs(L - Dis)


                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]
                                    

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
                                
                                C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
                                
                                C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
                                Dis = abs(C1U2new - C5U2new)

                                E11U2= abs(H - Dis)

                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]


                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('C1','C4', ))
                                
                                C1U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][2]
                                
                                C4U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c4[0])][0][1] + coc4[(c4[0])][2]
                                Dis = abs(C1U3new - C4U3new)

                                E11U3= abs(W - Dis)


                                V12=(E11U2/H)/(E11U1/L)
                                V13=(E11U3/W)/(E11U1/L)



                        ## Elastic modulus E22 ##                                
                        if E11==False or onlyPBC == True:
                                E11='N/A'
                                V12='N/A'
                                V13='N/A'

                        if E22==True and onlyPBC == False:

                                for i in mdb.models[modelName].loads.keys():
                                        del mdb.models[modelName].loads[i]
                                for i in mdb.models[modelName].boundaryConditions.keys():
                                        del mdb.models[modelName].boundaryConditions[i]

                                region = a.sets['RP5']
                                mdb.models[modelName].DisplacementBC(name='E22-1', createStepName='Step-1', 
                                    region=region, u1=UNSET, u2=Dispy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)


                                regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                    createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                    sectionPoints=DEFAULT, rebar=EXCLUDE)

                                import os, glob

                                mdb.Job(name='job-E22', model= modelName, description='', type=ANALYSIS, 
                                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                    scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                mdb.jobs['job-E22'].submit(consistencyChecking=OFF)
                                mdb.jobs['job-E22'].waitForCompletion()
                                o3 = session.openOdb(name='%s' % (path+'\job-E22.odb'))
                             
                                odb = session.odbs['%s' % (path+'\job-E22.odb')]

                                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                odbName=session.viewports[session.currentViewportName].odbDisplay.name


                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                    NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('RP5', ))

                                forceE22 = 0
                                for i in session.xyDataObjects.keys():
                                    forceE22=forceE22+(session.xyDataObjects[i][0][1])

                                stressE22 = abs(forceE22/(W*L))


                                E22 = stressE22/(Dispy/H)                                




                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]


                                
                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
                                
                                C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
                                
                                C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
                                Dis = abs(C1U1new - C2U1new)

                                E22U1= abs(L - Dis)


                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]
                                    

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
                                
                                C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
                                
                                C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
                                Dis = abs(C1U2new - C5U2new)

                                E22U2= abs(H - Dis)

                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]


                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('C1','C4', ))
                                
                                C1U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][2]
                                
                                C4U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c4[0])][0][1] + coc4[(c4[0])][2]
                                Dis = abs(C1U3new - C4U3new)

                                E22U3= abs(W - Dis)


                                V21=(E22U1/L)/(E22U2/H)
                                V23=(E22U3/W)/(E22U2/H)


                        ## Elastic modulus E33 ##
                        if E22==False or onlyPBC == True:
                                E22='N/A'
                                V21='N/A'
                                V23='N/A'


                        if E33==True and onlyPBC == False:
                                for i in mdb.models[modelName].loads.keys():
                                        del mdb.models[modelName].loads[i]
                                for i in mdb.models[modelName].boundaryConditions.keys():
                                        del mdb.models[modelName].boundaryConditions[i]


                                region = a.sets['RP6']
                                mdb.models[modelName].DisplacementBC(name='E33-1', createStepName='Step-1', 
                                    region=region, u1=UNSET, u2=UNSET, u3=Dispz, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                    createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                    sectionPoints=DEFAULT, rebar=EXCLUDE)

                                import os, glob

                                mdb.Job(name='job-E33', model= modelName, description='', type=ANALYSIS, 
                                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                    scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                mdb.jobs['job-E33'].submit(consistencyChecking=OFF)
                                mdb.jobs['job-E33'].waitForCompletion()
                                o3 = session.openOdb(name='%s' % (path+'\job-E33.odb'))
                             
                                odb = session.odbs['%s' % (path+'\job-E33.odb')]

                                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                odbName=session.viewports[session.currentViewportName].odbDisplay.name

                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                    NODAL, ((COMPONENT, 'RF3'), )), ), nodeSets=('RP6', ))

                                forceE33 = 0
                                for i in session.xyDataObjects.keys():
                                    forceE33=forceE33+(session.xyDataObjects[i][0][1])

                                stressE33 = abs(forceE33/(H*L))


                                E33 = stressE33/(Dispz/W)                                

                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]


                                
                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=('C1','C2', ))
                                
                                C1U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][0]
                                
                                C2U1new = session.xyDataObjects['U:U1 PI: %s N: %s' % (upperName,c2[0])][0][1] + coc2[(c2[0])][0]
                                Dis = abs(C1U1new - C2U1new)

                                E33U1= abs(L - Dis)


                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]
                                    

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U2'), )), ), nodeSets=('C1','C5', ))
                                
                                C1U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][1]
                                
                                C5U2new = session.xyDataObjects['U:U2 PI: %s N: %s' % (upperName,c5[0])][0][1] + coc5[(c5[0])][1]
                                Dis = abs(C1U2new - C5U2new)

                                E33U2= abs(H - Dis)

                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]


                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
                                    NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('C1','C4', ))
                                
                                C1U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c1[0])][0][1] + coc1[(c1[0])][2]
                                
                                C4U3new = session.xyDataObjects['U:U3 PI: %s N: %s' % (upperName,c4[0])][0][1] + coc4[(c4[0])][2]
                                Dis = abs(C1U3new - C4U3new)

                                E33U3= abs(W - Dis)


                                V31=(E33U1/L)/(E33U3/W)
                                V32=(E33U2/H)/(E33U3/W)


                        if E33==False or onlyPBC == True:
                                E33='N/A'
                                V31='N/A'
                                V32='N/A'

                        ## Creating constraints for shear moduli ##
                        if G12==True or G13==True or G23==True:
                                if onlyPBC == False:
                                        for i in mdb.models[modelName].constraints.keys():
                                                del mdb.models[modelName].constraints[i]

                                for i,k in zip(tops,bots):
                                    mdb.models[modelName].Equation(name='G-1-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 1), (-1.0, 'bots%s'%k, 1),(-1.0, 'RP4', 1)))
                                for i,k in zip(tops,bots):
                                    mdb.models[modelName].Equation(name='G-2-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 2), (-1.0, 'bots%s'%k, 2),(-1.0, 'RP1', 2)))
                                for i,k in zip(tops,bots):
                                    mdb.models[modelName].Equation(name='G-3-tops-bots%s'%i, terms=((1.0, 'tops%s'%i, 3), (-1.0, 'bots%s'%k, 3),(-1.0, 'RP6', 3)))

                               
                                for i,k in zip(lefts,rights):
                                    mdb.models[modelName].Equation(name='G-1-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 1), (-1.0, 'rights%s'%k, 1),(-1.0, 'RP5', 1)))
                                for i,k in zip(lefts,rights):
                                    mdb.models[modelName].Equation(name='G-2-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 2), (-1.0, 'rights%s'%k, 2),(-1.0, 'RP6', 2)))
                                for i,k in zip(lefts,rights):
                                    mdb.models[modelName].Equation(name='G-3-lefts-rights%s'%i, terms=((1.0, 'lefts%s'%i, 3), (-1.0, 'rights%s'%k, 3),(-1.0, 'RP2', 3)))

                                for i,k in zip(fronts,backs):
                                    mdb.models[modelName].Equation(name='G-1-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 1), (-1.0, 'backs%s'%k, 1),(-1.0, 'RP3', 1)))
                                for i,k in zip(fronts,backs):
                                    mdb.models[modelName].Equation(name='G-2-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 2), (-1.0, 'backs%s'%k, 2),(-1.0, 'RP4', 2)))
                                for i,k in zip(fronts,backs):
                                    mdb.models[modelName].Equation(name='G-3-fronts-backs%s'%i, terms=((1.0, 'fronts%s'%i, 3), (-1.0, 'backs%s'%k, 3),(-1.0, 'RP5', 3)))


                                mdb.models[modelName].Equation(name='G-1-c12', terms=((1.0, 'c6', 1), (-1.0, 'c2', 1),(1.0, 'RP4', 1)))
                                mdb.models[modelName].Equation(name='G-1-c23', terms=((1.0, 'c2', 1), (-1.0, 'c3', 1),(-1.0, 'RP5', 1)))
                                mdb.models[modelName].Equation(name='G-1-c34', terms=((1.0, 'c3', 1), (-1.0, 'c4', 1),(1.0, 'RP3', 1)))
                                mdb.models[modelName].Equation(name='G-1-c45', terms=((1.0, 'c4', 1), (-1.0, 'c8', 1),(-1.0, 'RP4', 1)))
                                mdb.models[modelName].Equation(name='G-1-c56', terms=((1.0, 'c8', 1), (-1.0, 'c5', 1),(1.0, 'RP5', 1)))
                                mdb.models[modelName].Equation(name='G-1-c67', terms=((1.0, 'c5', 1), (-1.0, 'c1', 1),(1.0, 'RP4', 1)))
                                mdb.models[modelName].Equation(name='G-1-c78', terms=((1.0, 'c1', 1), (-1.0, 'c7', 1),(-1.0, 'RP3', 1),(-1.0, 'RP4', 1),(-1.0, 'RP5', 1)))

                                    
                                mdb.models[modelName].Equation(name='G-2-c12', terms=((1.0, 'c6', 2), (-1.0, 'c2', 2),(1.0, 'RP1', 2)))
                                mdb.models[modelName].Equation(name='G-2-c23', terms=((1.0, 'c2', 2), (-1.0, 'c3', 2),(-1.0, 'RP6', 2)))
                                mdb.models[modelName].Equation(name='G-2-c34', terms=((1.0, 'c3', 2), (-1.0, 'c4', 2),(1.0, 'RP4', 2)))
                                mdb.models[modelName].Equation(name='G-2-c45', terms=((1.0, 'c4', 2), (-1.0, 'c8', 2),(-1.0, 'RP1', 2)))
                                mdb.models[modelName].Equation(name='G-2-c56', terms=((1.0, 'c8', 2), (-1.0, 'c5', 2),(1.0, 'RP6', 2)))
                                mdb.models[modelName].Equation(name='G-2-c67', terms=((1.0, 'c5', 2), (-1.0, 'c1', 2),(1.0, 'RP1', 2)))
                                mdb.models[modelName].Equation(name='G-2-c78', terms=((1.0, 'c1', 2), (-1.0, 'c7', 2),(-1.0, 'RP1', 2),(-1.0, 'RP4', 2),(-1.0, 'RP6', 2)))

        
                                mdb.models[modelName].Equation(name='G-3-c12', terms=((1.0, 'c6', 3), (-1.0, 'c2', 3),(1.0, 'RP6', 3)))
                                mdb.models[modelName].Equation(name='G-3-c23', terms=((1.0, 'c2', 3), (-1.0, 'c3', 3),(-1.0, 'RP2', 3)))
                                mdb.models[modelName].Equation(name='G-3-c34', terms=((1.0, 'c3', 3), (-1.0, 'c4', 3),(1.0, 'RP5', 3)))
                                mdb.models[modelName].Equation(name='G-3-c45', terms=((1.0, 'c4', 3), (-1.0, 'c8', 3),(-1.0, 'RP6', 3)))
                                mdb.models[modelName].Equation(name='G-3-c56', terms=((1.0, 'c8', 3), (-1.0, 'c5', 3),(1.0, 'RP2', 3)))
                                mdb.models[modelName].Equation(name='G-3-c67', terms=((1.0, 'c5', 3), (-1.0, 'c1', 3),(1.0, 'RP6', 3)))
                                mdb.models[modelName].Equation(name='G-3-c78', terms=((1.0, 'c1', 3), (-1.0, 'c7', 3),(-1.0, 'RP2', 3),(-1.0, 'RP5', 3),(-1.0, 'RP6', 3)))
                 

                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    mdb.models[modelName].Equation(name='G-1-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 1), (-1.0, 'btedge%s'%k, 1),(-1.0, 'RP3', 1)))
                                    mdb.models[modelName].Equation(name='G-1-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 1), (-1.0, 'bbedge%s'%j, 1),(-1.0, 'RP4', 1)))
                                    mdb.models[modelName].Equation(name='G-1-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 1), (-1.0, 'fbedge%s'%l, 1),(1.0, 'RP3', 1)))
                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    mdb.models[modelName].Equation(name='G-2-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 2), (-1.0, 'btedge%s'%k, 2),(-1.0, 'RP4', 2)))
                                    mdb.models[modelName].Equation(name='G-2-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 2), (-1.0, 'bbedge%s'%j, 2),(-1.0, 'RP1', 2)))
                                    mdb.models[modelName].Equation(name='G-2-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 2), (-1.0, 'fbedge%s'%l, 2),(1.0, 'RP4', 2)))
                                for i,k,j,l in zip(ftedge,btedge,bbedge,fbedge):
                                    mdb.models[modelName].Equation(name='G-3-ftedge-btedge%s'%i, terms=((1.0, 'ftedge%s'%i, 3), (-1.0, 'btedge%s'%k, 3),(-1.0, 'RP5', 3)))
                                    mdb.models[modelName].Equation(name='G-3-btedge-bbedge%s'%k, terms=((1.0, 'btedge%s'%k, 3), (-1.0, 'bbedge%s'%j, 3),(-1.0, 'RP6', 3)))
                                    mdb.models[modelName].Equation(name='G-3-bbedge-fbedge%s'%j, terms=((1.0, 'bbedge%s'%j, 3), (-1.0, 'fbedge%s'%l, 3),(1.0, 'RP5', 3)))


                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    mdb.models[modelName].Equation(name='G-1-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 1), (-1.0, 'bledge%s'%k, 1),(-1.0, 'RP3', 1)))
                                    mdb.models[modelName].Equation(name='G-1-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 1), (-1.0, 'bredge%s'%j, 1),(-1.0, 'RP5', 1)))
                                    mdb.models[modelName].Equation(name='G-1-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 1), (-1.0, 'fredge%s'%l, 1),(1.0, 'RP3', 1)))
                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    mdb.models[modelName].Equation(name='G-2-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 2), (-1.0, 'bledge%s'%k, 2),(-1.0, 'RP4', 2)))
                                    mdb.models[modelName].Equation(name='G-2-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 2), (-1.0, 'bredge%s'%j, 2),(-1.0, 'RP6', 2)))
                                    mdb.models[modelName].Equation(name='G-2-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 2), (-1.0, 'fredge%s'%l, 2),(1.0, 'RP4', 2)))
                                for i,k,j,l in zip(fledge,bledge,bredge,fredge):
                                    mdb.models[modelName].Equation(name='G-3-fledge-bledge%s'%i, terms=((1.0, 'fledge%s'%i, 3), (-1.0, 'bledge%s'%k, 3),(-1.0, 'RP5', 3)))
                                    mdb.models[modelName].Equation(name='G-3-bledge-bredge%s'%k, terms=((1.0, 'bledge%s'%k, 3), (-1.0, 'bredge%s'%j, 3),(-1.0, 'RP2', 3)))
                                    mdb.models[modelName].Equation(name='G-3-bredge-fredge%s'%j, terms=((1.0, 'bredge%s'%j, 3), (-1.0, 'fredge%s'%l, 3),(1.0, 'RP5', 3)))
                                    

                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    mdb.models[modelName].Equation(name='G-1-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 1), (-1.0, 'lbedge%s'%k, 1),(-1.0, 'RP4', 1)))
                                    mdb.models[modelName].Equation(name='G-1-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 1), (-1.0, 'rbedge%s'%j, 1),(-1.0, 'RP5', 1)))
                                    mdb.models[modelName].Equation(name='G-1-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 1), (-1.0, 'rtedge%s'%l, 1),(1.0, 'RP4', 1)))                                    
                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    mdb.models[modelName].Equation(name='G-2-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 2), (-1.0, 'lbedge%s'%k, 2),(-1.0, 'RP1', 2)))
                                    mdb.models[modelName].Equation(name='G-2-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 2), (-1.0, 'rbedge%s'%j, 2),(-1.0, 'RP6', 2)))
                                    mdb.models[modelName].Equation(name='G-2-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 2), (-1.0, 'rtedge%s'%l, 2),(1.0, 'RP1', 2)))
                                for i,k,j,l in zip(ltedge,lbedge,rbedge,rtedge):
                                    mdb.models[modelName].Equation(name='G-3-ltedge-lbedge%s'%i, terms=((1.0, 'ltedge%s'%i, 3), (-1.0, 'lbedge%s'%k, 3),(-1.0, 'RP6', 3)))
                                    mdb.models[modelName].Equation(name='G-3-lbtedge-rbedge%s'%k, terms=((1.0, 'lbedge%s'%k, 3), (-1.0, 'rbedge%s'%j, 3),(-1.0, 'RP2', 3)))
                                    mdb.models[modelName].Equation(name='G-3-rbedge-rtbedge%s'%j, terms=((1.0, 'rbedge%s'%j, 3), (-1.0, 'rtedge%s'%l, 3),(1.0, 'RP6', 3)))

                        ## Shear modulus G12 ##
                        if G12==True and onlyPBC == False:
                                for i in mdb.models[modelName].loads.keys():
                                        del mdb.models[modelName].loads[i]
                                for i in mdb.models[modelName].boundaryConditions.keys():
                                        del mdb.models[modelName].boundaryConditions[i]

                                region = a.sets['RP4']
                                mdb.models[modelName].DisplacementBC(name='G12-1', createStepName='Step-1', 
                                    region=region, u1=Dispx, u2=Dispy, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)


                                region = a.sets['RP5']
                                mdb.models[modelName].DisplacementBC(name='G12-2', createStepName='Step-1', 
                                    region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                region = a.sets['RP6']
                                mdb.models[modelName].DisplacementBC(name='G12-3', createStepName='Step-1', 
                                    region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                    createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                    sectionPoints=DEFAULT, rebar=EXCLUDE)

                                import os, glob

                                mdb.Job(name='job-G12', model= modelName, description='', type=ANALYSIS, 
                                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                    scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                mdb.jobs['job-G12'].submit(consistencyChecking=OFF)

                                mdb.jobs['job-G12'].waitForCompletion()


                                o3 = session.openOdb(name='%s' % (path+'\job-G12.odb'))


                                odb = session.odbs['%s' % (path+'\job-G12.odb')]

                                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                odbName=session.viewports[session.currentViewportName].odbDisplay.name


                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                    NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP4', ))

                                forceG12 = 0
                                for i in session.xyDataObjects.keys():
                                    forceG12=forceG12+(session.xyDataObjects[i][0][1])

                                stressG12 = abs(forceG12/(L*W))

                                G12 = stressG12/((Dispx/H)+(Dispy/L))

                        ## Shear modulus G13 ##                                
                        if G12==False or onlyPBC == True:
                                G12='N/A'


                        if G13==True and onlyPBC == False:
                                for i in mdb.models[modelName].loads.keys():
                                        del mdb.models[modelName].loads[i]
                                for i in mdb.models[modelName].boundaryConditions.keys():
                                        del mdb.models[modelName].boundaryConditions[i]
         

                                region = a.sets['RP5']
                                mdb.models[modelName].DisplacementBC(name='G13-1', createStepName='Step-1', 
                                    region=region, u1=Dispx, u2=UNSET, u3=Dispz, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                region = a.sets['RP4']
                                mdb.models[modelName].DisplacementBC(name='G13-2', createStepName='Step-1', 
                                    region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                region = a.sets['RP6']
                                mdb.models[modelName].DisplacementBC(name='G13-3', createStepName='Step-1', 
                                    region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)


                                regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                    createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                    sectionPoints=DEFAULT, rebar=EXCLUDE)

                                import os, glob

                                mdb.Job(name='job-G13', model= modelName, description='', type=ANALYSIS, 
                                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                    scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                mdb.jobs['job-G13'].submit(consistencyChecking=OFF)

                                mdb.jobs['job-G13'].waitForCompletion()


                                o3 = session.openOdb(name='%s' % (path+'\job-G13.odb'))


                                odb = session.odbs['%s' % (path+'\job-G13.odb')]

                                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                odbName=session.viewports[session.currentViewportName].odbDisplay.name



                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                    NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=('RP5', ))

                                forceG13 = 0
                                for i in session.xyDataObjects.keys():
                                    forceG13=forceG13+(session.xyDataObjects[i][0][1])

                                stressG13 = abs(forceG13/(H*L))


                                G13 = stressG13/((Dispx/W)+(Dispz/L))
                                
                        if G13==False or onlyPBC == True:
                                G13='N/A'

                        ## Shear modulus G23 ##
                        if G23==True and onlyPBC == False:
                                for i in mdb.models[modelName].loads.keys():
                                        del mdb.models[modelName].loads[i]
                                for i in mdb.models[modelName].boundaryConditions.keys():
                                        del mdb.models[modelName].boundaryConditions[i]
         

                                region = a.sets['RP6']
                                mdb.models[modelName].DisplacementBC(name='G23-1', createStepName='Step-1', 
                                    region=region, u1=UNSET, u2=Dispy, u3=Dispz, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                region = a.sets['RP4']
                                mdb.models[modelName].DisplacementBC(name='G23-2', createStepName='Step-1', 
                                    region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)

                                region = a.sets['RP5']
                                mdb.models[modelName].DisplacementBC(name='G23-3', createStepName='Step-1', 
                                    region=region, u1=0, u2=0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                                    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
                                    localCsys=None)


                                regionDef=mdb.models[modelName].rootAssembly.sets['c1']
                                mdb.models[modelName].HistoryOutputRequest(name='H-Output-2', 
                                    createStepName='Step-1', variables=('RT', ), region=regionDef, 
                                    sectionPoints=DEFAULT, rebar=EXCLUDE)

                                import os, glob


                                mdb.Job(name='job-G23', model= modelName, description='', type=ANALYSIS, 
                                    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                                    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                                    scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)
                                mdb.jobs['job-G23'].submit(consistencyChecking=OFF)

                                mdb.jobs['job-G23'].waitForCompletion()

                                o3 = session.openOdb(name='%s' % (path+'\job-G23.odb'))

                                odb = session.odbs['%s' % (path+'\job-G23.odb')]

                                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                                odbName=session.viewports[session.currentViewportName].odbDisplay.name

                                for i in session.xyDataObjects.keys():
                                    del session.xyDataObjects['%s' % (i)]

                                session.odbData[odbName].setValues(activeFrames=(('Step-1', (1, )), ))
                                session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
                                    NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('RP6', ))

                                forceG23 = 0
                                for i in session.xyDataObjects.keys():
                                    forceG23=forceG23+(session.xyDataObjects[i][0][1])

                                stressG23 = abs(forceG23/(L*H))


                                G23 = stressG23/((Dispy/W)+(Dispz/H))

                        if G23==False or onlyPBC == True:
                                G23='N/A'
                        
                        density = 0
                        if mass != None:
                                density = mass/(L*W*H)
                        
                        print ('----------------------------------------------------')
                        print ('----------------------------------------------------')
                        print ('The homogenised elastic properties:')
                        print ('E11=%s Stress units' % (E11))
                        print ('V12=%s ratio' % (V12))
                        print ('V13=%s ratio' % (V13))
                        print ('E22=%s Stress units' % (E22))
                        print ('V21=%s ratio' % (V21))
                        print ('V23=%s ratio' % (V23))
                        print ('E33=%s Stress units' % (E33))
                        print ('V31=%s ratio' % (V31))
                        print ('V32=%s ratio' % (V32))
                        print ('G12=%s Stress units' % (G12))
                        print ('G13=%s Stress units' % (G13))
                        print ('G23=%s Stress units' % (G23))
                        print ('----------------------------------------------------')
                        print ('Total mass=%s Mass units' % (mass))
                        print ('Homogenised density=%s Density units' % (density))
                        print ('----------------------------------------------------')
                        print ('Processing duration %s seconds' % (time.time()-start))
                        print ('----------------------------------------------------')
                        
                        filename = ('%s_elastic_properties.txt' % part)
                        print ('The homogenised elastic properties are saved in ABAQUS Work Directory under %s' % filename)
                        f = open(filename,'w')
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('Property','Value','Unit'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('E11',E11,'Stress units'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V12',V12,'ratio'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V13',V13,'ratio'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('E22',E22,'Stress units'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V21',V21,'ratio'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V23',V23,'ratio'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('E33',E33,'Stress units'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V31',V31,'ratio'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('V32',V32,'ratio'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('G12',G12,'Stress units'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('G13',G13,'Stress units'))
                        f.write('{0:^10}{1:^20}{2:^20}\n'.format('G23',G23,'Stress units'))

                        f.write ('Total mass=%s Mass units \n' % (mass))
                        f.write ('Homogenised density=%s Density units \n' % (density))
                   
                        f.write ('processing duration %s Seconds' % (time.time()-start))

                        f.close()

                        print ('Citation: Omairey S, Dunning P, Sriramula S (2018) Development of an ABAQUS plugin tool for periodic RVE homogenisation.')
                        print ('Engineering with Computers. https://doi.org/10.1007/s00366-018-0616-4')


                        filename = ('%s_elastic_properties(easycopy).txt' % part)
                        f = open(filename,'w')
                        f.write('{0:^10}\n'.format(E11))
                        f.write('{0:^10}\n'.format(E22))
                        f.write('{0:^10}\n'.format(E33))
                        f.write('{0:^10}\n'.format(G12))
                        f.write('{0:^10}\n'.format(G13))
                        f.write('{0:^10}\n'.format(G23))
                        f.write('{0:^10}\n'.format(V12))
                        f.write('{0:^10}\n'.format(V13))
                        f.write('{0:^10}\n'.format(V21))
                        f.write('{0:^10}\n'.format(V23))
                        f.write('{0:^10}\n'.format(V31))
                        f.write('{0:^10}\n'.format(V32))
                        f.write ('{0:^10}\n' .format(mass))
                        f.write ('{0:^10}\n' .format(density))                   
                        f.write ('{0:^10}\n' .format((time.time()-start)))

                        f.close()

                        print ('----------------------------------------------------')
                        if onlyPBC == True:
                                print ('EasyPBC created Period Boundary Conditions only. For further investigation, used relevant Reference Points to apply loads/displacements based on your needs. Details on the use of Reference Points are illustrated in Table 1 of the referred paper.')

                        for i in session.xyDataObjects.keys():
                            del session.xyDataObjects['%s' % (i)]
                        print ('---------------------------------------')
                        print ('--------- End of EasyPBC (3D) ---------')
                        print ('---------------------------------------')

                        if len(session.odbData.keys()) >= 1:                              
                                odb.close(odb, write=TRUE)
                                a = mdb.models[modelName].rootAssembly
                                session.viewports['Viewport: 1'].setValues(displayedObject=a)

                if error==True:
                        print ('Error(s) found during sets creation, please check the error No.(s) above with EasyPBC user guide.')
                        
                        a.SetFromNodeLabels(name='Error set', nodeLabels=((instanceName,errorset),))


