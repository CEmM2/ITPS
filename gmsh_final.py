#%%
import pygame
pygame.init()
pygame.mixer.init()
sound = pygame.mixer.Sound('C:\\temp\\cash.wav')
sound.set_volume(0.6)   
sound.play()     
#%%
import gmsh
import os
import sys
import numpy as np
import math
import csv
import meshio
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#check last simulated model

dir='C:\\temp\simulated_model.txt'
if os.path.exists(dir):
    with open(dir, 'r') as f:
        model_num=int(f.readline())+1
else:
    model_num=1

print('next run',model_num)

dir='C:/temp/model'+str(model_num)
if not os.path.exists(dir):
    os.mkdir(dir)
# #Parameter defenition




theta_f=np.zeros((2,4)) #fiber angle matrix
with open('C:/temp/geom_data.csv', newline='') as csvfile:
    data = csv.DictReader(csvfile)
    
    for i,row in enumerate(data):
        if i+1==model_num:
            iso=web=np.double(row['iso'])
            web=np.double(row['web'])
            theta_f[0,0]=np.double(row['theta(0,0)'])
            theta_f[0,1]=np.double(row['theta(0,1)'])
            theta_f[0,2]=np.double(row['theta(0,2)'])
            theta_f[0,3]=np.double(row['theta(0,3)'])
            theta_f[1,0]=np.double(row['theta(1,0)'])
            theta_f[1,1]=np.double(row['theta(1,1)'])
            theta_f[1,2]=np.double(row['theta(1,2)'])
            theta_f[1,3]=np.double(row['theta(1,3)'])
            seed_size=np.double(row['seed'])
            fibers=(row['fibers']=='TRUE')
            print('data for model ',model_num,' read from file')
            print('web',web)
            print('web',theta_f[0,0])
            print('fibers',fibers)
#%%
el_type='tri'
xlen=50
fs_s=1
fs=32
theta=60 # hexagon angle
iso_t=iso*web #isolation thicknessS
web_t=(web-iso_t)/2
web_b=iso_t+web_t
depth=0.3 #model depth
zdepth=depth
clearance=2.0*web
# seed_size=web*0.5
# seed_size=1
fiber_width=(web-iso_t)/2   #fiber width

#solution defenition-True for simulation execution
mech=True
Term=True
    
#
cot=math.cos(math.radians(theta))
sint=math.sin(math.radians(theta))
tant=math.tan(math.radians(theta))
a=(xlen-4*(web/tant))/(8*(1+cot))
xmin=0
xmax=4*(2*a*(1+cot)+web/tant)
ymin=-3*(2*a*sint+web)
ymax=(2*a*sint+web)
rect=np.zeros((4,2)) 
rect[0,:]= (xmin,ymin)
rect[1,:]= (xmin,ymax)                     
rect[2,:]= (xmax,ymax)                         
rect[3,:]= (xmax,ymin)  
##############################gmsh#
def shapley_fiber(points,rect):
    polygon = Polygon(points)
    square=Polygon(rect)
    poly=polygon.intersection(square)

    
    coordinates_array = np.asarray(poly.exterior.coords)
    
    return coordinates_array[0:-1][:]
def edge_fix(points,tanf,xmin,xmax,ymin,ymax):
    index=-1
    counter=[]
    for i in points:
        index+=1
        #left border
        if i[0]<xmin:
            i[1]+=tanf*(xmin-i[0])
            i[0]=xmin
            counter.append(index)
        if i[1]<ymin:
            i[0]+=(ymin-i[1])/tanf
            i[1]=ymin
            counter.append(index)
        if i[0]>xmax:
            
            i[1]-=tanf*(i[0]-xmax)
            i[0]=xmax
            counter.append(index)
        if i[1]>ymax:
           i[0]-=(i[1]-ymax)/tanf
           i[1]=ymax
           counter.append(index)
        if i[0]<xmin:
            i[1]+=tanf*(xmin-i[0])
            i[0]=xmin
            counter.append(index)

        if i[1]<ymin :
            i[1]=ymin
        if i[1]>ymax :
            i[1]=ymax
    # points=np.delete(points, counter, axis=0)     
    # print('indexes changed ',counter)
    counter1=[]
    
    for i in range(len(points[:,0])):
        for j in range(i+1,len(points[:,0])):
             
            if (round((points[i,0]-points[j,0]),2)==0 and round((points[i,1]-points[j,1]),2)==0):
                counter1.append(j)
    # print('deleted',counter1,'rows')                   
    points = np.delete(points, counter1, axis=0) 
    return points
def poly_surf(point_list):
    #creating points
    length=len(point_list[:,1])
    pl=[]
    ln=[]

    for i in(point_list):
        point_num=gmsh.model.occ.addPoint(i[0], i[1],0, lc)
        pl.append(point_num)
    # print('point list',pl)    
    #creating lines
    for j in range(length-1):
        line_num=gmsh.model.occ.addLine(pl[j],pl[j+1])
        ln.append(line_num)
    line_num=gmsh.model.occ.addLine(pl[length-1],pl[0])
    ln.append(line_num)
    # print('lines',ln)
    #creating curve loop
    cl=gmsh.model.occ.addCurveLoop(ln)
    # print('curve loop index',cl)
    # creating plane surface
    gmsh.model.occ.addPlaneSurface([cl])

def mirror_cell(surf1,surf2,p_x,p_y,p_z,rot_x,rot_y):
    if rot_x == True:
        #rotation-X
        gmsh.model.occ.rotate(surf1,p_x,p_y,p_z,1,0,0,-math.pi)
        gmsh.model.occ.rotate(surf2,p_x,p_y,p_z,1,0,0,-math.pi)
        
    elif rot_y == True:
        #rotation-Y
        gmsh.model.occ.rotate(surf1,p_x,p_y,p_z,0,1,0,-math.pi)
        gmsh.model.occ.rotate(surf2,p_x,p_y,p_z,0,1,0,-math.pi)
        # gmsh.model.occ.rotate(surf3,p_x,p_y,p_z,0,1,0,-math.pi)
        # print(surf1[0][1]) 
        # print(surf2[0][1]) 
        # print(surf3[0][1]) 
    
    return surf1,surf2
def extra_cut(xmin_in,ymin_in,zmin_in,xmax_in,ymax_in,zmax_in,xmin_out,ymin_out,zmin_out,xmax_out,ymax_out,zmax_out):

    inner = gmsh.model.getEntitiesInBoundingBox(xmin_in, ymin_in,
                                                    zmin_in,xmax_in,
                                                    ymax_in, zmax_in, 2)
    outer = gmsh.model.getEntitiesInBoundingBox(xmin_out, ymin_out,
                                                    zmin_out,xmax_out,
                                                    ymax_out, zmax_out, 2)
    print('inner surfaces')
    print(len(inner))
    print('all surfaces')
    print(len(outer))
    diff=set(outer) - set(inner)
    print(len(diff))
    for i in diff:
        print(i)
        gmsh.model.occ.remove([i],True)
def periodic(xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,eps):
    symin = gmsh.model.getEntitiesInBoundingBox(xmin, ymin, zmin, xmax,
                                            ymax,zmax, 1)
    pn=[]
    translation = [1, 0, 0, dx, 0, 1, 0, dy, 0, 0, 1, 0, 0, 0, 0, 1]
    counter=0
    for i in symin:
        # Then we get the bounding box of each left surface
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(i[0], i[1])
        # We translate the bounding box to the right and look for surfaces inside
        # it:
        symax = gmsh.model.getEntitiesInBoundingBox(xmin - eps + dx, ymin - eps,
                                                    zmin - eps, xmax + eps + dx,
                                                    ymax + eps+dy, zmax + eps, 1)
        # For all the matches, we compare the corresponding bounding boxes...
        for j in symax:
            xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
                j[0], j[1])
            
            # ...and if they match, we apply the periodicity constraint
            if (abs(xmin2 - xmin-dx) < eps and abs(xmax2 - xmax-dx) < eps
                    and abs(ymin2 - ymin-dy) < eps and abs(ymax2 - ymax-dy) < eps
                    and abs(zmin2 - zmin) < eps and abs(zmax2 - zmax) < eps):
                counter+=1
                gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], translation)
                pn.append(gmsh.model.mesh.get_periodic_nodes(1,i[1]))
    
    print('number of paired edges -',counter)
    return pn



lc=seed_size
eps=1e-5
gmsh.initialize()

gmsh.model.add("Mesh_model")

#Quarter cell defenition
#containing rect
points=np.zeros((4,2))
points[0,:]=[0,-3*( web+2*a*sint)]
points[1,:]=[0,web+2*a*sint]
points[2,:]=[4*(2*a*(1+cot)+web/tant),web+2*a*sint]
points[3,:]=[4*(2*a*(1+cot)+web/tant),-3*( web+2*a*sint)]
poly_surf(points)

#bottom web

points=np.zeros((8,2))
points[0,:]=[0,0]
points[1,:]=[0,web_t]
points[2,:]=[a+web_b/tant,web_t]
points[3,:]=[a*(1+2*cot)+web_b/tant,2*a*sint+web_t]
points[4,:]=[2*a*(1+cot)+web/tant,2*a*sint+web_t]
points[5,:]=[2*a*(1+cot)+web/tant,2*a*sint]
points[6,:]=[a*(1+2*cot)+web/tant,2*a*sint]
points[7,:]=[a+web/tant,0]
poly_surf(points)

# top web
points=np.zeros((8,2))
points[0,:]=[0,web_t+iso_t]
points[1,:]=[0,web]
points[2,:]=[a,web]
points[3,:]=[a+2*a*cot,web+2*a*sint]
points[4,:]=[2*a*(1+cot)+web/tant,2*a*sint+web]
points[5,:]=[2*a*(1+cot)+web/tant,2*a*sint+web_t+iso_t]
points[6,:]=[a+web_t/tant+2*a*cot,2*a*sint+web_t+iso_t]
points[7,:]=[a+web_t/tant,web_t+iso_t]
poly_surf(points)





surf1,surf2=mirror_cell(gmsh.model.occ.copy([(2, 2)]),gmsh.model.occ.copy([(2, 3)]),0,0,0,True,False) #mirror bottom
surf3,surf4=mirror_cell(gmsh.model.occ.copy(surf1),gmsh.model.occ.copy(surf2),2*a*(1+cot)+web/tant,0,0,False,True) #mirror right
surf5,surf6=mirror_cell(gmsh.model.occ.copy(surf3),gmsh.model.occ.copy(surf4),2*a*(1+cot)+web/tant,0,0,True,False)
#######mirror right all cell
surf7,surf8=mirror_cell(gmsh.model.occ.copy([(2, 2)]),gmsh.model.occ.copy([(2, 3)]),2*(2*a*(1+cot)+web/tant),0,0,False,True) #mirror bottom
surf9,surf10=mirror_cell(gmsh.model.occ.copy(surf1),gmsh.model.occ.copy(surf2),2*(2*a*(1+cot)+web/tant),0,0,False,True) #mirror bottom
surf11,surf12=mirror_cell(gmsh.model.occ.copy(surf5),gmsh.model.occ.copy(surf6),2*(2*a*(1+cot)+web/tant),0,0,False,True) #mirror bottom
surf13,surf14=mirror_cell(gmsh.model.occ.copy(surf3),gmsh.model.occ.copy(surf4),2*(2*a*(1+cot)+web/tant),0,0,False,True) #mirror bottom
#mirror bottom all
surf15,surf16=mirror_cell(gmsh.model.occ.copy([(2, 2)]),gmsh.model.occ.copy([(2, 3)]),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
surf17,surf18=mirror_cell(gmsh.model.occ.copy(surf1),gmsh.model.occ.copy(surf2),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
surf19,surf20=mirror_cell(gmsh.model.occ.copy(surf3),gmsh.model.occ.copy(surf4),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
surf21,surf22=mirror_cell(gmsh.model.occ.copy(surf5),gmsh.model.occ.copy(surf6),0,-1*(2*a*sint+web),0,True,False) #mirror bottom

surf23,surf24=mirror_cell(gmsh.model.occ.copy(surf7),gmsh.model.occ.copy(surf8),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
surf25,surf26=mirror_cell(gmsh.model.occ.copy(surf9),gmsh.model.occ.copy(surf10),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
surf27,surf28=mirror_cell(gmsh.model.occ.copy(surf11),gmsh.model.occ.copy(surf12),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
surf29,surf30=mirror_cell(gmsh.model.occ.copy(surf13),gmsh.model.occ.copy(surf14),0,-1*(2*a*sint+web),0,True,False) #mirror bottom
#fibers building
#fibers


#finiding the crossing coords between a cirvle with radius r located at the center of cell and fibers
def circle_sol(tanf,x_shift,r,center_x,center_y):
    sol_x=np.roots((1+tanf**2,2*x_shift*tanf**2,tanf**2*x_shift**2-r**2))
    sol_y=tanf*(sol_x+x_shift)
    sol_x=sol_x+center_x
    sol_y=sol_y+center_y
    p1=(sol_x[0],sol_y[0])
    p2=(sol_x[1],sol_y[1])

    if(np.linalg.norm(p1)>=np.linalg.norm(p2)):
        temp=p1
        p1=p2
        p2=temp
            
    return p1,p2

# fixing the fibers so they are rectengular for easier meshing 
# trimming the long fiber to be equal to the smaller one
def rect_fibers(p1,p2,p3,p4,theta_f):
    p1=np.array(p1)
    p2=np.array(p2)
    p3=np.array(p3)
    p4=np.array(p4)

    if(np.linalg.norm(p1-p2)<=np.linalg.norm(p3-p4)):
        #p1p2 is smaller fiber 
        add=fiber_width*np.array([np.cos(np.radians(90-theta_f)),-np.sin(np.radians(90-theta_f))])
        p4=p1+add
        p3=p2+add
    else :
        #p3p4 is smaller fiber
        add=fiber_width*np.array([np.cos(np.radians(90-theta_f)),-np.sin(np.radians(90-theta_f))])
        p1=p4-add
        p2=p3-add
    
    return p1,p2,p3,p4
#this function finds points location for edge cells for section assigment return all points and crossing divided by quarters

#parameter defenition
if fibers ==False:
    sides=0
else :
    sides=2

r=2*a*sint*0.95
b=0
center_x=2*(2*a*(1+cot)+web/tant)
center_y=(2*a*sint+web)


y_shift=0

n_fibers=50
center_x=0
center_y=-2*a*sint-web
xcoord=[1]

#running threw normal middle cells 
#range(1,4)
for l in range(1,4):
    center_x=0+(2*a*(1+cot)+web/tant)*l
    for s in range(2):
        if(s==1 and l==2):
            break
        if (l == 2):
            s=1
        cell=np.array([s,l])
        center_y=-2*a*sint-web+2*(2*a*sint+web)*s-np.remainder(l,2)*(2*a*sint+web)
        if(l==2):
            center_y-=2*(2*a*sint+web)
        print('cell num'+str(cell))
        tanf=math.tan(math.radians(theta_f[s,l]))
        m=tanf
        
        for j in range(sides):
            #building fibers in 2 steps - to the right and then to the left
            k=0
            if(j ==1): # fiber num starts from 1 and shift is to the oposite direction
                clearance=-clearance
                k=1
            for i in range(k,n_fibers): #running thru fibers until no crossing point was found between inner circle and fiber
                #first edge
                if (theta_f[s,l]==0): 
                    x_shift=0
                else:
                    x_shift=(fiber_width/2)/np.sin(np.radians(theta_f[s,l]))+(i*clearance)/np.sin(np.radians(theta_f[s,l]))
                p1,p2=circle_sol(tanf,x_shift,r,center_x,center_y)       
                
                #second edge
                if (theta_f[s,l]==0): 
                    x_shift=0
                else:

                    x_shift=-(fiber_width/2)/np.sin(np.radians(theta_f[s,l]))+(i*clearance)/np.sin(np.radians(theta_f[s,l]))
                y_shift=0
            
                p4,p3=circle_sol(tanf,x_shift,r,center_x,center_y) 
                #if no crossing was found - break
                if (np.all(np.isreal(np.array([p1,p2,p3,p4]))) == False):
                    break 
                #making fiber rectengular
                p1,p2,p3,p4=rect_fibers(p1,p2,p3,p4,theta_f[s,l])
                
                points=np.zeros((4,2))
                points[0,:]=p1
                points[1,:]=p2
                points[2,:]=p3
                points[3,:]=p4
                poly_surf(points)
                fs+=1
        print ('fs after i',s,'j',l,'cell',fs)          

#making (0,0) cell-cell is composed of 4 corners - 4 domains

center_x=0
center_y=-2*a*sint-web
tanf=math.tan(math.radians(theta_f[0,0]))
print('cell [0,0] ')
cell=np.array([0,0])
for q in range(2):
    for s in range(2):
        
        center_x=0+s*4*(2*a*(1+cot)+web/tant)
        center_y=-2*a*sint-web-q*(8*a*sint+4*web)+2*(2*a*sint+web)
        #matrix section-before the fibers creation for easier assigment
        


        
        for j in range(sides):
                    
                    
                    #building fibers in 2 steps - to the right and then to the left
                    k=0
                    if(j ==1):# fiber num starts from 1 and shift is to the oposite direction
                        clearance=-clearance
                        k=1
                    for i in range(k,n_fibers):#running thru fibers until no crossing point was found between inner circle and fiber
                        #first edge
                        x_shift=(fiber_width/2)/np.sin(np.radians(theta_f[0,0]))+(i*clearance)/np.sin(np.radians(theta_f[0,0]))
                        p1,p2=circle_sol(tanf,x_shift,r,center_x,center_y)       
                        
                        #second edge
                        x_shift=-(fiber_width/2)/np.sin(np.radians(theta_f[0,0]))+(i*clearance)/np.sin(np.radians(theta_f[0,0]))
                        y_shift=0
                    
                        p4,p3=circle_sol(tanf,x_shift,r,center_x,center_y) 
                        #if no crossing was found - break
                        if (np.all(np.isreal(np.array([p1,p2,p3,p4]))) == False):
                            break 
                        #making fiber rectengular
                        p1,p2,p3,p4=rect_fibers(p1,p2,p3,p4,theta_f[0,0])
                        x_shift=(i*clearance)/np.sin(np.radians(theta_f[0,0]))
                        
                        
                        
                        
                        
                        
                        #sketching fiber
                        # fiber_sketcher(np.array([p1,p2,p3,p4]),sketch_face,sketch_edge,1)
                        points=np.zeros((4,2))
                        points[0,:]=p1
                        points[1,:]=p2
                        points[2,:]=p3
                        points[3,:]=p4

                        xmin=0
                        ymin=-3*(2*a*sint+web)
                    
                        xmax=4*(2*a*(1+cot)+web/tant)
                        ymax=(2*a*sint+web)
                        pt=points
                        p=shapley_fiber(pt,rect)
                        points= edge_fix(points,tanf,xmin,xmax,ymin,ymax)
                        
                        # print(len(points[:,1]))
                        # if (k==0 and len(points[:,1])==4):
                        #     points= np.insert(points,2, np.array((center_x,center_y)), 0)
                            
                        if len(points[:,0])<3 :
                            print('breaking')
                            break                            
                        # if i ==0:
                        #     np.append(points,[center_x,center_y],axis=0)
                        # poly_surf(points)
                        poly_surf(p)
                        fs+=1
#making (1,0 cell)-coposef of 2 halfes split by vertical line
center_x=0
center_y=-2*a*sint-web
tanf=math.tan(math.radians(theta_f[1,0]))
print('cell [1,0] ')
cell=np.array([1,0])
for s in range(2):
    center_x=0+s*4*(2*a*(1+cot)+web/tant)
    #matrix section-before the fibers creation for easier assigment
    
    for j in range(sides):
                
                
                #building fibers in 2 steps - to the right and then to the left
                k=0
                if(j ==1):# fiber num starts from 1 and shift is to the oposite direction
                    clearance=-clearance
                    k=1
                for i in range(k,n_fibers):#running thru fibers until no crossing point was found between inner circle and fiber
                    #first edge
                    x_shift=(fiber_width/2)/np.sin(np.radians(theta_f[1,0]))+(i*clearance)/np.sin(np.radians(theta_f[1,0]))
                    p1,p2=circle_sol(tanf,x_shift,r,center_x,center_y)       
                    
                    #second edge
                    x_shift=-(fiber_width/2)/np.sin(np.radians(theta_f[1,0]))+(i*clearance)/np.sin(np.radians(theta_f[1,0]))
                    y_shift=0
                
                    p4,p3=circle_sol(tanf,x_shift,r,center_x,center_y) 
                    #if no crossing was found - break
                    if (np.all(np.isreal(np.array([p1,p2,p3,p4]))) == False):
                        break 
                    #making fiber rectengular
                    p1,p2,p3,p4=rect_fibers(p1,p2,p3,p4,theta_f[1,0])
                    x_shift=(i*clearance)/np.sin(np.radians(theta_f[1,0]))
                    points=np.zeros((4,2))
                    points[0,:]=p1
                    points[1,:]=p2
                    points[2,:]=p3
                    points[3,:]=p4


                   
                    # xmin=0
                    # ymin=-3*(2*a*sint+web)
                    
                    # xmax=4*(2*a*(1+cot)+web/tant)
                    # ymax=(2*a*sint+web)
                    pt=points
                    p=shapley_fiber(pt,rect)
                    points= edge_fix(points,tanf,xmin,xmax,ymin,ymax)
                        
                    if len(points[:,0])<3 :
                        print('breaking')
                        break
                    poly_surf(p)
                    fs+=1
                   
                    
                   

#making (0,2 cell)- cell containing to halfes splitted by horizontal line
center_x=2*(2*a*(1+cot)+web/tant)
center_y=-6*a*sint-3*web
tanf=math.tan(math.radians(theta_f[0,2]))
print('cell [0,2] ')
cell=np.array([1,0])
for s in range(2):
    #matrix section-before the fibers creation for easier assigment
    center_y=center_y+s*(8*a*sint+4*web)
    
    for j in range(sides):
        
                #building fibers in 2 steps - to the right and then to the left
                k=0
                if(j ==1):# fiber num starts from 1 and shift is to the oposite direction
                    clearance=-clearance
                    k=1
                for i in range(k,n_fibers):#running thru fibers until no crossing point was found between inner circle and fiber
                    #first edge
                    x_shift=(fiber_width/2)/np.sin(np.radians(theta_f[0,2]))+(i*clearance)/np.sin(np.radians(theta_f[0,2]))
                    p1,p2=circle_sol(tanf,x_shift,r,center_x,center_y)       
                    
                    #second edge
                    x_shift=-(fiber_width/2)/np.sin(np.radians(theta_f[0,2]))+(i*clearance)/np.sin(np.radians(theta_f[0,2]))
                    y_shift=0
                
                    p4,p3=circle_sol(tanf,x_shift,r,center_x,center_y) 
                
                    if (np.all(np.isreal(np.array([p1,p2,p3,p4]))) == False):
                        break 
                    
                    
                    p1,p2,p3,p4=rect_fibers(p1,p2,p3,p4,theta_f[0,2])
                    x_shift=(i*clearance)/np.sin(np.radians(theta_f[0,2]))
                    points=np.zeros((4,2))
                    points[0,:]=p1
                    points[1,:]=p2
                    points[2,:]=p3
                    points[3,:]=p4
                    
                   
                  
                    # xmin=0
                    # ymin=-3*(2*a*sint+web)
                    
                    # xmax=4*(2*a*(1+cot)+web/tant)
                    # ymax=(2*a*sint+web)
                    pt=points
                    p=shapley_fiber(pt,rect)
                    points= edge_fix(points,tanf,xmin,xmax,ymin,ymax)
                    
                    if len(points[:,0])<3 :
                        # print('breaking')
                        break   
                    print('points after edge fix=',points)    
                    poly_surf(p)
                    fs+=1

#cutting using fiberous material sections
dimtags=gmsh.model.occ.get_entities(dim=2)

gmsh.model.occ.fragment([(2,1)], [(2, i) for i in range(2, len(dimtags)+1)])


gmsh.model.occ.synchronize()


xmin=-eps
xmax=eps
ymin=-3*(2*a*sint+web)-eps
ymax=(2*a*sint+web)+eps
zmin=-eps
zmax=eps
dx=4*(2*a*(1+cot)+web/tant)
dy=0
rln=periodic(xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,eps)
#top bottom sides

xmin=-eps
xmax=4*(2*a*(1+cot)+web/tant)+eps
ymin=-3*(2*a*sint+web)-eps
ymax=-3*(2*a*sint+web)+eps
zmin=-eps
zmax=eps
dx=0
dy=4*(2*a*sint+web)
tbn=periodic(xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,eps)

####extrude all surfaces in z direction
h=depth
dimtags=gmsh.model.occ.get_entities(dim=2)

gmsh.model.occ.synchronize()
if (el_type=='quad'):
    gmsh.option.setNumber("Mesh.RecombineAll", 3)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    for i in zip(dimtags):
        gmsh.model.occ.extrude([i[0]], 0, 0, h, [1],recombine=True)
        
else:
    for i in zip(dimtags):
        gmsh.model.occ.extrude([i[0]], 0, 0, h, [1])

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.model.mesh.remove_duplicate_nodes()
allnodes=gmsh.model.mesh.get_nodes()
allele=gmsh.model.mesh.get_elements()
n_vol=gmsh.model.occ.get_entities(dim=3)

print('nele',(allele[0]))
print('n_nodes',len(allnodes[1]))
gmsh.write("Mesh_model_2.inp")


# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()


# #
gmsh.finalize()
vol,v=np.shape(n_vol)
print(np.shape(n_vol))
print(vol)

import numpy as np
import os

#%%
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
    
    print('line numbers',line_numbers)
    counter = 0
  
    with open(original_file, 'r') as read_obj, open(new_file, 'w+') as write_obj:
        for line in read_obj:
            
            if counter in line_numbers:
                print('counter',counter)
                write_obj.write(line)
            
            counter+=1
# %%

file_name='C:\\temp\Mesh_model_2.inp'
string_to_search='T3D2'
ans1=search_string_in_file(file_name, string_to_search)
if (el_type == 'quad'):
    string_to_search='C3D8'
else:
    string_to_search='C3D4'
ans2=search_string_in_file(file_name, string_to_search)
print('first index',ans1[0][0]-1,'last index',ans2[0][0])
delete_multiple_lines(file_name, np.arange(ans1[0][0]-1,ans2[0][0]-1))

# %% copy paste 

# file_name='C:\\temp\Job-1.inp'
# string_to_search='type=T3D2'
# ans1=search_string_in_file(file_name, string_to_search)

# string_to_search='type=C3D4'
# ans2=search_string_in_file(file_name, string_to_search)

# delete_multiple_lines(file_name, np.arange(ans1[0][0]-1,ans2[0][0]))
# %%


#read node data from file to connect nodes for pbc
# %%
# # print('all nodes 1',len(allnodes[1]))
# l=len(allnodes[1])
# all=np.array(allnodes[1])
# all=all.reshape(int(l/3),3)
# # print('all',all[1,1])
# index_arr=np.arange(1,int(l/3+1))
# index_arr=index_arr.reshape(int(l/3),1)
# all_nodes=np.append(all,index_arr,axis=1)
# # print('all w index',all)
# # %%


    
# #first step-pairing x dir nodes-edges negx_negz with posx_negz , edges negx_posz with posx_posz 
    
# xmin=0
# xmax=4*(2*a*(1+cot)+web/tant)
# ymin=-3*(2*a*sint+web)
# ymax=-3*(2*a*sint+web)
# zmin=0
# zmax=h

# nodes_x=[]
# nodes_y=[]
# nodes_z=[]
# dx=xmax-xmin
# dy=ymax-ymin
# dz=zmax-zmin
# face_neg_x=np.asarray(np.where(all_nodes[:,0]==xmin))

# face_neg_x=face_neg_x.flatten()
# print('neg x face',face_neg_x)
# print('nodes',all_nodes[face_neg_x,2])
# edge_negx_negz=[]
# edge_negx_posz=[]
# for i in face_neg_x:
    
#     if abs(all_nodes[i,2]-zmin)< eps:
#         edge_negx_negz.append([i])
#     elif abs(all_nodes[i,2]-zmax)< eps:
#         edge_negx_posz.append([i])   
# edge_negx_negz=np.int8(np.where(all_nodes[face_neg_x,2]==zmin))
# edge_negx_posz=np.int8(np.where(all_nodes[face_neg_x,2]==zmax))
# print(len(face_neg_x))
# print(len(edge_negx_negz))

# print('pos nodes',all_nodes[edge_negx_posz,:])

# face_posx=np.where(all_nodes[:,0]==xmax)
# edge_posx_negz=np.int8(np.where(all_nodes[face_posx,2]==zmin))
# edge_posx_posz=np.int8(np.where(all_nodes[face_posx,2]==zmax))
# face_neg_y=np.where(all_nodes[:,1]==ymin)
# edge_negy_negz=np.int8(np.where(all_nodes[face_neg_y,2]==zmin))
# edge_negy_posz=np.int8(np.where(all_nodes[face_neg_y,2]==zmax))
# face_pos_y=np.where(all_nodes[:,1]==ymax)
# edge_posy_negz=np.int8(np.where(all_nodes[face_neg_y,2]==zmin))
# edge_posy_posz=np.int8(np.where(all_nodes[face_neg_y,2]==zmax))
# #nodeslist
# print('index 1',len(edge_negx_negz))
# print('index 2',len(edge_negx_posz))

# e1=np.int8(edge_negx_negz[1,:])
# e2=np.int8(edge_posx_negz[1,:])
# print('check',all_nodes[e1,:])
# edge1=all_nodes[e1,:]
# edge2=all_nodes[e2,:]
# print('edge1',edge1)
# print('edge2',edge2)

# columnIndex = 1
# # # Sort 2D numpy array by 2nd Column
# edge1 =edge1[edge1[:,columnIndex].argsort()]
# edge2 =edge2[edge2[:,columnIndex].argsort()]
# print('edge1',edge1)
# print('edge2',edge2)

# find_periodic(all,xmin,xmax,ymin,ymax,zmin,zmax)

# %%
# f = open("Mesh_model_2.inp", "r")
# count=0
# for x in f:
#   count+=1
#   if(count>3 and count<10):
#     print(x)
#   if(x =='******* E L E M E N T S *************'):
#       break

xmin=0
xmax=4*(2*a*(1+cot)+web/tant)
ymin=-3*(2*a*sint+web)
ymax=(2*a*sint+web)
dx=xmax-xmin
dy=ymax-ymin
dz=h
# print(n_vol[-1][1],'',vol)
# print(str(vol)+'\n')
# print(fs+1)
# # %%
if (str(vol)+'\n'!=str(n_vol[-1][1])+'\n' ):
    print('problem')
    trans=vol-99
    fs+=trans
    fs_s=1+trans
with open('C:\\temp\section_data.txt','w') as d:
    d.write(str(fs_s) +'\n')
    d.write(str(fs+1) +'\n')
    d.write(str(vol)+'\n')
    d.write(str(dx)+'\n')
    d.write(str(dy)+'\n')
    d.write(str(dz)+'\n')
    d.write(str(model_num)+'\n')
    d.write('web='+str(web)+'\n')
    d.write('a='+str(a)+'\n')
    d.write('iso_t='+str(iso_t)+'\n')
    d.write('Theta_f='+str(theta_f[0,:])+'  '+str(theta_f[1,:])+'\n')
    d.write('fiber width='+str(fiber_width)+'\n')
    d.write('mesh size='+str(lc))
    d.write


print('finished meshing model num',model_num)
#%%
# import shutil
# from subprocess import run
# tes=shutil.which("abaqus")
# print(tes)
# #run('C:\Windows\system32\cmd.exe /k "C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Dassault Systemes SIMULIA Abaqus CAE 2018\ABC.lnk" ' ,check=True,shell=True)
# run(tes+" cae noGUI=mech_sim.py", shell=True)
# print("finished calculation")
# %%
# import os
# os.path.isfile('abaqus')


# %%
