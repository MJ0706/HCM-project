import sys
sys.path.append("/mnt/home/mojumder/HPCC_HCM")
import dolfin as dolfin
import vtk_py as vtk_py
import vtk
import os
import numpy as np
#from pyquaternion import Quaternion

######################################################################################
###################################################################

#Laxis = np.array([-0.610921976481633, +0.699597562144933, +0.370590865635734 ]) #P_1
#Laxis = np.array([-0.441013657716372, 0.810654622737842, 0.385157157977616 ]) #P_2
#Laxis = np.array([0.6174884547303, -0.711899587654135, -0.334525014569403 ]) #P_2_v2
#Laxis = np.array([-0.330766037414272, 0.883420264871849, 0.331907312524165]) #P_3
#Laxis = np.array([0.054157, -0.106292, 0.982859])
Laxis = np.array([0,0,1])

###################################################################

refinemesh=0.5
isepiflip = False
isendoflip = True
zoffset = 12.75#12.75 #8.0
angle = np.arccos(Laxis[2])
raxis = Laxis #np.cross(Laxis, np.array([0,0,1]))
raxis = raxis / np.linalg.norm(raxis)

#volumedatadir = './'
endostr = "endo_time_0_5"
epistr = "epi_time_0_8"

#endofilename = volumedatadir + endostr +  '.stl'
endofilename =  endostr +  '.stl'
endopdata = vtk_py.readSTL(endofilename)

#epifilename = volumedatadir + epistr + '.stl'
epifilename =  epistr + '.stl'
epipdata = vtk_py.readSTL(epifilename)

rotatedepipdata = vtk_py.rotatePData_w_axis(epipdata, angle, raxis)
rotatedendopdata = vtk_py.rotatePData_w_axis(endopdata, angle, raxis) 

height =  min(rotatedendopdata.GetBounds()[5], rotatedepipdata.GetBounds()[5]) - zoffset
#height =  max(rotatedendopdata.GetBounds()[5], rotatedepipdata.GetBounds()[5]) - zoffset

clipped_endo, clipped_epi = vtk_py.clipSurfacesForCutLVMesh(rotatedendopdata, rotatedepipdata, height, verbose=True)



vtk_py.writeSTL(clipped_epi, "clipped_epi_tmp.stl")
vtk_py.writeSTL(clipped_endo, "clipped_endo_tmp.stl")

filename = './case_P1'
 
vtk_py.createLVmesh(filename, refinemesh, "clipped_epi_tmp.stl", "clipped_endo_tmp.stl")

os.remove("clipped_epi_tmp.stl")
os.remove("clipped_endo_tmp.stl")
os.remove("LVtemp.geo.bak")

ugrid = vtk_py.readUGrid(filename +  '.vtk')
rotmat = vtk.vtkMatrix4x4()
rotmat.SetElement(0,0,1)
rotmat.SetElement(1,1,1)
rotmat.SetElement(2,2,1)
rotmat.SetElement(2,3,-ugrid.GetBounds()[5]*1)
rotmat.SetElement(3,3,1)

#rotmat = np.array([[0.1, 0, 0, 0], [0, 0.1, 0, 0], [0, 0, 0.1, 0], [0, 0, 0, 1]])
ugrid_transform = vtk_py.transform_mesh_w_4x4mat(filename + '.vtk', filename + '_rot.vtk', rotmat)
#ugrid_transform = vtk_py.readUGrid(filename +  '_rot.vtk')
print "ugrid=trans", ugrid_transform



