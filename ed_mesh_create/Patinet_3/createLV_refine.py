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

Laxis = np.array([0,0,1])

###################################################################

refinemesh=0.3
isepiflip = False
isendoflip = True
zoffset = 63.70 
angle = np.arccos(Laxis[2])
raxis = Laxis 
raxis = raxis / np.linalg.norm(raxis)

#volumedatadir = './'
endostr = "Endo_time_0_1_23"
epistr = "Epi_time_wo_peri_0_1_6"

#endofilename = volumedatadir + endostr +  '.stl'
endofilename =  endostr +  '.stl'
endopdata = vtk_py.readSTL(endofilename)

#epifilename = volumedatadir + epistr + '.stl'
epifilename =  epistr + '.stl'
epipdata = vtk_py.readSTL(epifilename)

rotatedepipdata = vtk_py.rotatePData_w_axis(epipdata, angle, raxis)
rotatedendopdata = vtk_py.rotatePData_w_axis(endopdata, angle, raxis) 

height =  min(rotatedendopdata.GetBounds()[5], rotatedepipdata.GetBounds()[5]) - zoffset

clipped_endo, clipped_epi = vtk_py.clipSurfacesForCutLVMesh(rotatedendopdata, rotatedepipdata, height, verbose=True)



vtk_py.writeSTL(clipped_epi, "clipped_epi_tmp.stl")
vtk_py.writeSTL(clipped_endo, "clipped_endo_tmp.stl")

filename = './case_P3'
 
vtk_py.createLVmesh(filename, refinemesh, "clipped_epi_tmp.stl", "clipped_endo_tmp.stl")

#os.remove("clipped_epi_tmp.stl")
#os.remove("clipped_endo_tmp.stl")
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



