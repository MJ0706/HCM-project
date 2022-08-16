import sys
#sys.path.append("/mnt/home/mojumder/HPCC_HCM/UCSF/mesh_create")
sys.path.append("/mnt/home/mojumder/HPCC_HCM")
import os
import vtk
import vtk_py as vtk_py
import dolfin as dolfin
import numpy as np
import math as math
from dolfin import *
from mpi4py import MPI as pyMPI

def partitionmeshforEP(ugrid, nsector=6, nz=2, meas=None, center=None, xaxis=[1,0]):

    if(center == None):
    	midx = 0.5*(ugrid.GetBounds()[0] + ugrid.GetBounds()[1])
    	midy = 0.5*(ugrid.GetBounds()[2] + ugrid.GetBounds()[3])
    	center = [midx, midy]

    
    cellcenter = vtk.vtkCellCenters()
    cellcenter.SetInputData(ugrid)
    cellcenter.Update()
    
    zpartition = np.linspace(ugrid.GetBounds()[4], ugrid.GetBounds()[5], nz+1)
    apartition = np.linspace(-math.pi, math.pi, nsector+1)
    
    regid = vtk.vtkIntArray()
    data = vtk.vtkFloatArray()

    for cellid in range(0, ugrid.GetNumberOfCells()):
     	x = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[0]
     	y = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[1]
     	z = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[2]
    	
    	# Determine position in z direction
    	zloc = np.argmax(zpartition>z)
    
    	# Determine position in theta direction
    	norm = np.linalg.norm([(x - midx), (y - midy)])
    	angle = np.arctan2((y - midy)/norm, (x - midx)/norm)
    	sloc = np.argmax(apartition>angle)
    
    	regloc = (zloc)*nsector + sloc 
	#regloc = (zloc-1)*nsector + sloc 
    	if regloc <=0:
		print cellid, zloc, norm, angle, sloc, regid

    	regid.InsertNextValue(regloc)
	#data.InsertNextValue(meas[regloc-1])
    
    regid.SetName("Regionid")
    #data.SetName("EP measurements")
    
    ugrid.GetCellData().AddArray(regid)
    #ugrid.GetCellData().AddArray(data)
    
    vtk_py.writeXMLUGrid(ugrid, "test.vtu")


 
isepiflip = False #True
isendoflip = True #False

casename = "case_P1"
LVangle = 70




meshfilename = casename + "_rot.vtk"

outdir = "./" +casename + "/"
directory = os.getcwd() + '/' +  casename + "/"
ugrid1 = vtk_py.readUGrid(meshfilename)
	

#print (ugrid)
	
mesh = vtk_py.convertUGridToXMLMesh(ugrid1)

if(1):
    x = mesh.coordinates()
    scaling_factor = 0.1
    x[:, :] *= scaling_factor
        #mesh.intersection_operator().clear() # up to version 1.2.0
    mesh.bounding_box_tree().build(mesh) # development version

	
print (mesh)
ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
comm2 = pyMPI.COMM_WORLD

xmlgrid = vtk_py.convertUGridToXMLMesh(ugrid)

Qelem = dolfin.FiniteElement("Quadrature", xmlgrid.ufl_cell(), degree=4, quad_scheme="default")
segFS = dolfin.FunctionSpace(xmlgrid, Qelem)
segid = dolfin.Function(segFS)
partitionmeshforEP(ugrid, nsector=6, nz=2)

cnt = 0
for cell in dolfin.cells(xmlgrid):
        idx = int(ugrid.GetCellData().GetArray("Regionid").GetTuple(cnt)[0])
        print idx, cnt
        xmlgrid.domains().set_marker((cell.index(), idx), 3)
        cnt += 1
    
matid = dolfin.MeshFunction("size_t", xmlgrid, 3, xmlgrid.domains())

dolfin.File(casename+"_subDomain"+".pvd") << matid


xmlgrid, xmlfacet, xmledges = vtk_py.extractFeNiCsBiVFacet(ugrid, geometry="LV")

VQuadelem = dolfin.VectorElement("Quadrature", 
                              xmlgrid.ufl_cell(), 
                              degree=4, 
                              quad_scheme="default")
VQuadelem._quad_scheme = 'default'

fiberFS = dolfin.FunctionSpace(xmlgrid, VQuadelem)

ef, es, en, eC, eL, eR = vtk_py.addLVfiber(xmlgrid, fiberFS, casename, LVangle, -LVangle, [] , isepiflip, isendoflip)




f = dolfin.HDF5File(xmlgrid.mpi_comm(), directory +casename+".hdf5", 'w')
f.write(xmlgrid, casename)
f.close()

f = dolfin.HDF5File(xmlgrid.mpi_comm(), directory +casename+".hdf5", 'a') 
f.write(xmlfacet, casename+"/"+"facetboundaries") 
f.write(xmledges, casename+"/"+"edgeboundaries") 
f.write(ef, casename+"/"+"eF") 
f.write(es, casename+"/"+"eS") 
f.write(en, casename+"/"+"eN")
f.write(eC, casename+"/"+"eC") 
f.write(eL, casename+"/"+"eL") 
f.write(eR, casename+"/"+"eR") 
f.write(matid, casename+"/"+"matid")
f.close()

dolfin.File(casename+"facetboundaries.pvd")  << xmlfacet
dolfin.File(casename+"Edges.pvd")  << xmledges


