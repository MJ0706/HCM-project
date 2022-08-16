# Partitioning Biventricular mesh for visualizing EP measurements
# This script partition the mesh into polar segments:
# Input: ugrid = vtk file containing BiVMesh
#      : center = origin of the polar coordinate system (if none, mesh centroid will be used)
#      : xaxis = direction corresponding to 0 degree in the polar coord system
#      : nsectors = number of circular partition
#      : nz = number of longitudinal partition
#      : meas = measurements of size nz*nsectors

import sys, pdb
sys.path.append("/mnt/home/mojumder/HPCC_HCM")
import vtk_py as vtk_py
import vtk
import numpy as np
import math as math
from dolfin import *
import pdb

def add_EpiSegments2hdf5(inPut, outPut, fname, caseName):
    
    #pdb.set_trace()

    mesh = Mesh()
    meshname = caseName
    f = HDF5File(mpi_comm_world(), inPut + fname+".hdf5", 'r') 
    f.read(mesh, meshname, False)

    facetboundaries = MeshFunction("size_t", mesh, 2)
    f.read(facetboundaries, meshname+"/"+"facetboundaries")

    matid = MeshFunction('size_t', mesh, 3)
    f.read(matid, meshname+"/"+"matid") 

    f.close()

    ugridLoc = vtk_py.convertXMLMeshToUGrid(mesh)

    nSec = 10; nzs = 4
    EPdelay = np.arange(nSec*nzs)*1
    ugridLoc = partitionmeshforEP(ugridLoc, nsector=nSec, nz=nzs, meas=EPdelay)

    outF = outPut
    casename = caseName

    AHAid = MeshFunction('size_t', mesh, 3)
    AHAid = defSubDomain_AHA(mesh = mesh, ugrid = ugridLoc, matId = matid)

    EpiBCid = FacetFunction('size_t', mesh)
    EpiBCid = defSubDomain_SurfaceSegment(mesh = mesh, ugrid = ugridLoc, matId = AHAid, 
                                 facetBdry = EpiBCid, epiMarks = facetboundaries)
    dolfin.File(outF+"/" + casename+"_epiFacets"+".pvd") << EpiBCid


    f = HDF5File(mpi_comm_world(), outPut+ fname+".hdf5", 'a') 
    #f.write(mesh, meshname)

    f.write(EpiBCid, casename+"/"+"EpiBCid_Corr")

    f.close()


def add_Segments2hdf5(inPut, outPut, fname, caseName):
    
    #pdb.set_trace()

    mesh = Mesh()
    meshname = caseName
    f = HDF5File(mpi_comm_world(), inPut + fname+".hdf5", 'r') 
    f.read(mesh, meshname, False)

    facetboundaries = MeshFunction("size_t", mesh, 2)
    f.read(facetboundaries, meshname+"/"+"facetboundaries")

    matid = MeshFunction('size_t', mesh, 3)
    f.read(matid, meshname+"/"+"matid") 

    f.close()

    ugridLoc = vtk_py.convertXMLMeshToUGrid(mesh)

    nSec = 10; nzs = 4
    EPdelay = np.arange(nSec*nzs)*1
    ugridLoc = partitionmeshforEP(ugridLoc, nsector=nSec, nz=nzs, meas=EPdelay)

    outF = outPut
    casename = caseName

    partId = MeshFunction('size_t', mesh, 3)
    partId = defSubDomain1(mesh = mesh, ugrid = ugridLoc)
    dolfin.File(outF+"/" + casename+"_partition"+".pvd") << partId

    PurKid = FacetFunction('size_t', mesh)
    PurKid = defSubDomain_Purk_3(mesh = mesh, ugrid = ugridLoc, matId = partId, 
                                 facetBdry = PurKid, epiMarks = facetboundaries)
    dolfin.File(outF+"/" + casename+"_Pukrfacets"+".pvd") << PurKid

    AHAid = MeshFunction('size_t', mesh, 3)
    AHAid = defSubDomain_AHA(mesh = mesh, ugrid = ugridLoc, matId = matid)
    dolfin.File(outF+"/" + casename+"_AHAsegmentation"+".pvd") << AHAid

    EpiBCid = FacetFunction('size_t', mesh)
    EpiBCid = defSubDomain_SurfaceSegment(mesh = mesh, ugrid = ugridLoc, matId = AHAid, 
                                 facetBdry = EpiBCid, epiMarks = facetboundaries)
    dolfin.File(outF+"/" + casename+"_epiFacets"+".pvd") << EpiBCid



    #f = HDF5File(mpi_comm_world(), outPut+ fname+".hdf5", 'w') 
    #f.write(mesh, meshname)
    #f.write(facetboundaries, casename+"/"+"facetboundaries")
    #f.write(fenics_edge_ref, casename+"/"+"edgeboundaries") 
    #f.write(matid, casename+"/"+"matid")

    f = HDF5File(mpi_comm_world(), outPut+ fname+".hdf5", 'a') 
    #f.write(mesh, meshname)

    f.write(partId, casename+"/"+"partId")
    f.write(PurKid, casename+"/"+"PurKid")
    f.write(EpiBCid, casename+"/"+"EpiBCid")
    f.write(AHAid, casename+"/"+"AHAid")

    f.close()
    

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
    
    	regloc = (zloc-1)*nsector + sloc 
    
    	regid.InsertNextValue(regloc)
	data.InsertNextValue(meas[regloc-1])
    
    regid.SetName("Regionid")
    data.SetName("EP measurements")
    
    ugrid.GetCellData().AddArray(regid)
    ugrid.GetCellData().AddArray(data)
    
    #vtk_py.writeXMLUGrid(ugrid, "./test.vtu")
    return ugrid

	
class K(Expression):
    def __init__(self, k_0, k_1, materials, **kwargs):
        self.materials = materials
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, ufc_cell):
        if self.materials[ufc_cell.index] == 7:
            values[0] = self.k_0
        else:
            values[0] = self.k_1

def defMaterialProperty1(dolfin_mesh, mId):
    matProp = K(12.0, 0.0, materials=mId, degree=0)
    dolfin.File(meshname+"_matProp"+".pvd") << interpolate(matProp, FunctionSpace(dolfin_mesh,'DG', 0))

def defMaterialProperty2(dolfin_mesh, mId):
    matProp = K(12.0, 0.0, materials=mId, degree=0)
    dolfin.File(meshname+"_matProp"+".pvd") << interpolate(matProp, FunctionSpace(dolfin_mesh,'DG', 0))

def defSubDomain1(mesh, ugrid, meshname = "CRT27_AS_smooth_fine"):
    subDomain1 = mesh.domains()

    cnt = 0
    for cell in cells(mesh):
        idx = int(ugrid.GetCellData().GetArray("Regionid").GetTuple(cnt)[0])
        #print idx, cnt
        subDomain1.set_marker((cell.index(), idx), 3)
        cnt += 1
    
    matId1 = MeshFunction("size_t", mesh, 3, subDomain1)
    #dolfin.File(meshname+"_partitionId1"+".pvd") << matId1
    #hfd5Fille.write(matId1, meshname+"/matId1")
    return matId1


def defSubDomain_AHA(mesh, ugrid, matId):
    subDomain1 = mesh.domains()

    V = FunctionSpace(mesh, "DG", 0)
    dm = V.dofmap()

    cnt = 0
    for cell in cells(mesh):
        idx = int(ugrid.GetCellData().GetArray("Regionid").GetTuple(cnt)[0])
        #mId = int(matId.array()[cnt-1]) # the DOFs are jumbled
        # https://fenicsproject.org/qa/1680/proper-way-to-turn-meshfunction-over-cells-into-a-function/ 
        mId2 = int(matId.array()[dm.cell_dofs(cell.index())]) 
        #print  cnt, idx, mId2
        #pdb.set_trace()
        
        # apical slice
        if idx in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] and mId2 in [0, 1]: 
            subDomain1.set_marker((cell.index(), 17), 3)

        elif idx in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] and mId2 in [2]: # RV
            subDomain1.set_marker((cell.index(), 20), 3)

        # mid apical slice
        elif idx in [11, 12, 13, 14] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 13), 3)
        elif idx in [15, 16] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 16), 3)
        elif idx in [17, 18] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 15), 3)
        elif idx in [19, 20] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 14), 3)

        # mid apical slice RV
        elif idx in [11, 12, 13, 14] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 21), 3)
        elif idx in [15, 16] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 21), 3) #21?
        elif idx in [17, 18] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 22), 3)
        elif idx in [19] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 22), 3)
        elif idx in [20] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 21), 3)

        # mid slice
        elif idx in [21, 22] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 8), 3)
        elif idx in [23, 24] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 7), 3)
        elif idx in [25] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 12), 3)
        elif idx in [26] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 11), 3)
        elif idx in [27, 28] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 10), 3)
        elif idx in [29, 30] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 9), 3)

        # mid slice RV
        elif idx in [21, 22] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 23), 3)
        elif idx in [23, 24] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 23), 3)
        #elif idx in [25] and mId2 in [2]:
        #    subDomain1.set_marker((cell.index(), 12), 3)
        #elif idx in [26] and mId2 in [2]:
        #    subDomain1.set_marker((cell.index(), 11), 3)
        elif idx in [27, 28] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 24), 3)
        elif idx in [29, 30] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 24), 3)

        # basal slice
        elif idx in [31, 32] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 2), 3)
        elif idx in [33, 34] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 1), 3)
        elif idx in [35] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 6), 3)
        elif idx in [36] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 5), 3)
        elif idx in [37, 38] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 4), 3)
        elif idx in [39, 40] and mId2 in [0, 1]:
            subDomain1.set_marker((cell.index(), 3), 3)

        # basal slice RV 
        elif idx in [31, 32] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 25), 3)
        elif idx in [33, 34] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 25), 3)
        #elif idx in [35] and mId2 in [2]:
        #    subDomain1.set_marker((cell.index(), 6), 3)
        #elif idx in [36] and mId2 in [2]:
        #    subDomain1.set_marker((cell.index(), 5), 3)
        elif idx in [37, 38] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 26), 3)
        elif idx in [39, 40] and mId2 in [2]:
            subDomain1.set_marker((cell.index(), 26), 3)

        #elif mId2 in [2]:
        #    subDomain1.set_marker((cell.index(), 0), 3)

        # bug 
        else: 
            subDomain1.set_marker((cell.index(), 99), 3)

        cnt += 1
    
    AHAId1 = MeshFunction("size_t", mesh, 3, subDomain1)
    #dolfin.File(meshname+"_partitionId1"+".pvd") << matId1
    #hfd5Fille.write(matId1, meshname+"/matId1")
    return AHAId1
        

def defSubDomain_Purk_3(mesh, ugrid, matId, facetBdry, epiMarks):

    purk_facets = dolfin.FacetFunction('size_t', mesh)
    purk_facets.set_all(0)

    V = FunctionSpace(mesh, "DG", 0)
    dm = V.dofmap()

    facetBdry.set_all(0)
    for cell in cells(mesh):
        cellID = cell.index()
        #print [facet.index() for facet in facets(cell)]
        mId2 = int(matId.array()[dm.cell_dofs(cellID)])
        for facet in facets(cell):
            #facetBdry[facet.index()] = mId2
            #https://fenicsproject.org/qa/6186/facet-exterior-always-false/
            if facet.exterior():
                epiID = epiMarks[facet.index()]
                #print epiID
                if epiID == 2: #3 # LV, supposed 2 b
                    facetBdry[facet.index()] = mId2

    return facetBdry

def defSubDomain_SurfaceSegment(mesh, ugrid, matId, facetBdry, epiMarks):

    purk_facets = dolfin.FacetFunction('size_t', mesh)
    purk_facets.set_all(0)

    V = FunctionSpace(mesh, "DG", 0)
    dm = V.dofmap()

    facetBdry.set_all(0)
    for cell in cells(mesh):
        cellID = cell.index()
        #print [facet.index() for facet in facets(cell)]
        mId2 = int(matId.array()[dm.cell_dofs(cellID)])
        for facet in facets(cell):
            if facet.exterior():
                epiID = epiMarks[facet.index()]
                #print epiID
                if epiID == 1: # EPI, supposed 2 b
                    facetBdry[facet.index()] = mId2
                elif epiID == 2: #3 # LV, supposed 2 b
                    facetBdry[facet.index()] = mId2 + 50 
                elif epiID == 3: #3 # RV, supposed 2 b
                    facetBdry[facet.index()] = mId2 + 100 
                elif epiID == 4: #3 # RV, supposed 2 b
                    facetBdry[facet.index()] = mId2 + 150 

    return facetBdry


def defSubDomain_EpiDelays(mesh, ugrid, matId, facetBdry, epiMarks):

    purk_facets = dolfin.FacetFunction('size_t', mesh)
    purk_facets.set_all(0)

    V = FunctionSpace(mesh, "DG", 0)
    dm = V.dofmap()

    # https://fenicsproject.org/qa/1741/tetrahedrons-list-mesh-facets-and-node-for-all-tetrahedron/
    #print max( [facet.index() for facet in facets(mesh)])
    #print max( [face.index() for face in faces(mesh)])
    #pdb.set_trace()

    '''
    for facet in facets(mesh):
        #cell_ID = facet.index() 
        #mId2 = int(matId.array()[dm.cell_dofs(cell_ID)]) 
        mId2 = int(matId.array()[cell_ID]) 
        print cell_ID, mId2
    '''
    facetBdry.set_all(0)
    for cell in cells(mesh):
        cellID = cell.index()
        #print [facet.index() for facet in facets(cell)]
        mId2 = int(matId.array()[dm.cell_dofs(cellID)])
        for facet in facets(cell):
            #facetBdry[facet.index()] = mId2
            #https://fenicsproject.org/qa/6186/facet-exterior-always-false/
            if facet.exterior():
                epiID = epiMarks[facet.index()]
                #print epiID
                if epiID == 1: # EPI, supposed 2 b
                    facetBdry[facet.index()] = mId2

    return facetBdry


def defSubDomain2(mesh, ugrid, meshname = "CRT27_AS_smooth_fine"):
    subDomain2 = mesh.domains()

    cnt = 0
    for cell in cells(mesh):
        idx = int(ugrid.GetCellData().GetArray("Regionid").GetTuple(cnt)[0])
        #print idx, cnt
        subDomain2.set_marker((cell.index(), idx + 100), 3)
        cnt += 1
    
    meshname = 'CRT27_AS_smooth_fine'
    matId2 = MeshFunction("size_t", mesh, 3, subDomain2)
    dolfin.File(meshname+"_partitionId2"+".pvd") << matId2
    #hfd5Fille.write(matId2, meshname+"/matId2")

'''
def getCPP_MatProp_string(mId_List):

    retString = 'True'
    for iD in mId_List:
        print iD
        print '(*materials)[cell.index] == %s' % iD
        adDString = '(*materials)[cell.index] == %s' % iD
        retString = retString + adDString
    return retString 
'''

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
def defCPP_Matprop(mesh, mId, k):

    #conditionalString = getCPP_MatProp_string(mId_List)

    cppcode = """
    class K : public Expression
    {
    public:

      void eval(Array<double>& values,
            const Array<double>& x,
            const ufc::cell& cell) const
      {
        if ( (*materials)[cell.index] == 5 || (*materials)[cell.index] == 6 )
          values[0] = k_0;
        else
          values[0] = k_1;
      }

      std::shared_ptr<MeshFunction<std::size_t>> materials;
      double k_0;
      double k_1;

    };
    """

    #pdb.set_trace()

    kappa = Expression(cppcode=cppcode, degree=0)
    kappa.materials = mId
    kappa.k_0 = k["kabNormal"] #0.1*146.e3
    kappa.k_1 = k["kNormal"]#146.e3

    #print getCPP_MatProp_string(mId_List)
    print k["kabNormal"], k["kNormal"]

    return kappa

    #dolfin.File(meshname+"_matProp2_Ischemia"+".pvd") <<  interpolate(kappa, FunctionSpace(mesh,'DG', 0))
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -




#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
def defCPP_Matprop_DIsch(mesh, mId, k):

    #conditionalString = getCPP_MatProp_string(mId_List)

    cppcode = """
    class K : public Expression
    {
    public:

      void eval(Array<double>& values,
            const Array<double>& x,
            const ufc::cell& cell) const
      {
        if ( (*materials)[cell.index] == 5 || (*materials)[cell.index] == 6 )
          values[0] = k_0;
        else
          values[0] = k_1;
      }

      std::shared_ptr<MeshFunction<std::size_t>> materials;
      double k_0;
      double k_1;

    };
    """

    #pdb.set_trace()

    kappa = Expression(cppcode=cppcode, degree=0)
    kappa.materials = mId
    kappa.k_0 = k["kIschemia"] 
    kappa.k_1 = k["kNormal"]

    #print getCPP_MatProp_string(mId_List)
    print k["kIschemia"], k["kNormal"]

    return kappa

    #dolfin.File(meshname+"_matProp2_Ischemia"+".pvd") <<  interpolate(kappa, FunctionSpace(mesh,'DG', 0))
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -




def defCPP_LBBB_Matprop(mesh, mId, meshname = "CRT27_AS_smooth_fine"):
    cppcode = """
    class K : public Expression
    {
    public:

      void eval(Array<double>& values,
            const Array<double>& x,
            const ufc::cell& cell) const
      {
        if ((*materials)[cell.index] == 14 || (*materials)[cell.index] == 15 || (*materials)[cell.index] == 24 || (*materials)[cell.index] == 25 || (*materials)[cell.index] == 26 )
          values[0] = k_0;
        else
          values[0] = k_1;
      }

      std::shared_ptr<MeshFunction<std::size_t>> materials;
      double k_0;
      double k_1;

    };
    """

    kappa = Expression(cppcode=cppcode, degree=0)
    kappa.materials = mId
    kappa.k_0 = 12.0
    kappa.k_1 = 1.0

    dolfin.File(meshname+"_matProp2_LBBB"+".pvd") <<  interpolate(kappa, FunctionSpace(mesh,'DG', 0))


def defCPP_Ischemia_Matprop(mesh, mId, meshname = "CRT27_AS_smooth_fine"):
    cppcode = """
    class K : public Expression
    {
    public:

      void eval(Array<double>& values,
            const Array<double>& x,
            const ufc::cell& cell) const
      {
        if ((*materials)[cell.index] <= 10 )
          values[0] = k_0;
        else
          values[0] = k_1;
      }

      std::shared_ptr<MeshFunction<std::size_t>> materials;
      double k_0;
      double k_1;

    };
    """

    kappa = Expression(cppcode=cppcode, degree=0)
    kappa.materials = mId
    kappa.k_0 = 12.0
    kappa.k_1 = 1.0

    dolfin.File(meshname+"_matProp2_Ischemia"+".pvd") <<  interpolate(kappa, FunctionSpace(mesh,'DG', 0))


def checkHDF5Partition(meshname = 'CRT27_AS_smooth_fine'):

    dolfin_mesh = Mesh()
    
    f = HDF5File(mpi_comm_world(), meshname+".hdf5", 'r') 
    f.read(dolfin_mesh, meshname, False)
    f.close()

    ugridLoc = vtk_py.convertXMLMeshToUGrid(dolfin_mesh)

    nSec = 10
    nzs = 3
    EPdelay = np.arange(nSec*nzs)*1
    ugridLoc = partitionmeshforEP(ugridLoc, nsector=nSec, nz=nzs, meas=EPdelay)

    matId = MeshFunction('size_t', dolfin_mesh, 3)
    matId = defSubDomain1(mesh = dolfin_mesh, ugrid = ugridLoc)
    defCPP_LBBB_Matprop(mesh = dolfin_mesh, mId = matId)
    defCPP_Ischemia_Matprop(mesh = dolfin_mesh, mId = matId)

    f = HDF5File(dolfin_mesh.mpi_comm(), meshname+".hdf5", 'a')
    f.write(matId, meshname+"/"+"partId")
    f.close()

    #defSubDomain2(mesh = dolfin_mesh, ugrid = ugridLoc)
    #defMaterialProperty1(dolfin_mesh, mId = matId)

		




