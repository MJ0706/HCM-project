from dolfin import * 

from edgetypebc import *
import numpy as np
import os as os

import pdb
import vtk_py as vtk_py
#from oops_ConstantDefinitions import Constant_definitions 
#from oops_printout import printout 

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
def printout(statement, mpi_comm): 
    if(MPI.rank(mpi_comm) == 0): 
        print statement
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
def update_mesh(mesh, displacement, boundaries):
    # https://fenicsproject.org/qa/13470/ale-move-class-and-meshes
    new_mesh = Mesh(mesh)
    new_boundaries = MeshFunction("size_t", new_mesh, 2)
    new_boundaries.set_values(boundaries.array())
    ALE.move(new_mesh, displacement)
    return new_mesh, new_boundaries
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


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
        if ((*materials)[cell.index] == 14 || (*materials)[cell.index] == 15 || (*materials)[cell.index] == 16 || (*materials)[cell.index] == 24 || (*materials)[cell.index] == 25 || (*materials)[cell.index] == 26 || (*materials)[cell.index] == 27 )
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
    kappa.k_0 = 0.0
    kappa.k_1 = 1.0

    return kappa

    #dolfin.File(meshname+"_matProp2_LBBB"+".pvd") <<  interpolate(kappa, FunctionSpace(mesh,'DG', 0))
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
class Constant_definitions(object):
    '''
    - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
    Constant definitions for the Nash Panfilov problem
    The default constants in this class has to be changed in order to change the constants.

    The point is to keep constant definitions away from the sciprt and solver 
    So options to change params has NOT been given 
    - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
    '''
    def __init__(self):
        self.parameters = self.default_parameters()		
		
    def default_parameters(self):

        d1 = 0.2 # diffusion constant for action potential 0.2 in Nash Panfilov 
        return {"alpha" : Constant(0.01), # 0.1 in 'a' in NashPanfilov, 0.01 in Gok_Kuhl 
                "gamma" : Constant(0.002), # 0.01 'epsilon' in NasPanlov, 0.002 in Gok_Kuhl  
                "b" : Constant(0.15), # 0.1 in Gok_Kuhl 
                "c" : Constant(8), # 'k' is NashPanfilov, 8 in Gok_Kuhl
                "mu1" : Constant(0.2), # 0.12 in NP ,# 0.2 in Gok_Kuhl 
                "mu2" : Constant(0.3), 
                "D1" : Constant(((d1,'0.0', '0.0'),('0.0',d1, '0.0'), ('0.0', '0.0', d1))), 
                "B" : Constant((0.0, 0.0, 0.0)), 
                "T" : Constant((0.0, 0.0, 0.0))
                }
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 



#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class biventricle_mesh(object): 
    """
    object for biventricle mesh
    input: mesh, facet, edge, matid, fibre file
    output: mesh 
    """

    def default_parameters(self):
        return {"directory" : "../CRT27/", 
                "casename" : "CRT27", 
                "fibre_quad_degree" : 4, 
                "outputfolder" : "../Outputs/",
                "topid": 4,
                "LVendoid": 2,
                "RVendoid": 3,
                "epiid": 1,
}

    def update_parameters(self, params):
        self.parameters.update(params) 

    def __init__(self, params, SimDet):
        
        self.mesh = Mesh() 
        self.parameters = self.default_parameters() 
        self.parameters.update(params) 

        directory = self.parameters["directory"]
        casename = self.parameters["casename"]
        outputfolder = self.parameters["outputfolder"]
        folderName = self.parameters["foldername"]

        meshfilename = directory + casename + ".hdf5"
        f = HDF5File(mpi_comm_world(), meshfilename, 'r')
        f.read(self.mesh, casename, False)

        self.facetboundaries = MeshFunction("size_t", self.mesh, 2)
        f.read(self.facetboundaries, casename+"/"+"facetboundaries")

	self.edgeboundaries = MeshFunction("size_t", self.mesh, 1)
        f.read(self.edgeboundaries, casename+"/"+"edgeboundaries")
        
        
        deg = self.parameters["fibre_quad_degree"]
        VQuadelem = VectorElement("Quadrature", 
                                  self.mesh.ufl_cell(), 
                                  degree=deg, 
                                  quad_scheme="default")
        VQuadelem._quad_scheme = 'default'
        
        self.fiberFS = FunctionSpace(self.mesh, VQuadelem)

        self.f0 = Function(self.fiberFS)
        self.s0 = Function(self.fiberFS)
        self.n0 = Function(self.fiberFS)

        #f.read(self.f0, casename+"/"+"eF")
        #f.read(self.s0, casename+"/"+"eS")
        #f.read(self.n0, casename+"/"+"eN")


        if SimDet["DTI_ME"] is True :
            f.read(self.f0, casename+"/"+"eF_proj_DTI")
            f.read(self.s0, casename+"/"+"eS_proj_DTI")
            f.read(self.n0, casename+"/"+"eN_proj_DTI")
        else: 
            f.read(self.f0, casename+"/"+"eF")
            f.read(self.s0, casename+"/"+"eS")
            f.read(self.n0, casename+"/"+"eN")

        #self.displacement = Function(VectorFunctionSpace(self.mesh,'CG',1))
        #f.read(self.displacement, casename+"/disp/vector_0")

	self.f0 = self.f0/sqrt(inner(self.f0, self.f0))
	self.s0 = self.s0/sqrt(inner(self.s0, self.s0))
	self.n0 = self.n0/sqrt(inner(self.n0, self.n0))


	self.matid = CellFunction('size_t', self.mesh)
    	if(f.has_dataset(casename+"/"+"matid")):
        	f.read(self.matid, casename+"/"+"matid")
	else:
		self.matid.set_all(0)


	self.AHAid = CellFunction('size_t', self.mesh)
    	if(f.has_dataset(casename+"/"+"AHAid")):
        	f.read(self.AHAid, casename+"/"+"AHAid")
	else:
		self.AHAid.set_all(0)

        EpiBCid = FacetFunction('size_t', self.mesh)
    	if(f.has_dataset(casename+"/"+"EpiBCid_Corr")):
        	f.read(EpiBCid, casename+"/"+"EpiBCid_Corr")
	else:
		EpiBCid.set_all(0)

        self.EpiBCid_me = EpiBCid

        f.close()

        self.topid = self.parameters["topid"]
        self.LVendoid = self.parameters["LVendoid"]
        self.RVendoid = self.parameters["RVendoid"]
        self.epiid = self.parameters["epiid"]

	
        dx = dolfin.dx(self.mesh, subdomain_data=self.matid)
        #dx = dolfin.dx(self.mesh, subdomain_data=self.AHAid)
        ds = dolfin.ds(self.mesh, subdomain_data=self.facetboundaries)
        self.dx = dx
        self.ds = ds

        print 'Mesh size is : %f '%(self.mesh.num_cells())
        r_qmin, r_qmax = MeshQuality.radius_ratio_min_max(self.mesh)
        print('Minimal radius ratio:', r_qmin)
        print('Maximal radius ratio:', r_qmax)

        '''
        # Show histogram using matplotlib
        hist = MeshQuality.radius_ratio_matplotlib_histogram(mesh)
        hist = hist.replace('    import matplotlib.pylab', '    import matplotlib\n    matplotlib.use(\'Agg\')\n    import matplotlib.pylab\n')
        hist = hist.replace('matplotlib.pylab.show()', 'matplotlib.pylab.savefig("mesh-quality.pdf")')
        print(hist)
        #exec(hist)
        '''
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class lv_mesh(object): 
    """
    object for lv mesh
    input: mesh, facet, edge, matid, fibre file
    output: mesh 
    """

    def default_parameters(self):
        return {"directory" : "../CRT27/", 
                "casename" : "CRT27", 
                "fibre_quad_degree" : 4, 
                "outputfolder" : "../Outputs/",
                "topid": 4,
                "LVendoid": 3,
                "epiid": 2,
}

    def update_parameters(self, params):
        self.parameters.update(params) 

    def __init__(self, params, SimDet):
        
        self.mesh = Mesh() 
        self.parameters = self.default_parameters() 
        self.parameters.update(params) 

        directory = self.parameters["directory"]
        casename = self.parameters["casename"]
        outputfolder = self.parameters["outputfolder"]
        #comm_common = self.parameters["common_communicator"]
        
        folderName = self.parameters["foldername"]

        meshfilename = directory + casename + ".hdf5"
        f = HDF5File(mpi_comm_world(), meshfilename, 'r')
        f.read(self.mesh, casename, False)

        self.facetboundaries = MeshFunction("size_t", self.mesh, 2)
        f.read(self.facetboundaries, casename+"/"+"facetboundaries")

	self.edgeboundaries = MeshFunction("size_t", self.mesh, 1)
        f.read(self.edgeboundaries, casename+"/"+"edgeboundaries")
        
        deg = self.parameters["fibre_quad_degree"]
        VQuadelem = VectorElement("Quadrature", 
                                  self.mesh.ufl_cell(), 
                                  degree=deg, 
                                  quad_scheme="default")
        VQuadelem._quad_scheme = 'default'
        
        self.fiberFS = FunctionSpace(self.mesh, VQuadelem)

        self.f0 = Function(self.fiberFS)
        self.s0 = Function(self.fiberFS)
        self.n0 = Function(self.fiberFS)
	
	f00 = Function(self.fiberFS)
	s00 = Function(self.fiberFS)
	n00 = Function(self.fiberFS)

        if SimDet["DTI_ME"] is True :
            f.read(self.f0, casename+"/"+"eF_proj_DTI")
            f.read(self.s0, casename+"/"+"eS_proj_DTI")
            f.read(self.n0, casename+"/"+"eN_proj_DTI")
        else: 
            f.read(self.f0, casename+"/"+"eF")
            f.read(self.s0, casename+"/"+"eS")
            f.read(self.n0, casename+"/"+"eN")
		
	self.f0_ori = self.f0
	self.s0_ori = self.s0
	self.n0_ori = self.n0
	#print self.f0.vector().array()

	File(outputfolder + folderName +"/Mesh/mesh.pvd") << self.mesh
	File(outputfolder + folderName +"/Mesh/facetbound.pvd") << self.facetboundaries
	File(outputfolder + folderName +"/Mesh/edgebound.pvd") << self.edgeboundaries
	outfolder = outputfolder + folderName +"Mesh/"
	outdirectory = ""
	vtkoutfile = outfolder+"_e0_fiber"
	vtk_py.convertQuadDataToVTK(self.mesh, self.fiberFS, self.f0, vtkoutfile, outdirectory)
	
	#self.f0 = project(f00, VectorFunctionSpace(self.mesh, "CG", 1))
	#self.n0 = project(n00, VectorFunctionSpace(self.mesh, "CG", 1))
	#self.s0 = project(s00, VectorFunctionSpace(self.mesh, "CG", 1))

	self.f0 = self.f0/sqrt(inner(self.f0, self.f0))
	self.s0 = self.s0/sqrt(inner(self.s0, self.s0))
	self.n0 = self.n0/sqrt(inner(self.n0, self.n0))

	self.matid = CellFunction('size_t', self.mesh)
    	if(f.has_dataset(casename+"/"+"matid")):
        	f.read(self.matid, casename+"/"+"matid")
	else:
		self.matid.set_all(0)

	self.AHAid = CellFunction('size_t', self.mesh)
    	if(f.has_dataset(casename+"/"+"AHAid")):
        	f.read(self.AHAid, casename+"/"+"AHAid")
	else:
		self.AHAid.set_all(0)

        EpiBCid = FacetFunction('size_t', self.mesh)
    	if(f.has_dataset(casename+"/"+"EpiBCid_Corr")):
        	f.read(EpiBCid, casename+"/"+"EpiBCid_Corr")
	else:
		EpiBCid.set_all(0)

        self.EpiBCid_me = EpiBCid

        f.close()

        self.topid = self.parameters["topid"]
        self.LVendoid = self.parameters["LVendoid"]
        self.epiid = self.parameters["epiid"]

	
        dx = dolfin.dx(self.mesh, subdomain_data=self.matid)
        #dx = dolfin.dx(self.mesh, subdomain_data=self.AHAid)
        ds = dolfin.ds(self.mesh, subdomain_data=self.facetboundaries)
        self.dx = dx
        self.ds = ds

        print 'Mesh size is : %f '%(self.mesh.num_cells())
        r_qmin, r_qmax = MeshQuality.radius_ratio_min_max(self.mesh)
        print('Minimal radius ratio:', r_qmin)
        print('Maximal radius ratio:', r_qmax)
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class DiffusingMedium(object): 
    """
    Object for anisotropic diffuion in FHN, and Monodomain equation
    Having a seperate object allows to add effects of Purkinjee and all that 

    Input: anisotropic and isotropic diffusion coefficients could be updated
    Myofibre angles 
    
    Output: Diffusion tensor (degree and family depends on the myofibre vectors) 
    """

    def __init__(self, params, mesh, mId, isLBBB):
        self.parameters = self.default_parameters()
        self.parameters.update(params)
        self.mesh = mesh
        self.mId = mId
        self.isLBBB = isLBBB

    def default_parameters(self):
        return {"D_iso" : Constant((('0.2','0.0', '0.0'),('0.0','0.2', '0.0'), ('0.0', '0.0', '0.2')))}

    def Dmat(self):
        f0 = self.parameters["fiber"]
	s0 = self.parameters["sheet"]
	n0 = self.parameters["sheet-normal"] 

        #D_tensor = self.parameters["Dmat"]

        D_iso = self.parameters["D_iso"]
        d_ani = self.parameters["d_ani"]        
        
        if self.isLBBB == True:
            d_ani = defCPP_LBBB_Matprop(mesh = self.mesh, mId = self.mId)

        Dij = d_ani*f0[i]*f0[j] + D_iso[i,j]        
        D_tensor = as_tensor(Dij, (i,j)) 

        return D_tensor
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class PV_Elas(object):
    """
    Elasticity equations object
    """

    def default_parameters(self):
        return {"outputfolder" : "../Outputs/",
                "foldername": "NashPanfilov_BiV_04"}

    def update_parameters(self, params):
        self.parameters.update(params) 
    
    #def __init__(self, biVMesh, params, phi_ref, phi_space):
    def __init__(self, biVMesh, params, SimDet):

        self.parameters = self.default_parameters() 
        self.parameters.update(params)
        
        self.deg = self.parameters["degree"]
        deg = self.deg

        self.isincomp = self.parameters["is_incompressible"]
        isincomp = self.isincomp

        mesh = biVMesh.mesh
    	self.isLV = SimDet["isLV"]

        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        #V = VectorFunctionSpace(mesh, 'CG', 2)
        Q = FunctionSpace(mesh,'CG',1)

        Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
        #Velem = VectorElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
        #Velem._quad_scheme = 'default'
        Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
        Qelem._quad_scheme = 'default'
        Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
        Relem._quad_scheme = 'default'
        Quadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=deg, quad_scheme="default")
        Quadelem._quad_scheme = 'default'
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -



        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        Telem2 = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg, shape=2*(3,), quad_scheme='default')
        Telem2._quad_scheme = 'default'
        for e in Telem2.sub_elements():
            e._quad_scheme = 'default'
        Telem4 = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg, shape=4*(3,), quad_scheme='default')
        Telem4._quad_scheme = 'default'
        for e in Telem4.sub_elements():
            e._quad_scheme = 'default'
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


        # Mixed Element for rigid body motion 
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -



        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        if(isincomp): 
	    if(self.isLV):
            	W = FunctionSpace(mesh, MixedElement([Velem, Qelem, Relem, VRelem]))
	    else:
            	W = FunctionSpace(mesh, MixedElement([Velem, Qelem, Relem, Relem, VRelem]))
            #W = FunctionSpace(mesh, MixedElement([(Velem+Belem), Qelem, Relem, Relem, VRelem]))
        else:
	    if(self.isLV):
            	W = FunctionSpace(mesh, MixedElement([Velem, Relem, VRelem]))
	    else:
            	W = FunctionSpace(mesh, MixedElement([Velem, Relem, Relem, VRelem]))

        Quad = FunctionSpace(mesh, Quadelem)

        TF = FunctionSpace(mesh, Telem2)

        self.W = W
        self.Q = Q

        self.TF = TF

        # this displacement is for FHN 
        #W_n = Function(W)
        #Ve = VectorFunctionSpace(mesh, 'CG', 2)
        #self.we_n = Function(Ve)
        self.we_n = Function(W.sub(0).collapse())

        #phi_me = Function(Q)
        #self.phi_me = phi_me
        #self.phi_ref = phi_ref
        #self.phi_ref_space = phi_space

    '''
    def interpolate_potential_ep2me(self, v_epMesh):

        lp = LagrangeInterpolator()
        #lp.interpolate(V_me, v_epMesh)
        lp.interpolate(self.phi_ref, v_epMesh)

    def tranfer_potential_vertexValues_ref2me(self):
    
        phi_ref_array = self.phi_ref.vector().array()
        phi_me_array = self.phi_me.vector().array()
        phi_me_new_array = phi_me_array
    
        d2v_me = dof_to_vertex_map(self.Q)
        d2v_ref = dof_to_vertex_map(self.phi_ref_space)

        print d2v_me, d2v_ref

        for idx, (p_me, p_ref) in enumerate(zip(phi_me_array, phi_ref_array)):
            phi_me_new_array[d2v_me[idx]] = phi_ref_array[d2v_ref[idx]]

        self.phi_me.vector()[:] = phi_me_new_array

    def tranfer_potential_ref2me(self):
    
        #print len(self.phi_me.vector()[:]), len(self.phi_ref.vector().array())
        phi_ref_array = self.phi_ref.vector().array()
        phi_me_array = self.phi_me.vector().array()
        phi_me_new_array = phi_me_array

        # zip acts as a set filter?
        # https://www.programiz.com/python-programming/methods/built-in/zip
        # No .... 
        for idx, (p_me, p_ref) in enumerate(zip(phi_me_array, phi_ref_array)):
            phi_me_new_array[idx] = phi_ref_array[idx]

        self.phi_me.vector()[:] = phi_me_new_array

        #self.phi_ref = phi_ref
        #self.phi_me = phi_me
    '''

    def get_Jn_Cn(self):
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
        # This is for operator split based coupling 
        
        we_n = self.we_n 
        #we_n = u        

        d_n = we_n.ufl_domain().geometric_dimension()
        I_n = Identity(d_n)
        F_n = I_n + grad(we_n)
        C_n = F_n.T*F_n

        Ic_n = tr(C_n)
        J_n = det(F_n)

        return (J_n, C_n)
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

    def get_Fn(self):
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
        # This is for operator split based coupling 
        
        we_n = self.we_n 
        #we_n = u        

        d_n = we_n.ufl_domain().geometric_dimension()
        I_n = Identity(d_n)
        F_n = I_n + grad(we_n)

        return F_n
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


    def set_BCs(self, bivMesh_o, SimDet):
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        # Using bubble element 
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        #baseconstraint = project(Expression(("0.0"), degree=2), W.sub(0).sub(2).collapse())
        #bctop = DirichletBC(W.sub(0).sub(2), baseconstraint, facetboundaries, topid)
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

        facetboundaries = bivMesh_o.facetboundaries
        edgeboundaries = bivMesh_o.edgeboundaries
        topid = bivMesh_o.topid 
        W = self.W
        
        bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries, topid)

        endoring = pick_endoring_bc(method="cpp")(edgeboundaries, 1)

        bcedge = DirichletBC(W.sub(0), Expression(("0.0", "0.0", "0.0"), degree = 0), endoring, method="pointwise")
    	if("springbc" in SimDet.keys() and SimDet["springbc"]):
        	bcs = [bctop]
	else:
        	bcs = [bctop]#, bcedge]
        #bcs = [] # changing top constraint

        return bcs 
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -



#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 

class FHN(object):
    """
    FHN equations object
    
    Essentially three variables
    phi : action potential, normalized, PDE
    r   : inhibitory variable, ODE 
    Ta  : active stress, ODE
    
    Better to have a seperate object for active stress Ta  
    """

    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -     
    def default_parameters(self):
        return {"outputfolder" : "../Outputs/",
                "foldername": "NashPanfilov_BiV_04"}
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 


    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def update_parameters(self, params):
        self.parameters.update(params) 
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 


    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def __init__(self, biVMesh, params):
        
        self.parameters = self.default_parameters() 
        self.parameters.update(params)

        mesh = biVMesh.mesh
        dx = biVMesh.dx
        ds = biVMesh.ds
        
        P1_fhn = FiniteElement('CG', mesh.ufl_cell(), 1, quad_scheme="default")
        P1_fhn._quad_scheme = 'default'

        W_fhn = FunctionSpace(mesh, MixedElement([P1_fhn, P1_fhn, P1_fhn]))
        
        # FHN test and trial functions 
        w_fhn = Function(W_fhn)
        dw_fhn = TrialFunction(W_fhn)
        wtest_fhn = TestFunction(W_fhn)
        
        self.W_fhn = W_fhn
        self.w_fhn = w_fhn
        self.dw_fhn = dw_fhn
        self.wtest_fhn = wtest_fhn

        w_n_fhn = self.set_ICs() 

        self.w_n_fhn = w_n_fhn
        phi_n, r_n, Ta_n_fhn = split(self.w_n_fhn)
        phi, r, Ta_fhn  = split(w_fhn)

        bcs_FHN = self.set_BCs()

        self.bcs_FHN = bcs_FHN 

        self.phi = phi 
        self.r = r 
        self.Ta_fhn = Ta_fhn 

        self.phi_n = phi_n 
        self.r_n = r_n 
        self.Ta_n_fhn = Ta_n_fhn 

        phi_test, r_test, Ta_test_fhn = split(wtest_fhn) 
        self.phi_test = phi_test 
        self.r_test = r_test 
        self.Ta_test_fhn = Ta_test_fhn 

        gradPhi_testFunc = VectorFunctionSpace(mesh, "CG", 1)
        gradPhi_test = Function(gradPhi_testFunc)
        self.gradPhi_test = gradPhi_test

        folderName = self.parameters["foldername"]
        outputfolder = self.parameters["outputfolder"]
        caseID = self.parameters["caseID"]

        vtkfile_phi  = File(outputfolder + folderName + caseID + 'solution_phi.pvd')
        vtkfile_r    = File(outputfolder + folderName + caseID + 'solution_r.pvd')
        vtkfile_Ta_fhn = File(outputfolder + folderName + caseID + 'solution_Ta_fhn.pvd')

        vtkfile_q = File(outputfolder + folderName + caseID + 'solution_q.pvd')
        fdataECG = open(outputfolder + folderName + caseID + 'IntegralCharge_.txt', "w", 0)
        self.fdataECG = fdataECG

        
        #fdataECG2 = open(outputfolder + folderName + caseID + 'IntegralCharge2_.txt', "w", 0)
        #self.fdataECG2 = fdataECG2

        self.vtkfile_phi = vtkfile_phi
        self.vtkfile_r = vtkfile_r
        self.vtkfile_Ta_fhn = vtkfile_Ta_fhn

        self.Dmat = self.parameters["Diffusion_Tensor"]

        stimPoint = self.parameters["Stimulus_Point"]
        #stimuLusString = '( (x[0] > -0.2 + {x0} ) && (x[0] < 0.2 + {x0}) && (x[1] > -0.2 + {x1}) && (x[1] < 0.2 + {x1}) && (x[2] > - 0.2 + {x2}) && (x[2] < 0.2 + {x2}) ) ? 0.3*isStim : 0.0'.format(x0= stimPoint[0], x1= stimPoint[1], x2= stimPoint[2])

        stimCondition = self.parameters["Stimulus_Condition"]
	
        if stimCondition == 'apex':
            stimuLusString = '(x[2] <= -6.44120193470301 ) ? 1*isStim : 0.0'
        elif stimCondition == 'LVFW':
            stimuLusString = '(x[0] >= 19.199910386985) ? 1*isStim : 0.0'
        else:
            #stimuLusString = '(x[0] >= 13.3) && (x[0] <= 13.6) && (x[1] >= 15.4) && (x[1] <= 15.8) && (x[2] >= -0.1) ?  30.0*isStim : 0.0'          
            stimuLusString = '( (x[0] > -0.3 + {x0} ) && (x[0] < 0.3 + {x0}) && (x[1] > -0.3 + {x1}) && (x[1] < 0.3 + {x1}) && (x[2] > - 0.3 + {x2}) && (x[2] < 0.3 + {x2}) ) ? 10.0*isStim : 0.0'.format(x0= stimPoint[0], x1= stimPoint[1], x2= stimPoint[2])

        print stimuLusString
        #pdb.set_trace()
        stimuLusExpression = stimuLusString
        #stimuLusExpression = self.parameters["Stimulus_Expression"]
        self.stimuLus = Expression(stimuLusExpression, degree = 1, isStim = 1.0) 
        self.f_1 = self.stimuLus

        #self.stimuLus = Expression('x[2] < -5.6 ? 3.0*isStim : 0.0', degree = 1, isStim = 1)
        # why not to interpolate or project
        # https://fenicsproject.org/qa/13372/update-expression-and-time-dependent-assembly 

        self.mesh = mesh
        self.dx = dx
        self.ds = ds

    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    
    
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def update_state(self, w_new):
        self.w_n_fhn = w_new 
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 


    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def set_BCs(self):
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
        # Boundary conditions

        W_fhn = self.W_fhn
        #mesh = self.mesh

        #FHN1_bot = Expression(("0.4*isBC"), degree = 1, isBC = 1.0)
        #FHN2_bot = Expression(("0.0*isBC"), degree = 1, isBC = 1.0)
        #bcBot_FHN1 = DirichletBC(W_fhn.sub(0), FHN1_bot, bcBottom)
        #bcBot_FHN2 = DirichletBC(W_fhn.sub(1), FHN2_bot, bcBottom)
        #bcs_FHN = [bcBot_FHN1, bcBot_FHN2] 
        bcs_FHN = [] 

        return bcs_FHN
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 



    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def set_ICs(self):

        W_fhn = self.W_fhn
        # Set initial conditions
        phi0 = interpolate(Expression('0.0', degree = 0), W_fhn.sub(0).collapse())
        r0 = interpolate(Expression('0.0', degree = 0), W_fhn.sub(1).collapse())
        Ta0_fhn = interpolate(Expression('0.0', degree = 0), W_fhn.sub(2).collapse())

        # Define test, trial, and required functions 
        w_n_fhn = Function(W_fhn) 
        assign(w_n_fhn, [phi0, r0, Ta0_fhn])

        return w_n_fhn 
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 



    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def update_state_w_n(self, w_n_fhn):

        self.w_n_fhn = w_n_fhn
        phi_n, r_n, Ta_n_fhn = split(self.w_n_fhn)

        self.phi_n = phi_n 
        self.r_n = r_n 
        self.Ta_n_fhn = Ta_n_fhn 

    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 



    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    def FHN_F_J(self, J_fhn=None, C_fhn=None): 
        
        phi = self.phi
        r = self.r
        Ta_fhn = self.Ta_fhn
        
        phi_n = self.phi_n
        r_n = self.r_n 
        Ta_n_fhn = self.Ta_n_fhn

        phi_test = self.phi_test
        r_test = self.r_test
        Ta_test_fhn = self.Ta_test_fhn
        mesh = self.mesh
        #print mesh.ufl_cell().geometric_dimension()


        if J_fhn == None:
            J_n = Constant(1.0)
        else:
            J_n = J_fhn

        if C_fhn == None:
            C_n = Identity(mesh.ufl_cell().geometric_dimension())
        else:
            C_n = C_fhn

        myC = Constant_definitions() 

        alpha = myC.parameters['alpha']
        g = myC.parameters['gamma']
        b = myC.parameters['b']
        c = myC.parameters['c']
        mu1 = myC.parameters['mu1']
        mu2 = myC.parameters['mu2']

        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
        def eps_FHN(phi, r):
            return (g + (mu1*r)/(mu2 + phi))

        def f_phi(phi, r):
            return -c*phi*(phi - alpha)*(phi - 1) - r*phi

        def f_r(phi, r):
            return eps_FHN(phi,r)*(- r -c*phi*(phi - b - 1.))
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 


        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 

        f_phi = f_phi(phi,r)
        f_r = f_r(phi,r)
        
        # Parameters for active stress 
        e0 = Constant(1.0)
        phi_mid = Constant(0.05)
        ephi = conditional(gt(phi, phi_mid), e0, 10*e0)
        kTa = Constant(2*84000.0)

        f_1 = self.f_1

        timeDt = self.parameters["time_step"]
        k = Constant(timeDt)

        Dmat = self.parameters["Diffusion_Tensor"]
        #isLBBB = self.parameters["LBBB"]
        dx = self.dx

        F_FHN = ( (phi - phi_n) / k)*phi_test*dx + 1/J_n*dot(J_n*inv(C_n)*Dmat*grad(phi), grad(phi_test))*dx \
            - f_phi*phi_test*dx \
            + ( (r - r_n) / k)*r_test*dx \
            - f_r*r_test*dx \
            + ( (Ta_fhn - Ta_n_fhn) / k)*Ta_test_fhn*dx \
            - ephi*( kTa*phi - Ta_fhn )*Ta_test_fhn*dx \
            - f_1*phi_test*dx 
            
        J_FHN = derivative(F_FHN, self.w_fhn, self.dw_fhn)
        
        return (F_FHN, J_FHN)
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 


    def write_state(self, t, tau):
        outputfolder = self.parameters["outputfolder"]
        foldername = self.parameters["foldername"]        


        vtkfile_phi = self.vtkfile_phi
        vtkfile_r = self.vtkfile_r
        vtkfile_Ta_fhn = self.vtkfile_Ta_fhn
        #vtkfile_q = self.vtkfile_q

        w_n_fhn = self.w_n_fhn
        phi_, r_, Ta_fhn_ = w_n_fhn.split(deepcopy=True)

        print('Time = : %0.2f, tau = : %0.2f' % (t, tau ) ) 
        print('Tau is arbitrary' ) 
        vtkfile_phi << phi_
        vtkfile_r << r_
        vtkfile_Ta_fhn << Ta_fhn_ 

        #https://fenicsproject.org/qa/10159/evaluating-integrals-on-domain
        qVec =  [assemble(grad(phi_)[i]*dx) for i in range(3)] 
        print qVec[:]
        qVecNorm = np.linalg.norm(qVec)
        print >> self.fdataECG, t, qVec[0], qVec[1], qVec[2], qVecNorm
        

    def advance_timeStepssolver_FHN(self, solver_FHN, fhn_steps): 

        #pdb.set_trace()
        for ii in range(0, fhn_steps):
            print 'Solving FHN'
            print 'Solving ...' 
            solver_FHN.solvenonlinear()
            print 'FHN equations solved'
            #fhn_obj.write_state(state_obj.t, tau)
            self.w_n_fhn.assign(self.w_fhn)
            self.update_state_w_n(self.w_n_fhn)

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -



#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class State_Variables(object):
    """
    State variables for LV and RV, ejecting, relaxing, and filling 
    """
    def __init__(self, mpi_comm, SimDet):

    	self.isLV = SimDet["isLV"]

        # Systole 
        self.isLVeject = 0
        self.isLVfill = 0
        self.isLVfilling = False
        self.isLVejecting = False
        self.isLVrelaxing = False

	if(not self.isLV):
        	self.isRVeject = 0 
        	self.isRVfill = 0
        	self.isRVfilling = False
        	self.isRVejecting = False
        	self.isRVrelaxing = False

        self.tstep = 0
        self.BCL = SimDet["HeartBeatLength"] 

        self.cycle = 0.0
        self.t = 0
        self.tstep = 0
        self.dt = Expression(("dt"), dt=0.0, degree=1)

        self.comm = mpi_comm

    def check_isnot_isovolumic(self):        
	if(self.isLV):
        	return bool(self.isLVeject) or bool(self.isLVfill)
	else:
        	return bool(self.isLVeject) or bool(self.isRVeject) or bool(self.isLVfill) or bool(self.isRVfill)


    def update_state(self, pv_o, circuit_o, t_a):

        if(self.isLVejecting and circuit_o.Qlv < 0.005):
            self.isLVeject = 0
            self.isLVejecting = False
            self.isLVrelaxing = True 

            pv_o.LVESV = pv_o.LVV_cav

        if(pv_o.LVP_cav*0.0075 > circuit_o.Pao and not(self.isLVrelaxing)):           
            self.isLVeject = 1 
            self.isLVejecting = True

	if(not self.isLV):
        	if(self.isRVejecting and circuit_o.Qrv < 0.005):
        	    self.isRVeject = 0
        	    self.isRVejecting = False
        	    self.isRVrelaxing = True

        	    pv_o.RVESV = pv_o.RVV_cav

        	if(pv_o.RVP_cav*0.0075 > circuit_o.Ppu and not(self.isRVrelaxing)):
        	    self.isRVeject = 1 
        	    self.isRVejecting = True

        if(pv_o.LVP_cav*0.0075 < circuit_o.Pmv and self.isLVrelaxing and pv_o.LVPfilling_step  == 0):
            self.isLVfill = 1
            self.isLVfilling = True
            
            #self.dt.dt = 10.0#1.0
            #LVVfilling_step = (LVEDV - LVESV)/(BCL - t_a.t_a)#*dt.dt
            pv_o.LVPfilling_step = (pv_o.LVEDP - pv_o.LVP_cav)/(self.BCL - t_a.t_a)*self.dt.dt
            pv_o.LVVfilling_step = (pv_o.LVEDV - pv_o.LVESV)/(self.BCL - t_a.t_a)*self.dt.dt

	if(not self.isLV):
        	if(pv_o.RVP_cav*0.0075 < circuit_o.Ptr and self.isRVrelaxing and pv_o.RVPfilling_step ==  0):
        	    self.isRVfill = 1
        	    self.isRVfilling = True
        	    
        	    #self.dt.dt = 10.0#1.0
        	    #RVVfilling_step = (RVEDV - RVESV)/(BCL - t_a.t_a)#*dt.dt
        	    pv_o.RVPfilling_step = (pv_o.RVEDP - pv_o.RVP_cav)/(self.BCL - t_a.t_a)*self.dt.dt
        	    pv_o.RVVfilling_step = (pv_o.RVEDV - pv_o.RVESV)/(self.BCL - t_a.t_a)*self.dt.dt   


    def print_time_details(self):
        #printout("Cycle number = " + str(cycle) + " cell time = " + str(t) +  " tstep = " + str(tstep) + " dt = " + str(dt.dt), comm) 
        string_to_print = "Cycle number = " + str(self.cycle) \
                          + " cell time = " + str(self.t) \
                          + " tstep = " + str(self.tstep) \
                          + " dt = " + str(self.dt.dt) 
        printout(string_to_print, self.comm)

    def print_LV_state(self):

        string_to_print = "isLVejecting: " + str(self.isLVejecting) \
                          + " isLVrelaxing: " + str(self.isLVrelaxing) \
                          + " isLVfilling: " + str(self.isLVfilling)
        printout(string_to_print, self.comm)
 
    def print_RV_state(self):

        string_to_print = "isRVejecting: " + str(self.isRVejecting) \
                          + " isRVrelaxing: " + str(self.isRVrelaxing) \
                          + " isRVfilling: " + str(self.isRVfilling)
        printout(string_to_print, self.comm)

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -



#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class Windkessel(object):
    def __init__(self, elas_o, mpi_comm, SimDet):
        self.Pao = SimDet['Pao'] #50 # LV ejection start
    	self.isLV = SimDet["isLV"]
	if(not self.isLV):
        	self.Ppu = SimDet['Ppu'] #12 # RV ejection start 

        self.Qlv = 0.0
        self.Qlv_prev = 0.0

        self.Qrv = 0.0
        self.Qrv_prev = 0.0

        self.Cmpl_Ao = 0.007*1.2  #0.01 
        self.Cmpl_Pa = 0.01 #0.005 

        self.Rper_sys = 100000*1.2 #100000 MRC latest results  
        self.Rper_pul = 60000 

        self.Rao = 5500*1.2 #800
        self.Rpa = 500

        self.Pper_sys = SimDet['Pper_sys']#600 #(Pa)
        self.Pper_pul = SimDet['Pper_pul']#200 # RV pressure 
        #200 pulmonary terminal pressure too high? RV pressure PAH #(Pa)

        self.Pmv = SimDet['Pmv'] #5 #mmHg
        self.Ptr = SimDet['Ptr'] #4 #mmHg

	if(self.isLV):
        	self.F = 0
        	self.Vdiff = 0 
        	self.dFdV = 0
        	self.Ch = 0
        else:
        	self.F = np.zeros((2,1))
        	self.Vdiff = np.zeros((2,1))
        	self.dFdV = np.array([[0,0],[0,0]])
        	self.Ch = np.array([[0,0],[0,0]]) 


        W_elas = elas_o.W
        self.w_cur = Function(W_elas)

        self.normVdiff = 1e8
        self.normFdiff = 1e8

        self.comm = mpi_comm

    def calculate_F_Q(self, state_o, pv_o):
        
        Pper_sys = self.Pper_sys
        Pper_pul = self.Pper_pul

        Rper_sys = self.Rper_sys
        Rper_pul = self.Rper_pul

        Cmpl_Ao = self.Cmpl_Ao
        Cmpl_Pa = self.Cmpl_Pa

        Qlv = self.Qlv
        Qlv_prev = self.Qlv_prev 

	if(not self.isLV):
        	Qrv = self.Qrv
        	Qrv_prev = self.Qrv_prev 

        Rao = self.Rao
        Rpa = self.Rpa

        isLVeject = state_o.isLVeject
        isLVfill = state_o.isLVfill 
	if(not self.isLV):
        	isRVeject = state_o.isRVeject
        	isRVfill = state_o.isRVfill


        LVP_cav = pv_o.LVP_cav
        LVV_cav = pv_o.LVV_cav

        LVP_cav_prev = pv_o.LVP_cav_prev
        LVV_cav_prev = pv_o.LVV_cav_prev

	if(not self.isLV):
        	RVP_cav = pv_o.RVP_cav
        	RVV_cav = pv_o.RVV_cav

        	RVP_cav_prev = pv_o.RVP_cav_prev
        	RVV_cav_prev = pv_o.RVV_cav_prev

        LVVfilling_step = pv_o.LVVfilling_step
        LVPfilling_step = pv_o.LVPfilling_step
	if(not self.isLV):
        	RVVfilling_step = pv_o.RVVfilling_step
        	RVPfilling_step = pv_o.RVPfilling_step 

        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        if(bool(isLVeject)):
            Qlv = ((LVP_cav - Pper_sys)/(Rper_sys) + Cmpl_Ao*(LVP_cav - LVP_cav_prev + Rao*Qlv_prev)/state_o.dt.dt)/(1.0 + Rao/Rper_sys + Cmpl_Ao*Rao/state_o.dt.dt)
	    if(self.isLV):
            	self.F = LVV_cav - LVV_cav_prev + state_o.dt.dt*Qlv
	    else:
            	self.F[0] = LVV_cav - LVV_cav_prev + state_o.dt.dt*Qlv
        elif(bool(isLVfill)):

            printout("LVPfilling_step = " + str(LVPfilling_step), self.comm)
	    if(self.isLV):
            	self.F = LVP_cav - LVP_cav_prev - LVPfilling_step
	    else:
            	self.F[0] = LVP_cav - LVP_cav_prev - LVPfilling_step

        else:	
            Qlv = 0.0
	    if(self.isLV):
            	self.F = LVV_cav - LVV_cav_prev
	    else:
            	self.F[0] = LVV_cav - LVV_cav_prev
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
	if(not self.isLV):
        	if(bool(isRVeject)):
        	    Qrv = ((RVP_cav - Pper_pul)/(Rper_pul) + Cmpl_Pa*(RVP_cav - RVP_cav_prev + Rpa*Qrv_prev)/state_o.dt.dt)/(1.0 + Rpa/Rper_pul + Cmpl_Pa*Rpa/state_o.dt.dt)
        	    self.F[1] = RVV_cav - RVV_cav_prev + state_o.dt.dt*Qrv
        	elif(bool(isRVfill)):

        	    printout("RVPfilling_step = " + str(RVPfilling_step), self.comm)
        	    self.F[1] = RVP_cav - RVP_cav_prev - RVPfilling_step
			
        	else:
        	    Qrv = 0.0
        	    self.F[1] = RVV_cav - RVV_cav_prev
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        
	if(self.isLV):
        	printout("Qlv = " + str(Qlv), self.comm)
	else:
        	printout("Qlv = " + str(Qlv) +  " Qrv = " + str(Qrv), self.comm)

        self.normFdiff = np.linalg.norm(self.F)


	if(self.isLV):
        	string_to_print = "Fdiff = " + str(self.normFdiff) \
                          + " F[0] = " + str(self.F) 
	else:
        	string_to_print = "Fdiff = " + str(self.normFdiff) \
                          + " F[0] = " + str(self.F[0]) \
                          + " F[1] = " + str(self.F[1])

        printout(string_to_print, self.comm)

        self.Qlv = Qlv
	if(not self.isLV):
        	self.Qrv = Qrv

        self.AorPres = LVP_cav - LVP_cav_prev + Rao*Qlv_prev
	if(not self.isLV):
        	self.PulPres = RVP_cav - RVP_cav_prev + Rpa*Qrv_prev

    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
    def set_dFdV(self, state_o, pv_o):

        Pper_sys = self.Pper_sys
        Pper_pul = self.Pper_pul

        Rper_sys = self.Rper_sys
        Rper_pul = self.Rper_pul

        Cmpl_Ao = self.Cmpl_Ao
        Cmpl_Pa = self.Cmpl_Pa

        Qlv = self.Qlv
        Qlv_prev = self.Qlv_prev 

        Qrv = self.Qrv
        Qrv_prev = self.Qrv_prev 

        Rao = self.Rao
        Rpa = self.Rpa

        isLVeject = state_o.isLVeject
        isLVfill = state_o.isLVfill 
	if(not self.isLV):
        	isRVeject = state_o.isRVeject
        	isRVfill = state_o.isRVfill 

        Ch = self.Ch


	if(self.isLV):
        	dFdV  = 1.0 + float(isLVeject)*Ch*(state_o.dt.dt/Rper_sys + Cmpl_Ao)/(1.0 + Rao/Rper_sys + Cmpl_Ao*Rao/state_o.dt.dt) + float(isLVfill)*Ch

        else:
	        dFdV11 = 1.0 + float(isLVeject)*Ch[0,0]*(state_o.dt.dt/Rper_sys + Cmpl_Ao)/(1.0 + Rao/Rper_sys + Cmpl_Ao*Rao/state_o.dt.dt) + float(isLVfill)*Ch[0,0]
	        dFdV12 = float(isLVeject)*Ch[0,1]*(state_o.dt.dt/Rper_sys + Cmpl_Ao)/(1.0 + Rao/Rper_sys + Cmpl_Ao*Rao/state_o.dt.dt)  + float(isLVfill)*Ch[0,1]
	        dFdV21 = float(isRVeject)*Ch[1,0]*(state_o.dt.dt/Rper_pul + Cmpl_Pa)/(1.0 + Rpa/Rper_pul + Cmpl_Pa*Rpa/state_o.dt.dt)  + float(isRVfill)*Ch[1,0]
	        dFdV22 = 1.0 + float(isRVeject)*Ch[1,1]*(state_o.dt.dt/Rper_pul + Cmpl_Pa)/(1.0 + Rpa/Rper_pul + Cmpl_Pa*Rpa/state_o.dt.dt)  + float(isRVfill)*Ch[1,1]
	        
        	dFdV = np.matrix([[dFdV11, dFdV12], [dFdV21, dFdV22]])
        self.dFdV = dFdV
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

    def update_Q(self):
        self.Qlv_prev = self.Qlv 
	if(not self.isLV):
        	self.Qrv_prev = self.Qrv 


    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
    def computeCompliance(self, solver_elas, pv_o, LVCavityvol, RVCavityvol, uflforms):
    
        LVP_cav = pv_o.LVP_cav
        LVV_cav = pv_o.LVV_cav

        LVP_cav_prev = pv_o.LVP_cav_prev
        LVV_cav_prev = pv_o.LVV_cav_prev

	if(not self.isLV):
        	RVP_cav = pv_o.RVP_cav
        	RVV_cav = pv_o.RVV_cav

        w_cur = self.w_cur
        w = solver_elas.parameters["w"]

        printout("######################################################", self.comm)
        printout("Calculating compliance ", self.comm)
        w_cur.vector()[:] = w.vector().array()[:]
        dV = 0.01
        LVCavityvol.vol = LVV_cav + dV
	if(not self.isLV):
        	RVCavityvol.vol = RVV_cav
        solver_elas.solvenonlinear()

	if(self.isLV):
        	dPlv_dVlv = (uflforms.LVcavitypressure() - LVP_cav)/dV
        	printout("dPlv_dVlv = " + str(dPlv_dVlv), self.comm)
	else:
        	dPrv_dVlv = (uflforms.RVcavitypressure() - RVP_cav)/dV
        	printout("dPlv_dVlv = " + str(dPlv_dVlv) +  " dPrv_dVlv = " + str(dPrv_dVlv), self.comm)

        w.vector()[:] = w_cur.vector().array()[:]
        LVCavityvol.vol = LVV_cav
	if(not self.isLV):
        	RVCavityvol.vol = RVV_cav + dV
        	solver_elas.solvenonlinear()

        	dPlv_dVrv = (uflforms.LVcavitypressure() - LVP_cav)/dV
        	dPrv_dVrv = (uflforms.RVcavitypressure() - RVP_cav)/dV
        	w.vector()[:] = w_cur.vector().array()[:]
        	printout("dPlv_dVrv = " + str(dPlv_dVrv) +  " dPrv_dVrv = " + str(dPrv_dVrv), self.comm)
        	printout("######################################################", self.comm)

        	RVCavityvol.vol = RVV_cav

	if(self.isLV):
        	self.Ch = dPlv_dVlv
	else:
        	self.Ch = np.array([[dPlv_dVlv, dPlv_dVrv], [dPrv_dVlv, dPrv_dVrv]])

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -



#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
class PV_Ventricles(object):
    def __init__(self, uflforms, mpi_comm, SimDet):

    	self.isLV = SimDet["isLV"]

        LVV_unload =  uflforms.LVcavityvol()
	if(not self.isLV):
        	RVV_unload =  uflforms.RVcavityvol()

        prescribedLVEDV_Shift = SimDet["LVEDV_Shift"]
	if(not self.isLV):
        	prescribedRVEDV_Shift = SimDet["RVEDV_Shift"]

        LVEDV = LVV_unload + prescribedLVEDV_Shift #30   
	if(not self.isLV):
        	RVEDV = RVV_unload + prescribedRVEDV_Shift #38  

        #LVEDV = LVV_unload + 30 
        #RVEDV = RVV_unload + 38 
        
        #LVEDV = 30 
        #RVEDV = 38 
	
        LVESV = 0
        RVESV = 0

        self.LVV_unload = LVV_unload
	if(not self.isLV):
        	self.RVV_unload = RVV_unload
  
        self.LVEDV = LVEDV
        self.LVESV = LVESV
	if(not self.isLV):
        	self.RVEDV = RVEDV
        	self.RVESV = RVESV

        LVEDP = 0.
        LVESP = 0.
	if(not self.isLV):
        	RVEDP = 0. 
        	RVESP = 0.

        self.LVEDP = LVEDP
        self.LVESP = LVESP
	if(not self.isLV):
        	self.RVEDP = RVEDP
        	self.RVESP = RVESP

        LVP_cav = 0.
        LVV_cav = 0.
	if(not self.isLV):
        	RVP_cav = 0.
        	RVV_cav = 0.

        LVP_cav_prev = 0.
        LVV_cav_prev = 0.
	if(not self.isLV):
        	RVP_cav_prev = 0.
        	RVV_cav_prev = 0.

        self.LVP_cav = LVP_cav
        self.LVV_cav = LVV_cav
	if(not self.isLV):
        	self.RVP_cav = RVP_cav
        	self.RVV_cav = RVV_cav

        self.LVP_cav_prev = LVP_cav_prev
        self.LVV_cav_prev = LVV_cav_prev
	if(not self.isLV):
        	self.RVP_cav_prev = RVP_cav_prev
        	self.RVV_cav_prev = RVV_cav_prev 

        LVVfilling_step = 0
        LVPfilling_step = 0
	if(not self.isLV):
        	RVVfilling_step = 0
        	RVPfilling_step = 0

        self.LVVfilling_step = LVVfilling_step
        self.LVPfilling_step = LVPfilling_step
	if(not self.isLV):
        	self.RVVfilling_step = RVVfilling_step
        	self.RVPfilling_step = RVPfilling_step 

        self.comm = mpi_comm 

    def update_pv(self, uflforms):
    	self.LVP_cav = uflforms.LVcavitypressure()
        self.LVV_cav = uflforms.LVcavityvol()
	if(not self.isLV):
    		self.RVP_cav = uflforms.RVcavitypressure()
        	self.RVV_cav = uflforms.RVcavityvol()


    def update_pv_prev(self):
    	self.LVP_cav_prev = self.LVP_cav
        self.LVV_cav_prev = self.LVV_cav
	if(not self.isLV):
        	self.RVP_cav_prev = self.RVP_cav
        	self.RVV_cav_prev = self.RVV_cav

    def print_PV_details(self):
	if(not self.isLV):
        	printout("LVP_cav = " + str(self.LVP_cav*0.0075) \
        	         + " LVV_cav = " + str(self.LVV_cav) \
        	         + " RVP_cav = " + str(self.RVP_cav*0.0075) \
        	         + " RVV_cav = " + str(self.RVV_cav), \
        	         self.comm)
	else:
		printout("LVP_cav = " + str(self.LVP_cav*0.0075) \
        	         + " LVV_cav = " + str(self.LVV_cav), \
        	         self.comm)

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -


class exportfiles(object): 

	def __init__(self, mpi_comm_me, mpi_comm_ep, IODet, SimDet):

    		self.isLV = SimDet["isLV"]
		self.outputfolder = IODet["outputfolder"]
		self.folderName = IODet["folderName"] + IODet["caseID"] + '/' 
        	self.comm_me = mpi_comm_me
        	self.comm_ep = mpi_comm_ep

		self.removeAllfiles()
		self.opendatafilestreams()
		self.openVTKfilestreams()
		self.openHDF5filestreams()

	def removeAllfiles(self):

		outputfolder = self.outputfolder
		folderName = self.folderName

    		os.system("rm "+ outputfolder + folderName + "*.pvd")
    		os.system("rm "+ outputfolder + folderName + "*.vtu")
    		os.system("rm "+ outputfolder + folderName + "*.pvtu")
    		os.system("rm "+ outputfolder + folderName + "*.txt")
    		os.system("rm "+ outputfolder + folderName + "*.hdf5")
    		os.system("rm -rf "+ outputfolder + folderName + "FHN")
    		os.system("rm -rf "+ outputfolder + folderName + "IMP")
    		os.system("rm -rf "+ outputfolder + folderName + "IMP2")
    		os.system("rm -rf "+ outputfolder + folderName + "IMP_J")
    		os.system("rm -rf "+ outputfolder + folderName + "IMP_Constraint")
    		os.system("rm -rf "+ outputfolder + folderName + "FHN_coarse")
    		os.system("rm -rf "+ outputfolder + folderName + "FHN_coarse2")
    		os.system("rm -rf "+ outputfolder + folderName + "FHN_coarseRef")
    		os.system("rm -rf "+ outputfolder + folderName + "deformation")
    		os.system("rm -rf "+ outputfolder + folderName + "deformation_loadED")
    		os.system("rm -rf "+ outputfolder + folderName + "deformation_unloadED")
    		os.system("rm -rf "+ outputfolder + folderName + "stretch")
    		os.system("rm -rf "+ outputfolder + folderName + "active")

	def exportVTKobj(self, filename, vtkobj):

		outputfolder = self.outputfolder
		folderName = self.folderName

    		File(outputfolder + folderName + filename) << vtkobj

		return

	def opendatafilestreams(self):

		outputfolder = self.outputfolder
		folderName = self.folderName
		comm = self.comm_ep


	        if(MPI.rank(comm) == 0):
		    self.fdatalog = open(outputfolder + folderName + "BiV_log.txt", "w", 0)
	            self.fdataPV = open(outputfolder + folderName + "BiV_PV.txt", "w", 0)
	            self.fdataQ = open(outputfolder + folderName + "BiV_Q.txt", "w", 0)
	            self.fdataP = open(outputfolder + folderName + "BiV_P.txt", "w", 0)
	            self.fdatatpt = open(outputfolder + folderName + "tpt.txt", "w", 0)
	   
	            self.fdataIMP = open(outputfolder + folderName + "BiV_IMP.txt", "w", 0)
	            self.fdataIMP2 = open(outputfolder + folderName + "BiV_IMP2.txt", "w", 0)
	            self.fdataIMP3 = open(outputfolder + folderName + "BiV_IMP_InC.txt", "w", 0)
	   
	            self.fdataStress = open(outputfolder + folderName + "BiV_fiberStress.txt", "w", 0)
	            self.fdataStrain = open(outputfolder + folderName + "BiV_fiberStrain.txt", "w", 0)
		    self.fdataStretch = open(outputfolder + folderName + "BiV_fiberStretch.txt", "w", 0)
	            self.fdataWork = open(outputfolder + folderName + "BiV_fiberWork.txt", "w", 0)
	   
	            self.fdata_C_Strain = open(outputfolder + folderName + "BiV_CStrain.txt", "w", 0)
	            self.fdata_L_Strain = open(outputfolder + folderName + "BiV_LStrain.txt", "w", 0)
	            self.fdata_R_Strain = open(outputfolder + folderName + "BiV_RStrain.txt", "w", 0)

		return


	def openVTKfilestreams(self):

		outputfolder = self.outputfolder
		folderName = self.folderName

	        self.vtkfile_phi  = File(outputfolder + folderName +  'FHN/phi.pvd')
    		self.vtkfile_phi_ref = File(outputfolder + folderName +  'FHN_coarseRef/phi_ref.pvd')
    		self.vtkfile_phi_me  = File(outputfolder + folderName +  'FHN_coarse/phi_me.pvd')
    		self.vtkfile_r  = File(outputfolder + folderName + 'FHN/r.pvd')

    		self.vtkfile_IMP  = File(outputfolder + folderName + 'IMP/IMP.pvd')
    		self.vtkfile_IMP2  = File(outputfolder + folderName + 'IMP2/IMP2.pvd')
    		self.vtkfile_IMP_Constraint  = File(outputfolder + folderName + 'IMP_Constraint/IMP_Constraint.pvd')

    		self.displacementfile = File(outputfolder + folderName +"deformation/u_disp.pvd")
    		self.displacementED_file = File(outputfolder + folderName +"/deformation_loadED/u_disp.pvd")

    		self.activationFile = File(outputfolder + folderName +"/active/activationTime.pvd")
    		self.vtkfile_Ecc = File(outputfolder + folderName +"/active/Ecc.pvd") 
    		self.vtkfile_Ell = File(outputfolder + folderName +"/active/Ell.pvd") 
    		self.vtkfile_Err = File(outputfolder + folderName +"/active/Err.pvd") 
    		self.vtkfile_Efiber = File(outputfolder + folderName +"/active/Eff.pvd") 
    		self.vtkfile_fstress = File(outputfolder + folderName +"/active/fstress.pvd") 
		
		return

	def openHDF5filestreams(self):

		outputfolder = self.outputfolder
		folderName = self.folderName
		comm_me = self.comm_me
		comm_ep = self.comm_ep

		# Set up hdf5
		self.hdf = HDF5File(comm_me, outputfolder + folderName + 'Data.h5', "w")

		return

	def closeHDF5filestreams(self):

		self.hdf.close();

	def writePV(self, MEmodel, t):

    		isLV = self.isLV
		comm = self.comm_ep

		LVP = MEmodel.GetLVP()*0.0075
		LVV = MEmodel.GetLVV()

		if(not isLV):
			RVP = MEmodel.GetRVP()*0.0075
			RVV = MEmodel.GetRVV()

		if(MPI.rank(comm) == 0):
			fdataPV = self.fdataPV
			if(not isLV):
            			print >> fdataPV, t, LVP, LVV, \
                     		 	RVP, RVV
			else:
				print >> fdataPV, t, LVP, LVV
 
		return

	def writeQ(self, MEmodel, Qarray, t):

    		isLV = self.isLV
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataQ = self.fdataQ
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdataQ, t, ' '.join(map(str, Qarray))
 
		return



	def writeP(self, MEmodel, Parray, t):

    		isLV = self.isLV
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataP = self.fdataP
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdataP, t, ' '.join(map(str, Parray))
 
		return

	def writetpt(self, MEmodel, tpt):

		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdatatpt = self.fdatatpt
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
			print >> fdatatpt, tpt
 
		return

	def writeIMP(self, MEmodel, t, fIMP):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataIMP = self.fdataIMP
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdataIMP, t, ' '.join(map(str, fIMP))

		return

	def writeIMP2(self, MEmodel, t, fIMP2):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataIMP2 = self.fdataIMP2
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdataIMP2, t, ' '.join(map(str, fIMP2))

		return

	def writeIMP3(self, MEmodel, t, fIMP3):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataIMP3 = self.fdataIMP3
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdataIMP3, t, ' '.join(map(str, fIMP3))

		return



	def writefStress(self, MEmodel, t, fStress):
		
		comm = self.comm_ep

        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
		if(MPI.rank(comm) == 0):
			fdataStress = self.fdataStress
            		#print >> fdataStress, t, ' '.join(map(str, fStress))
			print >> fdataStress, t,  fStress

		return

	def writefStrain(self, MEmodel, t, fStrain):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataStrain = self.fdataStrain
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdataStrain, t, ' '.join(map(str, fStrain))

		return

	def writefStretch(self, MEmodel, t, fStretch):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdataStretch = self.fdataStretch
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		#print >> fdataStretch, t, ' '.join(map(str, fStretch))
			print >> fdataStretch, t,  fStretch

		return


	def writeCStrain(self, MEmodel, t, CStrain):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdata_C_Strain = self.fdata_C_Strain
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdata_C_Strain, t, ' '.join(map(str, CStrain))

		return

	def writeLStrain(self, MEmodel, t, LStrain):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdata_C_Strain = self.fdata_C_Strain
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdata_C_Strain, t, ' '.join(map(str, LStrain))

		return

	def writeRStrain(self, MEmodel, t, RStrain):
		
		comm = self.comm_ep

		if(MPI.rank(comm) == 0):
			fdata_C_Strain = self.fdata_C_Strain
        	#if(MPI.rank(MEmodel.mesh_me.mpi_comm()) == 0):
            		print >> fdata_C_Strain, t, ' '.join(map(str, RStrain))

		return

	def printout(self, statement): 
		comm = self.comm_ep
    		if(MPI.rank(comm) == 0): 
			fdata_log = self.fdatalog
        		print >> fdata_log, statement



