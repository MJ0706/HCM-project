from dolfin import *
import math as math
import numpy as np
from ..mechanics.forms_MRC2 import Forms 
from ..mechanics.activeforms_MRC2 import activeForms 
from ..ep.EPmodel import EPmodel
from ..utils.oops_objects_MRC2 import State_Variables
from ..utils.nsolver import NSolver as NSolver

def run_isotonic_EM(IODet, SimDet):

    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]['representation'] = 'uflacs'
    parameters["form_compiler"]["quadrature_degree"] = 4

    outputfolder = IODet["outputfolder"]
    folderName = IODet["folderName"] + IODet["caseID"] + '/' 

    SimDet_ = {
	     "length" : 1,
    	     "width" : 1,
    	     "nelem_ep" : 1,
    	     "nelem_me" : 1,
    	     "nload" : 50,
    	     "applied_load" : 1.0,
    	     "dt" : 1.0,
    	     "ntpt" : 800
	     "HomogenousActivation": True,
	    };
    SimDet_.update(SimDet)

    # Create Mesh
    meshEP = createmesh(IODet, SimDet_, SimDet_["nelem_ep"])
    File(outputfolder + folderName + "subdomain_ep.pvd") <<  meshEP.subdomains
    File(outputfolder + folderName + "mesh_ep.pvd") << meshEP.mesh

    meshME = createmesh(IODet, SimDet_, SimDet_["nelem_me"])
    File(outputfolder + folderName + "subdomain_me.pvd") <<  meshME.subdomains
    File(outputfolder + folderName + "mesh_me.pvd") << meshME.mesh

    # Set Fiber
    f0 = Expression(("1.0","0.0","0.0"), degree=1)
    s0 = Expression(("0.0","1.0","0.0"), degree=1)
    n0 = Expression(("0.0","0.0","1.0"), degree=1)

    # Define state variables
    state_obj = State_Variables(meshEP.comm, SimDet) 
    state_obj.dt.dt = SimDet["dt"]

    AHAid_ep = CellFunction('size_t', meshEP.mesh)
    AHAid_ep.set_all(0)

    matid_ep = CellFunction('size_t', meshEP.mesh)
    matid_ep.set_all(0)

    # Define EP model and solver
    EPparams = {"EPmesh": meshEP.mesh,\
                "deg": 4,\
        	"facetboundaries": meshEP.subdomains,\
        	"f0": f0,\
        	"s0": s0,\
        	"n0": n0,\
        	"state_obj": state_obj,\
                "d_iso": SimDet["d_iso"],\
                "d_ani_factor": SimDet["d_ani_factor"],\
                "AHAid": AHAid_ep,\
                "matid": matid_ep,
		"Ischemia": SimDet["Ischemia"],
		"pacing_timing": SimDet["pacing_timing"],
		"ploc": SimDet["ploc"]
		}

    EPmodel_ = EPmodel(EPparams)
    EpiBCid_ep = EPmodel_.MarkStimulus()
    solver_FHN = EPmodel_.Solver()
    File(outputfolder + folderName + "EP_stimulus.pvd") << EpiBCid_ep

    MEparams = {"mesh": meshME.mesh,\
		"facetboundaries": meshME.subdomains, \
		"facetnormal": meshME.N,\
        	"f0": f0,\
        	"s0": s0,\
        	"n0": n0,\
		}

    MEmodel_ = mechanics(IODet, SimDet_, MEparams)
    solver_ME = MEmodel_.Solver()
	
    pload_arr = [0.0]

    F_EP = File(outputfolder + folderName + "EP.pvd")
    F_Act = File(outputfolder + folderName + "Act.pvd")
    F_ActF = File(outputfolder + folderName + "ActF.pvd")
    F_Disp = File(outputfolder + folderName + "Disp.pvd")
    ntpt = SimDet["ntpt"]
    dt = SimDet["dt"]

    cnt = 0
    potential_me = Function(FunctionSpace(MEmodel_.mesh,'CG',1))
    for tpt in np.arange(0, ntpt):

        state_obj.tstep = state_obj.tstep + state_obj.dt.dt
        state_obj.cycle = math.floor(state_obj.tstep/state_obj.BCL)
        state_obj.t = state_obj.tstep - state_obj.cycle*state_obj.BCL

        MEmodel_.t_a.vector()[:] = state_obj.t
    	print("Active contraction =", MEmodel_.t_a.vector().array()[0], " State obj t = ", state_obj.t)

	# Reset phi and r in EP at end of diastole
	if state_obj.t < state_obj.dt.dt:
		EPmodel_.reset();

	# Solve EP
	print("Solve EP")
        solver_FHN.solvenonlinear()
	EPmodel_.UpdateVar()

	# Solve Mechanics
	print("Solve ME")
	solver_ME.solvenonlinear()
	MEmodel_.UpdateVar()

        potential_ref = EPmodel_.interpolate_potential_ep2me_phi(V_me = Function(FunctionSpace(MEmodel_.mesh,'CG',1)))
        potential_ref.rename("v_ref", "v_ref")
        potential_me.vector()[:] = potential_ref.vector().array()[:]

        MEmodel_.activeforms.update_activationTime(potential_n = potential_me, comm = meshME.mesh.mpi_comm())

        if(cnt % SimDet["writeStep"] == 0.0):
		F_EP << EPmodel_.getphivar()	
		F_Act << potential_me
		F_Disp << MEmodel_.GetDisplacement()
		F_ActF << MEmodel_.GetSActive()
   
	cnt += 1	
    #tpt = np.arange(0, len(load_array))*dt

    return #tpt, active_load_array, load_array, lbda_arr, activeforms.matparams, uflforms.matparams

class createmesh(object):

    def __init__(self, IODet, SimDet, nelem):

        self.IODet = IODet 
	self.SimDet = SimDet
    
	length = SimDet["length"]
	width = SimDet["width"]
        outputfolder = IODet["outputfolder"]
        folderName = IODet["folderName"] + IODet["caseID"] + '/' 

    	# Create mesh 
    	mesh = BoxMesh(Point(0,0,0), Point(length,width,width), nelem*int(length/width), nelem, nelem)
    	print("Number of Elements in mesh = "+str(mesh.num_cells()))
	self.mesh = mesh

	self.N = FacetNormal(mesh)
	
	# Mark Facet
	class Right(SubDomain):
	    def inside(self, x, on_boundary):
	        return x[0] > length*(1.0 - DOLFIN_EPS)and on_boundary
	
	class Left(SubDomain):
	    def inside(self, x, on_boundary):
	        return x[0] < DOLFIN_EPS and on_boundary
	
	class Front(SubDomain):
	    def inside(self, x, on_boundary):
	        return x[1] > width*(1.0 - DOLFIN_EPS) and on_boundary
	
	class Back(SubDomain):
	    def inside(self, x, on_boundary):
	        return x[1] < DOLFIN_EPS and on_boundary
	
	class Top(SubDomain):
	    def inside(self, x, on_boundary):
	        return x[2] > width*(1.0 - DOLFIN_EPS) and on_boundary
	
	class Bot(SubDomain):
	    def inside(self, x, on_boundary):
	        return x[2] < DOLFIN_EPS and on_boundary
	
	# Create mesh functions over the cell facets
	sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, mesh.domains())
	
	# Mark all facets as sub domain 1
	sub_domains.set_all(0)
	
	# Mark left facets
	left = Left()
	left.mark(sub_domains, 1)
	
	right = Right()
	right.mark(sub_domains, 2)       
	
	front = Front()
	front.mark(sub_domains, 3)
	
	back = Back()
	back.mark(sub_domains, 4)       
	
	top = Top()
	top.mark(sub_domains, 5)
	
	bot = Bot()
	bot.mark(sub_domains, 6)
	
	self.subdomains = sub_domains	
	self.comm = mesh.mpi_comm()
	

class mechanics(object):

    def __init__(self, IODet, SimDet, MEparams):

	passive_mat_input = {}
	active_mat_input = {}
	if("Passive model" in SimDet.keys()):
	    	matmodel = SimDet["Passive model"]
	    	passive_mat_input.update({"material model": matmodel})
	if("Passive params" in SimDet.keys()):
	    	matparams= SimDet["Passive params"]
	    	passive_mat_input.update({"material params": matparams})
	if("Active model" in SimDet.keys()):
	    	matmodel = SimDet["Active model"]
	    	active_mat_input.update({"material model": matmodel})
	if("Active params" in SimDet.keys()):
	    	matparams= SimDet["Active params"]
	    	active_mat_input.update({"material params": matparams})

	mesh = MEparams["mesh"]
	sub_domains = MEparams["facetboundaries"]
	f0 = MEparams["f0"]
	s0 = MEparams["s0"]
	n0 = MEparams["n0"]
	N = MEparams["facetnormal"]
	#t_a = MEparams["t_a"]

        Quadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=4, quad_scheme="default")
        Quadelem._quad_scheme = 'default'
	t_a = Function(FunctionSpace(mesh, Quadelem))

	# Define Integration domain
	dx =  dolfin.dx(mesh, metadata = {"integration_order": 2})
	ds =  dolfin.ds(mesh, subdomain_data = sub_domains)
	
	# Define function space
	Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
	Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
	DGelem = FiniteElement("DG", mesh.ufl_cell(), 1, quad_scheme="default")
	
	W = FunctionSpace(mesh, MixedElement([Velem, Qelem]))
	V = FunctionSpace(mesh, Velem)
	QDG =  FunctionSpace(mesh, DGelem)
	
	# Define function
	w = Function(W)
	w_n = Function(W)
	wtest = TestFunction(W)
	dw = TrialFunction(W)
	du, dp = TrialFunctions(W)
	(u, p) = split(w)
	(v, q) = split(wtest)
	
	params= {"mesh": mesh,
	     "facetboundaries": sub_domains,
	     "facet_normal": N,
	     "mixedfunctionspace": W,
	     "mixedfunction": w,
	     "displacement_variable": u, 
	     "pressure_variable": p,
	     "fiber": f0,
	     "sheet": s0,
	     "sheet-normal": n0,
	     "growth_tensor": None,
	     "incompressible" : True,
	    }
	
	params.update(passive_mat_input)
	uflforms = Forms(params)
	
	SEF = uflforms.PassiveMatSEF()
	PK1pas = uflforms.PK1()
	Fe = uflforms.Fmat()
	n = uflforms.J()*inv(Fe.T)*N    
	
	#_a = Expression(("ta"), ta = 0.0, degree=0)
	activeparams = {"mesh": mesh,
	                "dx":dx,
	   	    	"deg":4,
	                "facetboundaries": sub_domains,
	                "facet_normal": N,
	                "displacement_variable": u, 
	                "pressure_variable": p,
	                "fiber": f0,
	                "sheet": s0,
	                "sheet-normal": n0,
	                "t_a": t_a,
	                "mesh": mesh, 
	                "Threshold_Potential": 0.9,
	                "growth_tensor": None,
	   	        "HomogenousActivation": SimDet_["HomogenousActivation"]
	           }
	activeparams.update(active_mat_input)
	
	activeforms = activeForms(activeparams)
	potential = Function(FunctionSpace(mesh,'CG',1))
	activeforms.update_activationTime(potential_n = potential, comm = mesh.mpi_comm())
	
	# Define Passive SEF
	SEF = uflforms.PassiveMatSEF()
	PK1pas = uflforms.PK1()
	Fe = uflforms.Fmat()
	n = uflforms.J()*inv(Fe.T)*N    
	
	# Define Active Stress
	Sactive = activeforms.PK2StressTensor()
	PK1act = Fe*Sactive
	
	# Boundary conditions
	uright = Expression("val", val = 0, degree=0)
	bcleft = DirichletBC(W.sub(0).sub(0), Constant((0.0)), sub_domains, 1)
	bcright = DirichletBC(W.sub(0).sub(0), uright, sub_domains, 2)
	bcback = DirichletBC(W.sub(0).sub(1), Constant((0.0)), sub_domains, 4)
	bcbot = DirichletBC(W.sub(0).sub(2), Constant((0.0)), sub_domains, 6)
	bcs = [bcleft, bcback, bcbot]
	pload = Expression("val", val = 0.0, degree=0)
	
	# Define Weak Form
	Fact = inner(PK1act, grad(v))*dx
	Fpas = derivative(SEF, w, wtest)*dx
	Fload =  -inner(pload*n, v)*ds(2)
	Ftotal = Fpas + Fact + Fload
	Jac = derivative(Ftotal, w, dw)
	
	# Solve variational problem
	# Optimization options for the form compiler
	ffc_options = {"optimize": True}
	solver_options = {"newton_solver":{"maximum_iterations":100, "absolute_tolerance":1e-8, "relative_tolerance":1e-7}}

	self.w = w
	self.w_n = w_n
	self.Ftotal = Ftotal
	self.Jac = Jac
	self.bcs = bcs
	self.mesh = mesh
	self.t_a = t_a
	self.activeforms = activeforms
	self.f0 = f0
	self.QDG = QDG

    def Solver(self):

        solverparams_ME = {"Jacobian": self.Jac,
                "F": self.Ftotal,
                "w": self.w,
                "boundary_conditions": self.bcs,
                "Type": 0,
                "mesh": self.mesh,
		"mode": 1
               }
        solver_ME = NSolver(solverparams_ME)

	return solver_ME

    def UpdateVar(self):

	self.w_n.assign(self.w)

    def GetDisplacement(self):

	u, p = self.w.split(deepcopy=True)
	u.rename("u_", "u_")

	return u

    def GetSActive(self):
 
	Sactive = self.activeforms.PK2StressTensor()
	Sactive_ = project(self.f0[i]*Sactive[i,j]*self.f0[j], self.QDG)
	Sactive_.rename("Sact", "Sact")

        return Sactive_

  


