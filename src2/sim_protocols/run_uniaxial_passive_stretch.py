# Simulate Uniaxial Passive Stretching of the tissue
from dolfin import *
from mshr import *
import fenicstools as ft
import numpy as np
from ..mechanics.forms_MRC2 import Forms 


def run_uniaxial_test(IODet, SimDet):

    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]['representation'] = 'uflacs'
    parameters["form_compiler"]["quadrature_degree"] = 4
    # Is displacement control

    outputfolder = IODet["outputfolder"]
    folderName = IODet["folderName"] + IODet["caseID"] + '/' 

    length = 1
    width = 1
    nelem = 1
    nload = 50
    maxdisp = 1.0
    isdispctrl = SimDet["isdispctrl"]
    if("length" in SimDet.keys()):
		length = SimDet["length"]
    if("width" in SimDet.keys()):
		width = SimDet["width"]
    if("nelem" in SimDet.keys()):
		nelem = SimDet["nelem"]
    if("nload" in SimDet.keys()):
		nload = SimDet["nload"]
    if("maxdisp" in SimDet.keys()):
		maxdisp = SimDet["maxdisp"]


    mat_input = {}
    if("Passive model" in SimDet.keys()):
		matmodel = SimDet["Passive model"]
		mat_input.update({"material model": matmodel})
    if("Passive params" in SimDet.keys()):
		matparams= SimDet["Passive params"]
		mat_input.update({"material params": matparams})

    
    # Create mesh 
    mesh = BoxMesh(Point(0,0,0), Point(length,width,width), nelem, nelem, nelem)
    print("Number of Elements = "+str(mesh.num_cells()))
    File(outputfolder + folderName + "mesh.pvd") << mesh
    N = FacetNormal(mesh)
    
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
    
    File(outputfolder + folderName + "subdomain.pvd") <<  sub_domains
    
    # Set Fiber
    f0 = Expression(("1.0","0.0","0.0"), degree=1)
    s0 = Expression(("0.0","1.0","0.0"), degree=1)
    n0 = Expression(("0.0","0.0","1.0"), degree=1)
    
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

    params.update(mat_input)
    uflforms = Forms(params)
    
    SEF = uflforms.PassiveMatSEF()
    PK1pas = uflforms.PK1()
    Fe = uflforms.Fmat()
    n = uflforms.J()*inv(Fe.T)*N    
    
    # Boundary conditions
    uright = Expression("val", val = 0, degree=0)
    bcleft = DirichletBC(W.sub(0).sub(0), Constant((0.0)), sub_domains, 1)
    bcright = DirichletBC(W.sub(0).sub(0), uright, sub_domains, 2)
    bcback = DirichletBC(W.sub(0).sub(1), Constant((0.0)), sub_domains, 4)
    bcbot = DirichletBC(W.sub(0).sub(2), Constant((0.0)), sub_domains, 6)
    if(isdispctrl):
    	bcs = [bcleft, bcright, bcback, bcbot]
    else:
    	bcs = [bcleft, bcback, bcbot]

    pload = Expression("val", val = 0.0, degree=0)
    
    # Define Weak Form
    Fpas = derivative(SEF, w, wtest)*dx
    Fload =  -inner(pload*n, v)*ds(2)
    if(isdispctrl):
    	Ftotal = Fpas 
    else:
    	Ftotal = Fpas + Fload
    Jac = derivative(Ftotal, w, dw)
    
    # Solve variational problem
    # Optimization options for the form compiler
    ffc_options = {"optimize": True}
    solver_options = {"newton_solver":{"maximum_iterations":100, "absolute_tolerance":1e-9, "relative_tolerance":1e-8}}
    
    maxload = 5	
    
    uright.val = 0
    w.vector()[:] = 0
    load_array = []
    
    # Preload 
    lbda_arr = [1.0]
    pload_arr = [0.0]
    for load in np.arange(0, nload):
    	if(isdispctrl):
    		print("Applied disp = ", uright.val)
    		dispstep = maxdisp/nload
    		uright.val += dispstep
    		solve(Ftotal == 0, w, bcs, J=Jac, form_compiler_parameters=ffc_options, solver_parameters=solver_options)
    		#pload_ = assemble((PK1pas[0,0])*ds(2))
    		#lbda = project(f0[i]*Fe[i,j]*f0[j], QDG).vector().array()[:]

		kappa = uflforms.matparams["kappa"]
		C = Fe.T*Fe
		I = Identity(u)
		H = kappa*I + (1-3*kappa)*as_tensor(f0[i]*f0[j], (i,j)) 
		lbda = project(sqrt(tr(H*C)), QDG).vector().array()[:]

		pload_ = assemble(tr(H*PK1pas)*ds(2))

    	else:
    		print("Applied load = ", pload.val)
    		loadstep = maxload/nload
    		pload.val += loadstep
    		solve(Ftotal == 0, w, bcs, J=Jac, form_compiler_parameters=ffc_options, solver_parameters=solver_options)
    		lbda = project(f0[i]*Fe[i,j]*f0[j], QDG).vector().array()[:]
    		pload_ = pload.val/lbda[0]
    
    	pload_arr.append(pload_)
    	lbda_arr.append(lbda[0])

    return pload_arr, lbda_arr, uflforms.matparams

