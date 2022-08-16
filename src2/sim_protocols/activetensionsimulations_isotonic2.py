from dolfin import *
from mshr import *
import fenicstools as ft
import numpy as np
from matplotlib import pylab as plt
from forms_MRC2 import Forms 
from activeforms_MRC2 import activeForms 
#from Holzapfel_uniaxial_stretching import Holzapfel_analytical_uniaxial as H_analytical
#from Time_varying import Time_varying as TV_analytical

#set_log_level(LogLevel.TRACE)
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]['representation'] = 'uflacs'
parameters["form_compiler"]["quadrature_degree"] = 4

"""# **Set Up Mesh**"""

# Create mesh 
#domain = Box(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0))
#mesh = generate_mesh(domain, 2)
mesh = BoxMesh(Point(0,0,0), Point(1,1,1), 1, 1, 1)
print("Number of Elements = "+str(mesh.num_cells()))
plot(mesh, title="Mesh")
File("mesh.pvd") << mesh

# Mark Facet
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1.0 - DOLFIN_EPS and on_boundary

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

class Front(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 1.0 - DOLFIN_EPS and on_boundary

class Back(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] > 1.0 - DOLFIN_EPS and on_boundary

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

plot(sub_domains, title="SubDomains")
File("subdomain.pvd") <<  sub_domains

# Set Fiber
f0 = Expression(("1.0","0.0","0.0"), degree=1)
s0 = Expression(("0.0","1.0","0.0"), degree=1)
n0 = Expression(("0.0","0.0","1.0"), degree=1)

# Define Integration domain
dx =  dx(mesh, metadata = {"integration_order": 2})
ds =  ds(mesh, subdomain_data = sub_domains)

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


a_pas = Expression(("val"), val = 50.059, degree=0)

params= {"mesh": mesh,
     "facetboundaries": sub_domains,
     "facet_normal": FacetNormal(mesh),
     "mixedfunctionspace": W,
     "mixedfunction": w,
     "displacement_variable": u, 
     "pressure_variable": p,
     "fiber": f0,
     "sheet": s0,
     "sheet-normal": n0,
     "growth_tensor": None,
     "incompressible" : True,
     #"material model": {"Name": "Guccione"}, 
     #"material params": {"bff": Constant(29.9), "bfx": Constant(18.8), "bxx": Constant(18.8), "C_param": Constant(100)},
     "material model": {"Name": "HolzapfelOgden"},
     #"material params": {"a": a_pas, "b": 8.023, "a_f": 18.472, "b_f": 16.026, "a_s": 12.481, "b_s": 11.120, "a_fs": 0.216, "b_fs": 11.436}
     "material params": {"a": 10.059, "b": 8.023, "a_f": a_pas, "b_f": 8.026, "a_s": 12.481, "b_s": 11.120, "a_fs": 0.216, "b_fs": 11.436}
    }


uflforms = Forms(params)
Tmax = Expression(("val"), val = 100.0e3, degree=0)

t_a = Expression(("ta"), ta = 0.0, degree=0)
activeparams = {"mesh": mesh,
                "dx":dx,
		"deg":4,
                "facetboundaries": sub_domains,
                "facet_normal": FacetNormal(mesh),
                "displacement_variable": u, 
                "pressure_variable": p,
                "fiber": f0,
                "sheet": s0,
                "sheet-normal": n0,
                "t_a": t_a,
                "mesh": mesh, 
                "Threshold_Potential": 0.9,
                "growth_tensor": None,
	        "material model": {"Name": "Time-varying"},
     	        "material params": {"tau" : 25, "t_trans" : 300, "B" : 4.75,  "t0" : 250,  "deg" : 4, "l0" : 0.50, "Tmax" : Tmax, \
				    "Ca0" : 4.35, "Ca0max" : 4.35, "lr" : 1.9},
	        #"material model": {"Name": "Guccione"},
	
		"HomogenousActivation": True,
               }

activeforms = activeForms(activeparams)
potential = Function(FunctionSpace(mesh,'CG',1))
activeforms.update_activationTime(potential_n = potential, comm = mesh.mpi_comm())


# Define UFL forms
d = u.ufl_domain().geometric_dimension()
I = Identity(d)
F = I + grad(u)
F = variable(F)
J = det(F);
E = 0.5*(F.T*F - I);
N = FacetNormal(mesh)
n = J*inv(F.T)*N

SEF = uflforms.PassiveMatSEF()
PK1pas = uflforms.PK1()
PK1sef = uflforms.PK1sef()
PK1vol = uflforms.PK1vol()
PK1W1 = uflforms.PK1W1()
PK1W4f = uflforms.PK1W4f()
PK1W4s = uflforms.PK1W4s()
PK1W8 = uflforms.PK1W8()
Fe = uflforms.Fmat()


# Define Stress
Sactive = activeforms.PK2StressTensor()
PK1act = F*Sactive


# Boundary conditions
uright = Expression("val", val = 0, degree=0)
bcleft = DirichletBC(W.sub(0).sub(0), Constant((0.0)), sub_domains, 1)
bcright = DirichletBC(W.sub(0).sub(0), uright, sub_domains, 2)
bcback = DirichletBC(W.sub(0).sub(1), Constant((0.0)), sub_domains, 4)
bcbot = DirichletBC(W.sub(0).sub(2), Constant((0.0)), sub_domains, 6)
bcs = [bcleft, bcback, bcbot]
pload = Expression("val", val = 0.0, degree=0)

# Define Weak Form
Fact = inner(F*Sactive, grad(v))*dx
Fpas = derivative(SEF, w, wtest)*dx
Fload =  -inner(pload*n, v)*ds(2)
Ftotal = Fpas + Fact + Fload
Jac = derivative(Ftotal, w, dw)

# Solve variational problem
# Optimization options for the form compiler
ffc_options = {"optimize": True}
solver_options = {"newton_solver":{"maximum_iterations":100, "absolute_tolerance":1e-9, "relative_tolerance":1e-8}}
dispfile = File("u.pvd") 


Tmaxarray = [10e3,  30e3,  50e3,  70e3, 510e3, 5000e3]#, 500e3, 600e3]#, 700e3, 900e3]#, 100e3, 200e3, 300e3]

for Tmax_ in Tmaxarray:

	Tmax.val = Tmax_	

	ploadarray = [1000]#, 10, 10]#np.linspace(100,5000,5)#[100,1000, 5000, 10000, 15000]
	a_pas_array = [100]#, 5000, 10000]
	
	for (a_pas_val, pl) in zip(a_pas_array, ploadarray):
	
		pload.val = 0
		w.vector()[:] = 0
		t_a.ta = 0
		a_pas.val = a_pas_val

		probepts = np.array([[1.0, 0.5, 0.5]])
		xdisp = ft.Probes(probepts.flatten(), V)
	
		# Preload 
		nload = 10000
		lbda_arr = []
		active_load_array = []
		passive_load_array = []
		passiveW1_load_array = []
		passiveW4f_load_array = []
		passiveW4s_load_array = []
		passiveW8_load_array = []
		passivevol_load_array = []
		total_load_array = []

		for load in np.arange(0, nload):
			pload.val += pl
			print("Applied load = ", pload.val)
			solve(Ftotal == 0, w, bcs, J=Jac, form_compiler_parameters=ffc_options, solver_parameters=solver_options)
			lbda = project(f0[i]*Fe[i,j]*f0[j], QDG).vector().array()[:]
			pload_ = pload.val/lbda[0]

			if(lbda[0] > 1.15):
				break;
			#xdisp(w.sub(0))

			# Plot Deformation
			#dispfile << w.sub(0)
		
		# Contraction
		dt = 1.0
		ntpt = 500
		for tpt in np.arange(0, ntpt):
			t_a.ta += dt
			print("Active contraction =", t_a.ta)
			solve(Ftotal == 0, w, bcs, J=Jac, form_compiler_parameters=ffc_options, solver_parameters=solver_options)
		
			#if(t_a.ta%50 == 0):
			#	dispfile << w.sub(0)
		
		        PK1active_ = project(f0[i]*(PK1act[i,j])*f0[j], QDG)
			pload_active = assemble((PK1active_)*ds(2))*lbda[0]
			active_load_array.append(pload_active)

			PK1passive_ = project(f0[i]*(PK1pas[i,j])*f0[j], QDG)
			pload_passive = assemble((PK1passive_)*ds(2))*lbda[0]
			passive_load_array.append(pload_passive)

			pload_passiveW1 = assemble((project(f0[i]*(PK1W1[i,j])*f0[j], QDG))*ds(2))*lbda[0]
			passiveW1_load_array.append(pload_passiveW1)

			pload_passiveW4f = assemble((project(f0[i]*(PK1W4f[i,j])*f0[j], QDG))*ds(2))*lbda[0]
			passiveW4f_load_array.append(pload_passiveW4f)

			pload_passiveW4s = assemble((project(f0[i]*(PK1W4s[i,j])*f0[j], QDG))*ds(2))*lbda[0]
			passiveW4s_load_array.append(pload_passiveW4s)

			pload_passiveW8 = assemble((project(f0[i]*(PK1W8[i,j])*f0[j], QDG))*ds(2))*lbda[0]
			passiveW8_load_array.append(pload_passiveW8)

			pload_passivevol = assemble((project(f0[i]*(PK1vol[i,j])*f0[j], QDG))*ds(2))*lbda[0]
			passivevol_load_array.append(pload_passivevol)

			total_load_array.append(pload_active + pload_passiveW1  + pload_passiveW4f +  pload_passiveW4s + pload_passiveW8 + pload_passivevol)

			xdisp(w.sub(0))

			lbda = project(f0[i]*Fe[i,j]*f0[j], QDG).vector().array()[:]
			lbda_arr.append(lbda[0]*activeforms.matparams["lr"])

		
		ux = xdisp.array()[0,:]
		tpt = np.arange(0, len( xdisp.array()[0,:]))
		dudx = np.gradient(ux - ux[0], dt)

		plt.figure(2)
		plt.plot(tpt, ux - ux[0], label= 'load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))

		plt.figure(4)
		plt.plot(tpt, ux, label= 'load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))

		plt.figure(5)
		plt.plot(tpt, lbda_arr, label='load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))

		plt.figure(3)
		plt.plot(tpt, active_load_array, '--', label='A load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))
		#plt.plot(tpt, passive_load_array, '-', label='P load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))
		#plt.plot(tpt, total_load_array, '-', label='T load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))

		plt.figure(6)
		plt.plot(tpt, passiveW1_load_array, '-*', label='W1 load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))


		plt.figure(7)
		plt.plot(tpt, passiveW4f_load_array, '-', label='W4f load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))
		#plt.plot(tpt, passiveW4s_load_array, '-', label='W4s load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))
		#plt.plot(tpt, passiveW8_load_array, '-', label='W8 load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))


		plt.figure(8)
		plt.plot(tpt, passivevol_load_array, '-', label='Vol load = ' + str(pl) + ' a pas = ' + str(a_pas.val) + ' Tmax = ' + str(Tmax.val))
		
	
	
plt.figure(2)
plt.xlabel("tpt")
plt.ylabel("ux - ux0")
plt.legend()
plt.savefig("isotonic_dispx3.png")

plt.figure(3)
plt.xlabel("tpt")
plt.ylabel("active stress")
#plt.ylim([-500,500])
plt.legend()
plt.savefig("isotonic_dispx3active.png")

plt.figure(4)
plt.xlabel("tpt")
plt.ylabel("ux")
plt.legend()
plt.savefig("isotonic_dispx3a.png")

plt.figure(5)
plt.xlabel("tpt")
plt.ylabel("ls")
plt.legend()
plt.savefig("isotonic_dispx3b.png")

plt.figure(6)
plt.xlabel("tpt")
plt.ylabel("W1 passive stress")
plt.legend()
plt.savefig("isotonic_dispx3c.png")

plt.figure(7)
plt.xlabel("tpt")
plt.ylabel("W4f passive stress")
plt.legend()
plt.savefig("isotonic_dispx3d.png")

plt.figure(8)
plt.xlabel("tpt")
plt.ylabel("Vol passive stress")
plt.legend()
plt.savefig("isotonic_dispx3e.png")




