# Electrophysiology model using Aliev Panilfov model
from dolfin import *                
import numpy as np
from ..utils.nsolver import NSolver as NSolver
from fenicstools import *

class EPmodel(object):

    def __init__(self, params):

        self.parameters = self.default_parameters()
        self.parameters.update(params)

        self.mesh_ep = self.parameters["EPmesh"]
    	P1_ep = FiniteElement('CG', self.mesh_ep.ufl_cell(), 1, quad_scheme="default")
    	P1_ep._quad_scheme = 'default'
    	P2_ep = FiniteElement('DG', self.mesh_ep.ufl_cell(), 0, quad_scheme="default")
    	P2_ep._quad_scheme = 'default'

    	self.W_ep = FunctionSpace(self.mesh_ep, MixedElement([P1_ep, P2_ep]))
    	self.w_ep = Function(self.W_ep)
    	self.dw_ep = TrialFunction(self.W_ep) 
    	self.wtest_ep = TestFunction(self.W_ep)
        self.w_n_ep = Function(self.W_ep)

	#self.F_FHN, self.J_FHN, self.f_1, self.f_2 = self.Problem()
	self.F_FHN, self.J_FHN, self.fstim_array = self.Problem()


    def default_parameters(self):
        return {
		"ploc": [[0.0, 0.0, 0.0, 100.0, 1]],
            	"pacing_timing": [[1, 20.0]],
                };


    def calculateDmat(self, f0, mesh, mId):

    	d_iso = self.parameters["d_iso"] #0.01 #0.02 
    	d_ani = d_iso*self.parameters["d_ani_factor"] #0.08 #0.1 #0.2

    	d_Ani = d_ani

    	dParam = {}
    	if self.parameters["Ischemia"]:
    	    dParam['kNormal'] = d_ani
    	    dParam['kIschemia'] = d_iso

    	    d_Ani = defCPP_Matprop_DIsch(mesh = mesh, mId = mId, k = dParam)

    	D_iso = Constant(((d_iso,'0.0', '0.0'),('0.0',d_iso, '0.0'), ('0.0', '0.0', d_iso)))
    	Dij = d_Ani*f0[i]*f0[j] + D_iso[i,j]        
    	D_tensor = as_tensor(Dij, (i,j)) 

    	return D_tensor


    def MarkStimulus(self):

        mesh_ep = self.mesh_ep
        ploc = self.parameters["ploc"]
        EpiBCid_ep = FacetFunction('size_t', mesh_ep)
	EpiBCid_ep.set_all(0)

	for ploc_ in ploc:
    		class Omega_0(SubDomain):
    			def inside(self, x, on_boundary):
        			return ((x[0] - ploc_[0])**2 + (x[1] - ploc_[1])**2 + (x[2] - ploc_[2])**2) <= ploc_[3]**2
    		subdomain_0 = Omega_0()
    		subdomain_0.mark(EpiBCid_ep, int(ploc_[4]))

	# Find maximum number of pace label
	self.max_pace_label = int(max(np.array(ploc)[:,4]))

	return EpiBCid_ep


    def Problem(self):

        state_obj = self.parameters["state_obj"]
	f0_ep = self.parameters["f0"]
        AHAid_ep = self.parameters["AHAid"]
	matid_ep = self.parameters["matid"]
        facetboundaries_ep = self.parameters["facetboundaries"]
        EpiBCid_ep = self.MarkStimulus()
	mesh_ep = self.mesh_ep

	W_ep = self.W_ep
	w_n_ep = self.w_n_ep
	w_ep = self.w_ep
	dw_ep = self.dw_ep
	wtest_ep = self.wtest_ep

        phi0 = interpolate(Expression('0.0', degree=0), W_ep.sub(0).collapse())
        r0 = interpolate(Expression('0.0', degree=0), W_ep.sub(1).collapse())

        assign(w_n_ep, [phi0, r0])

        phi_n, r_n = split(w_n_ep)
        phi, r = split(w_ep)

        bcs_ep = []
        
        phi_test, r_test = split(wtest_ep)

        alpha = Constant(0.01) 
        g = Constant(0.002)
        b = Constant(0.15)
        c = Constant(8)
        mu1 = Constant(0.2)
        mu2 = Constant(0.3)

        def eps_FHN(phi, r):
            return (g + (mu1*r)/(mu2 + phi))

        def f_phi(phi, r):
            return -c*phi*(phi - alpha)*(phi - 1) - r*phi

        def f_r(phi, r):
            return eps_FHN(phi,r)*(- r -c*phi*(phi - b - 1.))


        fhn_timeNormalizer = 12.9
        k = state_obj.dt.dt/fhn_timeNormalizer 
        f_phi = f_phi(phi,r)
        f_r = f_r(phi,r)

    	D_tensor = self.calculateDmat(f0_ep, mesh = mesh_ep, mId = AHAid_ep)
        Dmat = D_tensor

	assert len(self.parameters["pacing_timing"]) == self.max_pace_label,\
		 "Number of pacing timing not equal to number of ploc labels"


	self.fstim_array = []
	for p in np.arange(0, self.max_pace_label):
		self.fstim_array.append(Expression('iStim', iStim=0.001, degree=1))

	#print len(self.parameters["pacing_timing"])	
	#stop

        #self.f_1 = Expression('iStim', iStim=0.001, degree=1) # 0.001
        #self.f_2 = Expression('iStim', iStim=0.001, degree=1)

        dx_ep = dolfin.dx(mesh_ep, subdomain_data=matid_ep)
        ds_ep = dolfin.ds(mesh_ep, subdomain_data=facetboundaries_ep)
    	ds_ep_epi = dolfin.ds(mesh_ep, subdomain_data=EpiBCid_ep)

        #pacing_integral1 = []
        #pacing_integral1.append( self.f_1*phi_test*ds_ep_epi(1) )

        #pacing_integral2 = []

        self.F_FHN = ( (phi - phi_n) / k)*phi_test*dx_ep + dot( Dmat*grad(phi), grad(phi_test))*dx_ep \
                - f_phi*phi_test*dx_ep \
                + ( (r - r_n) / k)*r_test*dx_ep \
                - f_r*r_test*dx_ep \
                #- sum( pacing_integral1 ) \
                #- sum( pacing_integral2 )

	label = 1
	for fstim in self.fstim_array:
		self.F_FHN -= fstim*phi_test*ds_ep_epi(label)
		label += 1

        self.J_FHN = derivative(self.F_FHN, w_ep, dw_ep)

	return self.F_FHN, self.J_FHN, self.fstim_array, #self.f_1, self.f_2

    def Solver(self):

        solverparams_FHN = {"Jacobian": self.J_FHN,
                "F": self.F_FHN,
                "w": self.w_ep,
                "boundary_conditions": [],
                "Type": 0,
                "mesh": self.mesh_ep,
		"mode": 1
               }
        solver_FHN = NSolver(solverparams_FHN)

        return solver_FHN

    def UpdateVar(self):

	# Update EP variable
	self.w_n_ep.assign(self.w_ep)

	# Update Stimulus variable
        state_obj = self.parameters["state_obj"]
	pace_time_array = self.parameters["pacing_timing"]

	for fstim, pace_time in zip(self.fstim_array, pace_time_array):
        	fstim.iStim = 0.0
		time = pace_time[0]
		duration = pace_time[1]
        	if state_obj.t >= time and ( state_obj.t <= time + duration): 
            		fstim.iStim = 0.3


    def Reset(self):

	self.w_ep.assign(self.w_n_ep)

    def getphivar(self):
	
        phi_, r_ = self.w_n_ep.split(deepcopy=True)
        phi_.rename("phi_", "phi_")

	return phi_
	
    def getrvar(self):
	
        phi_, r_ = self.w_n_ep.split(deepcopy=True)
        r_.rename("r_", "r_")

	return r_
	
    def interpolate_potential_ep2me_phi(self, V_me):

        lp = LagrangeInterpolator()
        lp.interpolate(V_me, self.getphivar())

        return V_me

    def reset(self):

        phi0 = interpolate(Expression('0.0', degree=0), self.W_ep.sub(0).collapse())
        r0 = interpolate(Expression('0.0', degree=0), self.W_ep.sub(1).collapse())
	assign(self.w_n_ep, [phi0, r0])
	assign(self.w_ep, [phi0, r0])

	return;

	

