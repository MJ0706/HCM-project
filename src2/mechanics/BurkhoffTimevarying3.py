# Burkhoff et al. 1993: Time varying elastance type active contraction model with length dependency for peak tension but not relaxation
# This modification account for nonuniform spatial activation 
# to is a function of space, FE model using gauss points 

from dolfin import *                
import math
import numpy as np

class BurkhoffTimevarying(object):

    def __init__(self, params):

        self.parameters = self.default_parameters()
        self.parameters.update(params)

        deg = self.parameters["deg"] 
        bivMesh = self.parameters["mesh"]
        mesh = bivMesh     


	print "Using Burkhoff model"

        self.t_init = self.parameters["t_init"]
        self.isActive = self.parameters["isActive"]
         
        print "t_init vector is :" 
        print self.t_init.vector().array()

        print "isActive vector is :" 
        print self.isActive.vector().array()

        isActive_write_elem = FunctionSpace(mesh, "CG", 1)
        self.isActive_write_elem = isActive_write_elem

        t_init_write_elem = FunctionSpace(mesh, "CG", 1)
        self.t_init_write_elem = t_init_write_elem

    def default_parameters(self):
        return {"material params": {"tau" : 30,
				"t_trans" : 320,
				"B" : 4.75, #4.75 no strech dependency seen, 
				"t0" : 275, 
                		"deg" : 4,
                		"l0" : 1.6, 
        			"Tmax" : 135.7e3,
				"Ca0" : 4.35,
				"Ca0max" : 4.35,
				"lr" : 1.85,
				"kap" : 0.0,
				}
                	}

    def Getmatparam(self):
    
    	return self.parameters["material params"]


    def structTensor(self):
	k =  self.parameters["material params"]["kap"]	
	F = self.parameters["Fmat"]
	f0 = self.parameters["fiber"]
	Mij = as_tensor(f0[i]*f0[j],(i,j))

	I = Identity(F.ufl_domain().geometric_dimension()) 
	
	stensor = k*I +(1- 3*k)*Mij

	return stensor

    

    def ECa(self):

	Ca0max = self.parameters["material params"]["Ca0max"]
	B = self.parameters["material params"]["B"]

	ls_l0 = self.ls_l0()
	denom = sqrt(exp(B*(ls_l0)) - 1)

	ECa = Ca0max/denom
	
	return ECa

    def ls_l0(self):

	lmbda = self.lmbda()
	lr = self.parameters["material params"]["lr"]
	l0 = self.parameters["material params"]["l0"]

	F = self.parameters["Fmat"]
	Cmat = F.T*F

	stensor = self.structTensor()
	tensor2 = Cmat*stensor
	I4 = tr(tensor2)
	ls = sqrt(I4)*1.85
	#ls = lmbda*lr
	ls_l0 = conditional(le(ls, l0+.002), 0.002, ls - l0)

	return ls_l0
	

    def lmbda1(self):

	F = self.parameters["Fmat"]
	f0 = self.parameters["fiber"]
	Cmat = F.T*F
	lmbda = sqrt(dot(f0, Cmat*f0))

	return lmbda

    def lmbda(self):

	F = self.parameters["Fmat"]
	Cmat = F.T*F
	stensor = self.structTensor()
	tensor2 = Cmat*stensor
	I4 = tr(tensor2)
	lmbda = sqrt(I4)

	return lmbda



    def Ct(self):

	F = self.parameters["Fmat"]
	f0 = self.parameters["fiber"]
	t0 = self.parameters["material params"]["t0"]
	t_a = self.parameters["t_a"]

	Cmat = F.T*F
	lmbda = sqrt(dot(f0, Cmat*f0))

	w1 = self.w1()
	w2 = self.w2()

	Ct = w1+w2

	return Ct

    def w1(self):

	t0 = self.parameters["material params"]["t0"]
	t_a = self.parameters["t_a"] # current time 

        if("t_trans" in self.parameters["material params"].keys()):
        	t_trans = self.parameters["material params"]["t_trans"]
	else:
        	t_trans = 1.5*t0


        #isHomogenousActivation = self.parameters["HomogenousActivation"]

        #if isHomogenousActivation: # activation at 10 ms 
        #    t_init = self.t_init
        #    t_init.vector()[:] = 0.0*np.ones(len(t_init.vector().array()))
        #else:
        #    t_init = self.t_init # time of activation 

        t_init = self.t_init
        t_since_activation = t_a - t_init 

        xp4 = conditional( gt(t_since_activation, Constant(0.0)), 1.0, 0.0 )
        xp5 = conditional( lt(t_since_activation, Constant(9998.0)), 1.0, 0.0 )

        if("trans0" in self.parameters["material params"].keys()):
		trans0 = self.parameters["material params"]["trans0"]
		xp1 = conditional( lt(t_since_activation,trans0), 1.0, 0.0 ) # HACK LCLEE
	else:
		xp1 = conditional( lt(t_since_activation,t_trans), 1.0, 0.0 )

	w1 = xp5*xp4*xp1*0.5*(1 - cos(pi*t_since_activation/t0))

	return w1
	
    def w2(self):

	t0 = self.parameters["material params"]["t0"]
	t_a = self.parameters["t_a"]

        if("t_trans" in self.parameters["material params"].keys()):
        	t_trans = self.parameters["material params"]["t_trans"]
	else:
        	t_trans = 1.5*t0


	tr = self.parameters["material params"]["tau"]

        t_init = self.t_init # time of activation 
        t_since_activation = t_a - t_init 

      	xp2 = conditional(le(t_trans,t_since_activation), 1.0, 0.0)
        if("trans0" in self.parameters["material params"].keys()):
		A = 1.0;
	else:
		A = 0.5*(1 - cos(pi*t_trans/t0))
	w2 = xp2*A*exp(-1.0*(t_since_activation - t_trans)/tr)

	return w2

    def PK2Stress(self):

	Ca0 = self.parameters["material params"]["Ca0"]
	Tmax = self.parameters["material params"]["Tmax"]
	
	
	Ct = self.Ct()
	ECa = self.ECa()

	Sact = (Tmax*Ct*Ca0**2.0)/(Ca0**2.0 + ECa**2.0)

	return Sact


    def I4f_k(self):

	kappa = self.parameters["material params"]["kap"]

	F = self.parameters["Fmat"]
	Cmat = F.T*F
	I1 = tr(Cmat)

	f0 = self.parameters["fiber"]
	I4 = inner(Cmat*f0, f0)

	I4fk = kappa*I1 + (1.0-3.0*kappa)*I4

	I4f = sqrt(I4fk)

	return I4f

