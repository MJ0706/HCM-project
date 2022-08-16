# Guccione2 et al. 1993: Simple time varying elastance type active contraction model 
# This modification account for nonuniform spatial activation 
# to is a function of space, FE model using gauss points 

from dolfin import *                
import math
import numpy as np

class GuccioneAct(object):

    def __init__(self, params):

        self.parameters = self.default_parameters()
        self.parameters.update(params)
        
        deg = self.parameters["deg"] 
        bivMesh = self.parameters["mesh"]
        mesh = bivMesh     

	print "Using Guccione Active model"

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
				    "B" : 4.75, 
				    "t0" : 275, 
                		    "deg" : 4,
                		    "l0" : 1.58,
        			    "Tmax" : 135.7e3,
				    "Ca0" : 4.35,
				    "Ca0max" : 4.35,
				    "m": 1048,
				    "b": -1675}

                };

    def Getmatparam(self):
    
    	return self.parameters["material params"]
 

    def PK2Stress(self):

	Ca0 = self.parameters["material params"]["Ca0"]
	Tmax = self.parameters["material params"]["Tmax"]

	Ct = self.Ct()
	ECa = self.ECa()

	Sact = (Tmax*Ct*Ca0**2.0)/(Ca0**2.0 + ECa**2.0)

	return Sact


#    Jay's original 
    def ECa(self):

	Ca0max = self.parameters["material params"]["Ca0max"]
	B = self.parameters["material params"]["B"]

	ls_l0 = self.ls_l0()
	denom = sqrt(exp(B*(ls_l0)) - 1)

	ECa = Ca0max/denom
	
	return ECa


    def ls_l0(self):

	lmbda = self.lmbda()
	ls = lmbda*1.85
	l0 = self.parameters["material params"]["l0"]
	ls_l0 = conditional(le(ls, l0+.002), 0.002, ls - l0)

	return ls_l0
	


    def tr(self):

	F = self.parameters["Fmat"]
	f0 = self.parameters["fiber"]
	b = self.parameters["material params"]["b"]
	m = self.parameters["material params"]["m"]
	Cmat = F.T*F
	lmbda = self.lmbda()
	ls = lmbda*1.85

	tr = m*ls + b

	return tr

    def Ct(self):

	F = self.parameters["Fmat"]
	f0 = self.parameters["fiber"]
	t0 = self.parameters["material params"]["t0"]
	b = self.parameters["material params"]["b"]
	m = self.parameters["material params"]["m"]
	t_a = self.parameters["t_a"]

	Cmat = F.T*F
	lmbda = sqrt(dot(f0, Cmat*f0))
	ls = lmbda*1.85

	tr = self.tr()
	xp1 = conditional(lt(t_a,t0), 1.0, 0.0)
	w1 = self.w1()
      	xp2 = conditional(le(t0,t_a), 1.0, 0.0)
	xp3 = conditional(lt(t_a,t0+tr), 1.0, 0.0)
	w2 = self.w2()

	Ct = 0.5*(1 - cos(w1+w2))

	return Ct


    def w1(self):

	t0 = self.parameters["material params"]["t0"]
	t_a = self.parameters["t_a"] # current time 

        t_init = self.t_init
        t_since_activation = t_a - t_init 

	xp1 = conditional( lt(t_since_activation,t0), 1.0, 0.0 )
        xp4 = conditional( gt(t_since_activation, Constant(0.0)), 1.0, 0.0 )
        xp5 = conditional( lt(t_since_activation, Constant(9998.0)), 1.0, 0.0 )

	w1 = xp5*xp4*xp1*pi*t_since_activation/t0

	return w1
	
    def w2(self):

	t0 = self.parameters["material params"]["t0"]
	t_a = self.parameters["t_a"]
	tr = self.tr()

        t_init = self.t_init # time of activation 
        t_since_activation = t_a - t_init 

      	xp2 = conditional(le(t0,t_since_activation), 1.0, 0.0)
	xp3 = conditional(lt(t_since_activation,t0+tr), 1.0, 0.0) 

	w2 = xp2*xp3*pi*(t_since_activation - t0 + tr)/tr

	return w2

    def lmbda(self):

	F = self.parameters["Fmat"]
	f0 = self.parameters["fiber"]
	Cmat = F.T*F
	lmbda = sqrt(dot(f0, Cmat*f0))

	return lmbda



