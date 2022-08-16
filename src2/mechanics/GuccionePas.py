from dolfin import *

class GuccionePas(object):

    	def __init__(self, params):        

		self.parameters = self.default_parameters()
         	self.parameters.update(params)


        def default_parameters(self):

     		return {"material params": { 
					  "bff": Constant(29.9), 
					  "bfx": Constant(13.3), 
					  "bxx": Constant(26.6), 
					  "Cparam": Constant(100)}}

	def Getmatparam(self):

		return self.parameters["material params"]



	def Emat(self):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)
		F = self.parameters["F"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))

		return Emat

	
	def PassiveMatSEF(self):

		Ea = self.Emat()
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]
 		bff = self.parameters["material params"]["bff"]
 		bfx = self.parameters["material params"]["bfx"]
 		bxx = self.parameters["material params"]["bxx"]
		isincomp = self.parameters["incompressible"]
	

		if(isincomp):
			p = self.parameters["pressure_variable"]

		C = self.parameters["material params"]["C_param"]


		Eff = inner(f0, Ea*f0)
		Ess = inner(s0, Ea*s0)
		Enn = inner(n0, Ea*n0)
		Efs = inner(f0, Ea*s0)
		Efn = inner(f0, Ea*n0)
		Ens = inner(n0, Ea*s0)
		Esf = inner(s0, Ea*f0)
		Enf = inner(n0, Ea*f0)
		Esn = inner(s0, Ea*n0)
	
		QQ = bff*Eff**2.0 + bxx*(Ess**2.0 + Enn**2.0 + Ens**2.0 +  Esn**2.0) + bfx*(Efs**2.0 + Esf**2.0 + Efn**2.0 + Enf**2.0)

		Wp = C/2.0*(exp(QQ) -  1.0)

		return Wp



	def PK1(self):

		u = self.parameters["displacement_variable"]

		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]
 		bff = self.parameters["material params"]["bff"]
 		bfx = self.parameters["material params"]["bfx"]
 		bxx = self.parameters["material params"]["bxx"]
		C = self.parameters["material params"]["C_param"]
		#p = self.parameters["pressure_variable"]

		#### For some reason to use dolfin.diff, you need to declare everything starting from u #############################
		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		F = dolfin.variable(F)
		J = det(F)

		Ea = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))

		Eff = inner(f0, Ea*f0)
		Ess = inner(s0, Ea*s0)
		Enn = inner(n0, Ea*n0)
		Efs = inner(f0, Ea*s0)
		Efn = inner(f0, Ea*n0)
		Ens = inner(n0, Ea*s0)
		Esf = inner(s0, Ea*f0)
		Enf = inner(n0, Ea*f0)
		Esn = inner(s0, Ea*n0)
		
		QQ = bff*Eff**2.0 + bxx*(Ess**2.0 + Enn**2.0 + Ens**2.0 +  Esn**2.0) + bfx*(Efs**2.0 + Esf**2.0 + Efn**2.0 + Enf**2.0)
		Wp = C/2.0*(exp(QQ) -  1.0) #- p*(J - 1.0)

                PK1 = dolfin.diff(Wp,F)               
                return PK1


