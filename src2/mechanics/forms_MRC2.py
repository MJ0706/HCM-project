from dolfin import *
import sys

class Forms(object):

    	def __init__(self, params):        

		self.parameters = self.default_parameters()
         	self.parameters.update(params)

		Matmodel =  self.parameters["material model"]["Name"]
		assert (Matmodel == "Guccione" or Matmodel == "HolzapfelOgden"), "Material model not implemented"

		self.parameters.update({"F": self.Fe()})
		if(Matmodel == "Guccione"):
        		from GuccionePas import GuccionePas as Passive

		if(Matmodel == "HolzapfelOgden"):
        		from holzapfelogden import HolzapfelOgden as Passive

	        self.passiveforms = Passive(self.parameters)
		self.matparams = self.passiveforms.Getmatparam()

        def default_parameters(self):
	        return {"material model": {"Name": "Guccione"},
			"Kappa": 1e5,
			"incompressible" : True,
			};
	
	def PassiveMatSEF(self):

		Wp = self.passiveforms.PassiveMatSEF() + self.Wvolumetric()

		return Wp

	def PK1(self):

		PK1 =  self.PK1volumetric() + self.passiveforms.PK1() 
		return PK1

	def PK2(self):

		PK2 = inv(self.Fmat())*self.PK1()
		return PK2

	def sigma(self):

		sigma = 1.0/self.J()*self.PK1()*transpose(self.Fmat())
		return sigma


	def Fmat(self):
	
		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
	    	return F

	def Fe(self):
		Fg = self.parameters["growth_tensor"]
	    	F = self.Fmat() 
		if (Fg is None):
                	Fe = F
        	else:
                	Fe = as_tensor(F[i,j]*inv(Fg)[j,k], (i,k))

		return Fe
	

	
	def Emat(self):
	 
		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)
	    	F = self.Fe()
	    	return 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))


	def J(self):

	    	F = self.Fe()
		return det(F)	


	def Wvolumetric(self):

		isincomp = self.parameters["incompressible"]
		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		F = dolfin.variable(F)
		J = det(F)

		if(isincomp):
			p = self.parameters["pressure_variable"]
			Wvolumetric = -1.0*p*(J - 1.0)
		else:
			Kappa = self.parameters["Kappa"]
			Wvolumetric =  Kappa/2.0*(J - 1.0)**2.0

		return Wvolumetric

	def PK1volumetric(self):

		isincomp = self.parameters["incompressible"]
		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
		I = Identity(d)
		F = I + grad(u)
		F = dolfin.variable(F)
		J = det(F)
		

		if(isincomp):

			p = self.parameters["pressure_variable"]
			Wvolumetric = -1.0*p*(J - 1.0)
			PK1volumetric = dolfin.diff(Wvolumetric, F)
		else:
			Kappa = self.parameters["Kappa"]
			Wvolumetric =  Kappa/2.0*(J - 1.0)**2.0
			PK1volumetric = dolfin.diff(Wvolumetric, F)

		return PK1volumetric

		
	def LVcavityvol(self):
	
		u = self.parameters["displacement_variable"]
		N = self.parameters["facet_normal"]
		mesh = self.parameters["mesh"]
	    	X = SpatialCoordinate(mesh)
		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
	    	
	    	F = self.Fmat()
	    	
	    	vol_form = -Constant(1.0/3.0) * inner(det(F)*dot(inv(F).T, N), X + u)*ds(self.parameters["LVendoid"])
	    	
	    	return assemble(vol_form, form_compiler_parameters={"representation":"uflacs"})

	def RVcavityvol(self):
	
		u = self.parameters["displacement_variable"]
		N = self.parameters["facet_normal"]
		mesh = self.parameters["mesh"]
	    	X = SpatialCoordinate(mesh)
		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
	    	
	    	F = self.Fmat()
	    	
	    	vol_form = -Constant(1.0/3.0) * inner(det(F)*dot(inv(F).T, N), X + u)*ds(self.parameters["RVendoid"])
	    	
	    	return assemble(vol_form, form_compiler_parameters={"representation":"uflacs"})


        def LVcavitypressure(self):
    
		W = self.parameters["mixedfunctionspace"]
		w = self.parameters["mixedfunction"]
		mesh = self.parameters["mesh"]

		comm = W.mesh().mpi_comm()
    	        dofmap =  W.sub(self.parameters["LVendo_comp"]).dofmap()
        	val_dof = dofmap.cell_dofs(0)[0]

	        # the owner of the dof broadcasts the value
	        own_range = dofmap.ownership_range()
    
	        try:
	            val_local = w.vector()[val_dof][0]
	        except IndexError:
	            val_local = 0.0


    		pressure = MPI.sum(comm, val_local)

        	return pressure



        def RVcavitypressure(self):
    
		W = self.parameters["mixedfunctionspace"]
		w = self.parameters["mixedfunction"]
		mesh = self.parameters["mesh"]

		comm = W.mesh().mpi_comm()
    	        dofmap =  W.sub(self.parameters["RVendo_comp"]).dofmap()
        	val_dof = dofmap.cell_dofs(0)[0]

	        # the owner of the dof broadcasts the value
	        own_range = dofmap.ownership_range()
    
	        try:
	            val_local = w.vector()[val_dof][0]
	        except IndexError:
	            val_local = 0.0


    		pressure = MPI.sum(comm, val_local)

        	return pressure


	def LVV0constrainedE(self):


		mesh = self.parameters["mesh"]
		u = self.parameters["displacement_variable"]
		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
		dsendo = ds(self.parameters["LVendoid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
		area = self.parameters["LVendo_area"]
		pendo = self.parameters["lv_volconst_variable"] 
		V0= self.parameters["lv_constrained_vol"] 

	    	X = SpatialCoordinate(mesh)
		x = u + X

	    	F = self.Fmat()
		N = self.parameters["facet_normal"]
        	n = cofac(F)*N

        	V_u = - Constant(1.0/3.0) * inner(x, n)
		Wvol = (Constant(1.0)/area * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)

		return Wvol


	def RVV0constrainedE(self):


		mesh = self.parameters["mesh"]
		u = self.parameters["displacement_variable"]
		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
		dsendo = ds(self.parameters["RVendoid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
		pendo = self.parameters["rv_volconst_variable"] 
		V0= self.parameters["rv_constrained_vol"] 

	    	X = SpatialCoordinate(mesh)
		x = u + X

	    	F = self.Fmat()
		N = self.parameters["facet_normal"]
        	n = cofac(F)*N

		area = assemble(Constant(1.0) * dsendo, form_compiler_parameters={"representation":"uflacs"})
        	V_u = - Constant(1.0/3.0) * inner(x, n)
		Wvol = (Constant(1.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)

		return Wvol


	def fiberstress_wokappa(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1

		f0 = self.parameters["fiber"]

		return f0[i]*Sca[i,j]*f0[j]

	def sheetstress(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1

		s0 = self.parameters["sheet"]

		return s0[i]*Sca[i,j]*s0[j]

	def normalstress(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1

		n0 = self.parameters["sheet-normal"]

		return n0[i]*Sca[i,j]*n0[j]


	def fiberstress_1(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1
		H = self.passiveforms.structTensor() 
		fstress = tr(Sca*H)

		return fstress

	def fiberstress_2(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()
		f0 = self.parameters["fiber"]
		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1
		H = self.passiveforms.structTensor() 
		fstress = f0[i]*Sca[i,j]*f0[j]

		return fstress

	def fiberstress_3(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()
		f0 = self.parameters["fiber"]
		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1
		H = self.passiveforms.structTensor() 
		fstress = f0[i]*Sca[i,j]*f0[j]
		kappa = self.passiveforms.parameters["material params"]["kappa"]

		return (1-3*kappa)*fstress

	def fiberstress_4(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()
		f0 = self.parameters["fiber"]
		Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1
		H = self.passiveforms.structTensor() 
		fstress = f0[i]*Sca[i,j]*f0[j]
		kappa = self.passiveforms.parameters["material params"]["kappa"]

		return (1-2*kappa)*fstress

	def fiberstrain_1(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		f0 = self.parameters["fiber"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))
		H = self.passiveforms.structTensor()
		fstrain = tr(Emat*H)

		return fstrain #f0[i]*Emat[i,j]*f0[j]

	def fiberstrain_2(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		f0 = self.parameters["fiber"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))
		#H = self.passiveforms.structTensor()
		#fstrain = tr(Emat*H)

		return f0[i]*Emat[i,j]*f0[j]

	def fiberstrain_3(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		f0 = self.parameters["fiber"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))
		#H = self.passiveforms.structTensor()
		#fstrain = tr(Emat*H)
		kappa = self.passiveforms.parameters["material params"]["kappa"]
		fstrain = f0[i]*Emat[i,j]*f0[j]
		return (1-3*kappa)*fstrain

	def fiberstrain_4(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		f0 = self.parameters["fiber"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))
		#H = self.passiveforms.structTensor()
		#fstrain = tr(Emat*H)
		kappa = self.passiveforms.parameters["material params"]["kappa"]
		fstrain = f0[i]*Emat[i,j]*f0[j]
		return (1-2*kappa)*fstrain

	def sheetstrain(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		s0 = self.parameters["sheet"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))


		return s0[i]*Emat[i,j]*s0[j]

	def normalstrain(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		n0 = self.parameters["sheet-normal"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))


		return n0[i]*Emat[i,j]*n0[j]

	def fiberstrainalmansi(self, F_ref):

		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		f0 = self.parameters["fiber"]
		b = F*F.T
		B = inv(b)



		emat = 0.5*(I-B)
		H = self.passiveforms.structTensor()
		hi = H*F.T
		h = F*hi 
		fstrain = tr(emat*H)
		#fstrain = tr(emat*h.T)

		return fstrain



	def fiberstrain_wokappa(self, F_ref):


		u = self.parameters["displacement_variable"]
	    	d = u.ufl_domain().geometric_dimension()
	    	I = Identity(d)

                F = self.Fe()*inv(F_ref)
		f0 = self.parameters["fiber"]
		Emat = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))

		return f0[i]*Emat[i,j]*f0[j]


	def fiberwork(self, F_ref):

                F = self.Fe()*inv(F_ref)

		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T 
		Sca = inv(F)*PK1

		f0 = self.parameters["fiber"]
		Emat = self.Emat()

		H = self.passiveforms.structTensor()
		strain = tr(H*Emat)
		stress = tr(H*Sca)

		#f = f0[i]*Sca[i,j]*f0[j]
		#s = f0[i]*Emat[i,j]*f0[j]

		return strain*stress #dot(f,s)


	def IMP(self):
		u = self.parameters["displacement_variable"]
                J = self.J()
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]

		#F = self.Fmat()
		F = self.Fe()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T

		s = F*s0/sqrt(inner(F*s0, F*s0))
		n = F*n0/sqrt(inner(F*n0, F*n0))

		Ipressure = s[i]*Tca[i,j]*s[j]

                return Ipressure

	def IMP2(self):
		u = self.parameters["displacement_variable"]
                J = self.J()
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]

		#F = self.Fmat()
		F = self.Fe()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T

		s = F*s0/sqrt(inner(F*s0, F*s0))
		n = F*n0/sqrt(inner(F*n0, F*n0))

		Ipressure = 0.5*(s[i]*Tca[i,j]*s[j] + n[i]*Tca[i,j]*n[j])

                return Ipressure

	def IMPendo(self):

		u = self.parameters["displacement_variable"]
		N = self.parameters["facet_normal"]
		F = self.Fe()
		PK1 = self.PK1()

		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
		n = F*N/sqrt(inner(F*N, F*N))

		Ipressure = -n[i]*PK1[i,j]*N[j]*ds(self.parameters["LVendoid"])

		return Ipressure

	def IMPepi(self):

		u = self.parameters["displacement_variable"]
		N = self.parameters["facet_normal"]
		F = self.Fe()
		PK1 = self.PK1()

		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
		n = F*N/sqrt(inner(F*N, F*N))

		Ipressure = -n[i]*PK1[i,j]*N[j]*ds(self.parameters["epiid"])

		return Ipressure

	def areaendo(self):
		u = self.parameters["displacement_variable"]
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]


		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		J = det(F)

		s = F*s0/sqrt(inner(F*s0, F*s0))

		N = self.parameters["facet_normal"]
		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
		n = F*N/sqrt(inner(F*N, F*N))

		return (J*inv(F.T)*N)[i]*n[i]*ds(self.parameters["endoid"])

	def areaepi(self):
		u = self.parameters["displacement_variable"]
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]


		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		J = det(F)

		s = F*s0/sqrt(inner(F*s0, F*s0))

		N = self.parameters["facet_normal"]
		ds = dolfin.ds(subdomain_data = self.parameters["facetboundaries"])
		n = F*N/sqrt(inner(F*N, F*N))

		
		return (J*inv(F.T)*N)[i]*n[i]*ds(self.parameters["epiid"])



	def globalfiberstress(self):

		u = self.parameters["displacement_variable"]
		f0 = self.parameters["fiber"]

		F = self.Fmat()
		J = self.J()

		W = self.passiveforms.TotalPasMatSEF(diff=0) + self.Wvolumetric()

		#F = dolfin.variable(F)

		#PK1 = dolfin.diff(W,F)

		PK1 = self.PK1()
		mesh1 = self.parameters["mesh"]
		dx = dolfin.dx(mesh1) #(subdomain_data = self.parameters["facetboundaries"])
		#Tca = (1.0/J)*PK1*F.T
		#Sca = inv(F)*PK1
		kappa = self.passiveforms.parameters["material params"]["kappa"]

		d = u.geometric_dimension()
	    	I = Identity(d)
		H = kappa*I + (1-3*kappa)*as_tensor(f0[i]*f0[j], (i,j)) 
		fstress = tr(H*PK1)

		volume = assemble(Constant(1.0)*dx, form_compiler_parameters={"representation":"uflacs"})
		stress = assemble(fstress*dx, form_compiler_parameters={"representation":"uflacs"})/volume
		
		return stress #f0[i]*Sca[i,j]*f0[j]

	def globalfiberstretch(self):

		I4 = self.passiveforms.globalpassivefstretch()
		mesh1 = self.parameters["mesh"]
		dx = dolfin.dx(mesh1)

		volume = assemble(Constant(1.0)*dx, form_compiler_parameters={"representation":"uflacs"})
		stretch = assemble(I4*dx, form_compiler_parameters={"representation":"uflacs"})/volume
		
		return stretch #f0[i]*Sca[i,j]*f0[j]


	def FiberCauchyStress(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		#Sca = inv(F)*PK1
		H = self.passiveforms.structTensor() 
		hi = H*F.T
		h = F*hi
		fstress = tr(Tca*H)

		return fstress

	def FiberPk2Stress(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		#Tca = (1.0/J)*PK1*F.T
		Sca = inv(F)*PK1
		H = self.passiveforms.structTensor() 
		hi = H*F.T
		h = F*hi
		fstress = tr(Sca*H)

		return fstress


	def SheetCauchyStress(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		#Sca = inv(F)*PK1
		ef0 = self.parameters["sheet"]
		ef = (F*ef0)/sqrt(inner(F*ef0, F*ef0))

		fstress = dot(Tca*ef, ef)

		return fstress

	def SheetnormalCauchyStress(self):

		F = self.Fmat()
		J = self.J()
		PK1 = self.PK1()

		Tca = (1.0/J)*PK1*F.T
		#Sca = inv(F)*PK1
		ef0 = self.parameters["sheet-normal"]
		ef = (F*ef0)/sqrt(inner(F*ef0, F*ef0))

		fstress = dot(Tca*ef, ef)

		return fstress


	


