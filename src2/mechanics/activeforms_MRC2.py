from dolfin import *
import sys
import math
import numpy as np

class activeForms(object):

    def __init__(self, params):        
    
    	self.parameters = self.default_parameters()
	#print self.parameters
     	self.parameters.update(params)
	#print self.parameters
    
    	Matmodel =  self.parameters["material model"]["Name"]
    	assert (Matmodel == "Guccione" or Matmodel == "Time-varying"), "Material model not implemented"
    
    	#self.parameters.update({"F": self.Fe()})
    	if(Matmodel == "Guccione"):
    		from GuccioneAct import GuccioneAct as Active
    
    	if(Matmodel == "Time-varying"):
    		from BurkhoffTimevarying3 import BurkhoffTimevarying as Active

        deg = self.parameters["deg"] 
        bivMesh = self.parameters["mesh"]
        mesh = bivMesh     

        Quadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=deg, quad_scheme="default")
        Quadelem._quad_scheme = 'default'
        self.Quadelem = Quadelem
        self.Quad = FunctionSpace(mesh, Quadelem)

        self.t_init = self.get_t_init()
        self.isActive = self.get_isActive()
 
	self.parameters.update({"t_init": self.t_init,
				"isActive": self.isActive,
				"Fmat": self.Fmat()
			       })
    
        self.activeforms = Active(self.parameters)
    	self.matparams = self.activeforms.Getmatparam()
    
    def default_parameters(self):
         return {"material model": {"Name": "Time-varying"}
    		};

    def get_isActive(self):
    
        deg = self.parameters["deg"]
        mesh = self.parameters["mesh"]
         
        Quad = self.Quad
        isActive = Function(Quad)
    
        isActive_array = isActive.vector().array()
        isActive_array = np.zeros(len(isActive_array))

        isActive.vector()[:] = isActive_array

        return isActive

    def Get_t_a(self):

	return self.parameters["t_a"].ta

    def update_activationTime(self, potential_n, comm):
        '''
        if V[i] >= V_thres[i] and isActivation[i] == 0 
            then t0[i] = t_a 
        ''' 

        mesh = self.parameters["mesh"]

        V_thres = self.parameters["Threshold_Potential"] 
        current_ta = self.parameters["t_a"]
        current_ta_array = interpolate(current_ta, self.Quad).vector().array()

        t_init_array = self.t_init.vector().array() 
        isActive_array = self.isActive.vector().array()

        phi_n_interp = interpolate(potential_n, self.Quad)

        V_n_array = phi_n_interp.vector().array() 
        
        tol_isActive = 1E-1
        comm.Barrier()
        for idx, (vn, isactive, tinit, cur_t) in enumerate(zip(V_n_array, isActive_array, t_init_array, current_ta_array)): 
          
            if ( abs(isactive) <= tol_isActive ) and vn >= V_thres : 
                isActive_array[idx] = 1.0 
                t_init_array[idx] = current_ta_array[idx] 

        comm.Barrier()

	# If homogeneous ---------------------------------------------------
	isHomogenousActivation = self.parameters["HomogenousActivation"]
	#print "isHomogeneous = ", isHomogenousActivation
        if isHomogenousActivation: # activation at 10 ms 
            t_init_array = 0.0*np.ones(len(self.t_init.vector().array()))
	#-------------------------------------------------------------------


        self.t_init.vector()[:] = t_init_array
        self.isActive.vector()[:] = isActive_array

    def get_t_init(self):
        
        Quad = self.Quad

        t_init = Function(Quad)

        t_init_array = t_init.vector().array()
        t_init_array = 9999.0*np.ones(len(t_init_array))

        t_init.vector()[:] = t_init_array

        return t_init 


    def restart_t_init(self):
        self.t_init = self.get_t_init()
        self.isActive = self.get_isActive()



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


    def GetPact(self):
	return self.activeforms.PK1Stress()


    def PK2StressTensor(self):

	F = self.Fmat()
	f0 = self.parameters["fiber"]
	Mij = f0[i]*f0[j]

	H = self.activeforms.structTensor()
	#Pact = self.PK1Stress()
	Sact = self.activeforms.PK2Stress()

	Pact_tensor = Sact*H #as_tensor(Mij, (i,j)) 

	return Pact_tensor


    def fiberstress_wokappa(self):
	
	PK2 = self.PK2StressTensor()
        F = self.Fe()
	f0 = self.parameters["fiber"]
	J = det(F)

	#Tca = (1.0/J)*PK1*F.T
	#Sca = inv(F)*PK1

	return f0[i]*PK2[i,j]*f0[j]

    def fiberstress(self):
	
	PK2 = self.PK2StressTensor()

	fstress = tr(PK2)

	return fstress 

    def sheetstress(self):
	
	PK2 = self.PK2StressTensor()
        F = self.Fe()
	s0 = self.parameters["sheet"]
	J = det(F)

	#Tca = (1.0/J)*PK1*F.T
	#Sca = inv(F)*PK1

	return s0[i]*PK2[i,j]*s0[j]

    def sheetnormalstress(self):
	
	PK2 = self.PK2StressTensor()
        F = self.Fe()
	n0 = self.parameters["sheet-normal"]
	J = det(F)

	#Tca = (1.0/J)*PK1*F.T
	#Sca = inv(F)*PK1

	return n0[i]*PK2[i,j]*n0[j]


    def CalculateFiberNaturalStrain(self, F_, F_ref, e_fiber, VolSeg):
        
        I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = F_*inv(F_ref)
        # Right Cauchy Green
        C = F.T*F
        
        C_fiber = inner(C*e_fiber, e_fiber)
	E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        E_fiber_BiV = [assemble(E_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return E_fiber_BiV, E_fiber 


    def CalculateFiberNaturalStrain_wkappa(self, F_, F_ref, e_fiber, VolSeg):
        
        I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = F_*inv(F_ref)
        # Right Cauchy Green
        C = F.T*F
        H = self.activeforms.structTensor()

        #C_fiber = tr(H*C)
	C_fiber = tr(C*H)
	E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        E_fiber_BiV = [assemble(E_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return E_fiber_BiV, E_fiber 


    def CalculateLagrangianStrain_wkappa(self, F_, F_ref, e_fiber, VolSeg):
        
        I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = F_*inv(F_ref)
        # Left Cauchy Green
        b = F*F.T
        k = self.activeforms.parameters["material params"]["kap"]
	f0 = self.activeforms.parameters["fiber"]
	ef = F*f0
	Mij = as_tensor(ef[i]*ef[j],(i,j))

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

	h = k*b + (1-3*k)*Mij
	E_fiber = inner(h*e_fiber, e_fiber)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        E_fiber_BiV = [assemble(E_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return E_fiber_BiV, E_fiber 


    def CalculateLagrangianStrain_wkappa2(self, F_, F_ref, e_fiber, VolSeg):
        
        I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = F_*inv(F_ref)
        # Right Cauchy Green
        C = F.T*F
        H = self.activeforms.structTensor()
	hi = H*F.T
	h = F*hi 
	E_fiber = inner(h*e_fiber, e_fiber)

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        E_fiber_BiV = [assemble(E_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return E_fiber_BiV, E_fiber 


    def lmbda_ff_ED(self, F_ref, VolSeg):

	F = self.Fmat()*inv(F_ref)
	f0 = self.parameters["fiber"]
	Cmat = F.T*F
	stensor = self.activeforms.structTensor()
	tensor2 = Cmat*stensor.T
	
	lmbda_ff_ED = sqrt(tr(tensor2)) #sqrt(dot(f0, Cmat*f0))

	mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        l_fiber_lv = [assemble(lmbda_ff_ED * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

	return l_fiber_lv, lmbda_ff_ED

    def lmbda_ss_ED(self, F_ref, VolSeg):

	F = self.Fmat()*inv(F_ref)
	s0 = self.parameters["sheet"]
	Cmat = F.T*F
	stensor = self.activeforms.structTensor()
	tensor2 = stensor*Cmat
	lmbda_ss_ED = sqrt(dot(s0, Cmat*s0))

	mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        l_sheet_lv = [assemble(lmbda_ss_ED * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

	return l_sheet_lv, lmbda_ss_ED


    def lmbda_nn_ED(self, F_ref, VolSeg):

	F = self.Fmat()*inv(F_ref)
	n0 = self.parameters["sheet-normal"]
	Cmat = F.T*F
	stensor = self.activeforms.structTensor()
	tensor2 = stensor*Cmat
	lmbda_nn_ED = sqrt(dot(n0, Cmat*n0))

	mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        l_sheetnormal_lv = [assemble(lmbda_nn_ED * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

	return l_sheetnormal_lv, lmbda_nn_ED


    def cauchy(self):
	

	F = self.Fmat()

	Sact = self.PK2StressTensor() 
	J = det(F)

	P =F*Sact	
	sigma = (1.0/J)*P*F.T
	return sigma


    def CalculateFiberCauchyStress(self, VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()

	Sact = self.PK2StressTensor() 
	J = det(F)

	P =F*Sact	
	sigma = (1.0/J)*P*F.T
        H = self.activeforms.structTensor()
	hi = H*F.T
	h = F*hi 
	sigma_fiber = tr(sigma*H)

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        sigma_fiber_LV = [assemble(sigma_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return sigma_fiber 


    def CalculateFiberPK2Stress_1(self, VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()

	Sact = self.PK2StressTensor() 

        H = self.activeforms.structTensor()
	hi = H*F.T
	h = F*hi 
	sigma_fiber = tr(Sact*H)


        return sigma_fiber 


    def CalculateFiberPK2Stress_2(self, VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()

	Sact = self.PK2StressTensor() 

        H = self.activeforms.structTensor()
	f0 = self.parameters["fiber"]
	sigma_fiber = f0[i]*Sact[i,j]*f0[j]


        return sigma_fiber

    def CalculateFiberPK2Stress_3(self, VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()
        k = self.activeforms.parameters["material params"]["kap"]
	Sact = self.PK2StressTensor() 

        H = self.activeforms.structTensor()
	f0 = self.parameters["fiber"]
	sigma_fiber = f0[i]*Sact[i,j]*f0[j]


        return (1-3*k)*sigma_fiber


    def CalculateFiberPK2Stress_4(self, VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()
        k = self.activeforms.parameters["material params"]["kap"]
	Sact = self.PK2StressTensor() 

        H = self.activeforms.structTensor()
	f0 = self.parameters["fiber"]
	sigma_fiber = f0[i]*Sact[i,j]*f0[j]


        return (1-2*k)*sigma_fiber


    def CalculateFiberCauchyStress3(self, VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()

	Sact = self.PK2StressTensor() 
	J = det(F)

	P =F*Sact	
	sigma = (1.0/J)*P*F.T
        H = self.activeforms.structTensor()
	hi = H*F.T
	h = F*hi 
	sigma_fiber = tr(H*sigma)

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        sigma_fiber_LV = [assemble(sigma_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return sigma_fiber 


    def CalculateSheetCauchyStress(self,  VolSeg):
        
        #I = Identity(F_.ufl_domain().geometric_dimension())       

        # Pull back to ED
        F = self.Fmat()

	Sact = self.PK2StressTensor() 
	J = det(F)

	P =F*Sact	
	sigma = (1.0/J)*P*F.T
        H = self.activeforms.structTensor()
	e_fiber = self.parameters["sheet"]
	n  = (F*e_fiber)/ sqrt(inner(F*e_fiber, F*e_fiber))

	sigma_fiber = dot(sigma*n, n)

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        sigma_fiber_LV = [assemble(sigma_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return sigma_fiber 

    def CalculateSheetnormalCauchyStress(self, VolSeg):

        #F = self.Fmat() #F_*inv(F_ref)
        
        #I = Identity(F.ufl_domain().geometric_dimension())       

        # Pull back to ED

        # Right Cauchy Green
        #sigma = self.cauchy()


	F = self.Fmat()

	Sact = self.PK2StressTensor() 
	J = det(F)

	P =F*Sact	
	sigma = (1.0/J)*P*F.T
        H = self.activeforms.structTensor()
	e_fiber = self.parameters["sheet-normal"]
	n  = (F*e_fiber)/ sqrt(inner(F*e_fiber, F*e_fiber))

	sigma_fiber = dot(sigma*n, n)

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        sigma_fiber_LV = [assemble(sigma_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return sigma_fiber 


    def CalculateFiberCauchyStress2(self, VolSeg):

        F = self.Fmat() #F_*inv(F_ref)
        
        I = Identity(F.ufl_domain().geometric_dimension())       

        # Pull back to ED


        # Left Cauchy Green
        b = F*F.T
        k = self.activeforms.parameters["material params"]["kap"]
	f0 = self.parameters["fiber"]
	ef = (F*f0)/sqrt(inner(F*f0, F*f0))
	Mij = as_tensor(ef[i]*ef[j],(i,j))

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

	h = k*b + (1-3*k)*Mij

        # Right Cauchy Green
        sigma = self.cauchy()
	sigma_fiber = tr(sigma*h)

        #C_fiber = tr(H*C)
	#E_fiber = 0.5*(1 - 1/C_fiber)
	#E_fiber = 0.5*(C_fiber-1)

        mesh = self.parameters["mesh"]
        dx = self.parameters["dx"]

        sigma_fiber_LV = [assemble(sigma_fiber * dx(ii), form_compiler_parameters={"representation":"uflacs"}) for ii in VolSeg]

        return sigma_fiber 



	



