from dolfin import *


class HolzapfelOgden(object):

    	r"""
    	Orthotropic model by Holzapfel and Ogden

    	.. math::

    	   \mathcal{W}(I_1, I_{4f_0})
    	   = \frac{a}{2 b} \left( e^{ b (I_1 - 3)}  -1 \right)
    	   + \frac{a_f}{2 b_f} \left( e^{ b_f (I_{4f_0} - 1)_+^2} -1 \right)
    	   + \frac{a_s}{2 b_s} \left( e^{ b_s (I_{4s_0} - 1)_+^2} -1 \right)
    	   + \frac{a_fs}{2 b_fs} \left( e^{ b_fs I_{8fs}^2} -1 \right)
    	where

    	.. math::

    	   (\cdot)_+ = \max\{x,0\}


    	.. rubric:: Reference

    	[1] Holzapfel, Gerhard A., and Ray W. Ogden.
    	"Constitutive modelling of passive myocardium:
    	a structurally based framework for material characterization.
    	"Philosophical Transactions of the Royal Society of London A:
    	Mathematical, Physical and Engineering Sciences 367.1902 (2009): 3445-3475.

    	"""
    	name = "holzapfel_ogden"

    	def __init__(self, params):        

		self.parameters = self.default_parameters()
         	self.parameters.update(params)

    	def default_parameters(self):
        	"""
        	Default matereial parameter for the Holzapfel Ogden model

        	Taken from Table 1 row 3 of [1]
        	"""

        	return {"material params": {
        	    			      "a": 0.059,
        	    			      "b": 0.023,
        	    			      "a_f": 18.472,
        	    			      "b_f": 16.026,
        	    			      "a_s": 2.481,
        	    			      "b_s": 11.120,
        	    			      "a_fs": 0.216,
        	    			      "b_fs": 11.436}}
	def Getmatparam(self):

		return self.parameters["material params"]

	def I1(self):

		F = self.parameters["F"]
		Cmat = F.T*F
		I1 = tr(Cmat)

		return I1

	def I4f(self):

		f0 = self.parameters["fiber"]
		F = self.parameters["F"]
		Cmat = F.T*F
		I4f = inner(Cmat*f0, f0)

		return I4f

	def I4s(self):

		s0 = self.parameters["sheet"]
		F = self.parameters["F"]
		Cmat = F.T*F
		I4s = inner(Cmat*s0, s0)

		return I4s

	def I8(self):

		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		F = self.parameters["F"]
		Cmat = F.T*F
		I8 = inner(Cmat*f0, s0)

		return I8





    	def W_1(self, diff=0, *args, **kwargs):
        	r"""
        	Isotropic contribution.

        	If `diff = 0`, return

        	.. math::

        	   \frac{a}{2 b} \left( e^{ b (I_1 - 3)}  -1 \right)

        	If `diff = 1`, return

        	.. math::

        	   \frac{a}{b} e^{ b (I_1 - 3)}

        	If `diff = 2`, return

        	.. math::

        	   \frac{a b}{2}  e^{ b (I_1 - 3)}

        	"""

        	a = self.parameters["material params"]["a"]
        	b = self.parameters["material params"]["b"]

        	if diff == 0:
        	    try:
        	        if float(a) > dolfin.DOLFIN_EPS:
        	            if float(b) > dolfin.DOLFIN_EPS:
        	                return a / (2.0 * b) * (dolfin.exp(b * (self.I1() - 3)) - 1.0)
        	            else:
        	                return a / 2.0 * (self.I1() - 3)
        	        else:
        	            return 0.0
        	    except Exception:
        	        return a / (2.0 * b) * (dolfin.exp(b * (self.I1() - 3)) - 1)
        	elif diff == 1:
        	    return a / 2.0 * dolfin.exp(b * (self.I1() - 3))
        	elif diff == 2:
        	    return a * b / 2.0 * dolfin.exp(b * (self.I1() - 3))

    	def W_4s(self, diff=0, use_heaviside=False, *args, **kwargs):
	        r"""
        	Anisotropic contribution.

        	If `diff = 0`, return

        	.. math::

        	   \frac{a_f}{2 b_f} \left( e^{ b_f (I_{4f_0} - 1)_+^2} -1 \right)

        	If `diff = 1`, return

        	.. math::

        	   a_f (I_{4f_0} - 1)_+ e^{ b_f (I_{4f_0} - 1)^2}

        	If `diff = 2`, return

        	.. math::

        	   a_f h(I_{4f_0} - 1) (1 + 2b(I_{4f_0} - 1))
        	   e^{ b_f (I_{4f_0} - 1)_+^2}

        	where

        	.. math::

        	   h(x) = \frac{\mathrm{d}}{\mathrm{d}x} \max\{x,0\}

        	is the Heaviside function.

        	"""
        	a = self.parameters["material params"]["a_s"]
        	b = self.parameters["material params"]["b_s"]

        	if self.I4s() == 0:
        	    return 0

        	if diff == 0:
        	    try:
        	        if float(a) > dolfin.DOLFIN_EPS:
        	            if float(b) > dolfin.DOLFIN_EPS:
        	                return (
        	                    a / (2.0 * b) * (dolfin.exp(b * self.subplus(self.I4s() - 1) ** 2) - 1.0)
        	                )
        	            else:
        	                return a / 2.0 * self.subplus(self.I4s() - 1) ** 2
        	        else:
        	            return 0.0
        	    except Exception:
        	        # Probably failed to convert a and b to float
        	        return a / (2.0 * b) * (dolfin.exp(b * self.subplus(self.I4s() - 1) ** 2) - 1.0)

        	elif diff == 1:
        	    return a * self.subplus(self.I4s() - 1) * dolfin.exp(b * pow(self.I4s() - 1, 2))
        	elif diff == 2:
        	    return (
        	        a
        	        * heaviside(self.I4s() - 1)
        	        * (1 + 2.0 * b * pow(self.I4s() - 1, 2))
        	        * dolfin.exp(b * pow(self.I4s() - 1, 2))
        	    )

    	def W_8(self, *args, **kwargs):
        	"""
        	Cross fiber-sheet contribution.
        	"""
        	a = self.parameters["material params"]["a_fs"]
        	b = self.parameters["material params"]["b_fs"]

        	try:
        	    if float(a) > dolfin.DOLFIN_EPS:
        	        if float(b) > dolfin.DOLFIN_EPS:
        	            return a / (2.0 * b) * (dolfin.exp(b * self.I8() ** 2) - 1.0)
        	        else:
        	            return a / 2.0 * self.I8() ** 2
        	    else:
        	        return 0.0
        	except Exception:
        	    return a / (2.0 * b) * (dolfin.exp(b * self.I8() ** 2) - 1.0)



    	def W_4f(self, diff=0, use_heaviside=False, *args, **kwargs):
	        r"""
        	Anisotropic contribution.

        	If `diff = 0`, return

        	.. math::

        	   \frac{a_f}{2 b_f} \left( e^{ b_f (I_{4f_0} - 1)_+^2} -1 \right)

        	If `diff = 1`, return

        	.. math::

        	   a_f (I_{4f_0} - 1)_+ e^{ b_f (I_{4f_0} - 1)^2}

        	If `diff = 2`, return

        	.. math::

        	   a_f h(I_{4f_0} - 1) (1 + 2b(I_{4f_0} - 1))
        	   e^{ b_f (I_{4f_0} - 1)_+^2}

        	where

        	.. math::

        	   h(x) = \frac{\mathrm{d}}{\mathrm{d}x} \max\{x,0\}

        	is the Heaviside function.

        	"""
        	a = self.parameters["material params"]["a_f"]
        	b = self.parameters["material params"]["b_f"]

        	if self.I4f() == 0:
        	    return 0

        	if diff == 0:
        	    try:
        	        if float(a) > dolfin.DOLFIN_EPS:
        	            if float(b) > dolfin.DOLFIN_EPS:
        	                return (
        	                    a / (2.0 * b) * (dolfin.exp(b * self.subplus(self.I4f() - 1) ** 2) - 1.0)
        	                )
        	            else:
        	                return a / 2.0 * self.subplus(self.I4f() - 1) ** 2
        	        else:
        	            return 0.0
        	    except Exception:
        	        # Probably failed to convert a and b to float
        	        return a / (2.0 * b) * (dolfin.exp(b * self.subplus(self.I4f() - 1) ** 2) - 1.0)

        	elif diff == 1:
        	    return a * self.subplus(self.I4f() - 1) * dolfin.exp(b * pow(self.I4f() - 1, 2))
        	elif diff == 2:
        	    return (
        	        a
        	        * heaviside(self.I4f() - 1)
        	        * (1 + 2.0 * b * pow(self.I4f() - 1, 2))
        	        * dolfin.exp(b * pow(self.I4f() - 1, 2))
        	    )


    	def PassiveMatSEF(self):

        	r"""
        	Strain-energy density function.

        	.. math::

        	   \mathcal{W} = \mathcal{W}_1 + \mathcal{W}_{4f}
        	   + \mathcal{W}_{\mathrm{active}}

        	where

        	.. math::

        	   \mathcal{W}_{\mathrm{active}} =
        	   \begin{cases}
        	     0 & \text{if acitve strain} \\
        	     \gamma I_{4f} & \text{if active stress}
        	   \end{cases}


        	:param F: Deformation gradient
        	:type F: :py:class:`dolfin.Function`

        	"""

        	#dim = get_dimesion(F)
        	W1 = self.W_1(diff=0)
        	W4f = self.W_4f(diff=0)
        	W4s = self.W_4s(diff=0)
        	W8fs = self.W_8(diff=0)

        	W = W1 + W4f + W4s + W8fs

        	return W

	def PK1(self):

		PK1 = self.PK1_W1(diff=0) + \
		      self.PK1_W4f(diff=0) + \
		      self.PK1_W4s(diff=0) + \
		      self.PK1_W8()

		return PK1

	def PK1_W8(self, *args, **kwargs):
        	"""
        	Cross fiber-sheet contribution.
        	"""
        	a = self.parameters["material params"]["a_fs"]
        	b = self.parameters["material params"]["b_fs"]

		u = self.parameters["displacement_variable"]
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		F = dolfin.variable(F)
		Cmat = F.T*F
		J = det(F)
		I8 = inner(Cmat*f0, s0)

        	try:
        	    if float(a) > dolfin.DOLFIN_EPS:
        	        if float(b) > dolfin.DOLFIN_EPS:
        	            W8 = a / (2.0 * b) * (dolfin.exp(b * I8 ** 2) - 1.0)
        	        else:
        	            W8 = a / 2.0 * I8 ** 2
        	    else:
        	        W8 = 0.0
        	except Exception:
        	    W8 = a / (2.0 * b) * (dolfin.exp(b * I8 ** 2) - 1.0)


		PK1 = dolfin.diff(W8,F)               

                return PK1


	def PK1_W4s(self, diff=0, *args, **kwargs):

        	a = self.parameters["material params"]["a_s"]
        	b = self.parameters["material params"]["b_s"]

		u = self.parameters["displacement_variable"]
		s0 = self.parameters["sheet"]
		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		F = dolfin.variable(F)
		Cmat = F.T*F
		J = det(F)
		I4s = inner(Cmat*s0, s0)

        	if I4s == 0:
        	    W4s =  0

        	if diff == 0:
        	    try:
        	        if float(a) > dolfin.DOLFIN_EPS:
        	            if float(b) > dolfin.DOLFIN_EPS:
        	                W4s =  (
        	                    a / (2.0 * b) * (dolfin.exp(b * self.subplus(I4s - 1) ** 2) - 1.0)
        	                )
        	            else:
        	                W4s = a / 2.0 * self.subplus(I4s - 1) ** 2
        	        else:
        	            return 0.0
        	    except Exception:
        	        # Probably failed to convert a and b to float
        	        W4s =  a / (2.0 * b) * (dolfin.exp(b * self.subplus(I4s - 1) ** 2) - 1.0)

        	elif diff == 1:
        	    W4s = a * self.subplus(I4s - 1) * dolfin.exp(b * pow(I4s - 1, 2))
        	elif diff == 2:
        	    W4s = (
        	        a
        	        * heaviside(I4s - 1)
        	        * (1 + 2.0 * b * pow(I4s - 1, 2))
        	        * dolfin.exp(b * pow(I4s - 1, 2))
        	    )

                PK1 = dolfin.diff(W4s,F)               


                return PK1


	def PK1_W4f(self, diff=0, *args, **kwargs):

        	a = self.parameters["material params"]["a_f"]
        	b = self.parameters["material params"]["b_f"]

		u = self.parameters["displacement_variable"]
		f0 = self.parameters["fiber"]
		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		F = dolfin.variable(F)
		Cmat = F.T*F
		J = det(F)
		I4f = inner(Cmat*f0, f0)

        	if I4f == 0:
        	    W4f =  0

        	if diff == 0:
        	    try:
        	        if float(a) > dolfin.DOLFIN_EPS:
        	            if float(b) > dolfin.DOLFIN_EPS:
        	                W4f =  (
        	                    a / (2.0 * b) * (dolfin.exp(b * self.subplus(I4f - 1) ** 2) - 1.0)
        	                )
        	            else:
        	                W4f = a / 2.0 * self.subplus(I4f - 1) ** 2
        	        else:
        	            return 0.0
        	    except Exception:
        	        # Probably failed to convert a and b to float
        	        W4f =  a / (2.0 * b) * (dolfin.exp(b * self.subplus(I4f - 1) ** 2) - 1.0)

        	elif diff == 1:
        	    W4f = a * self.subplus(I4f - 1) * dolfin.exp(b * pow(I4f - 1, 2))
        	elif diff == 2:
        	    W4f = (
        	        a
        	        * heaviside(I4f - 1)
        	        * (1 + 2.0 * b * pow(I4f - 1, 2))
        	        * dolfin.exp(b * pow(I4f - 1, 2))
        	    )

                PK1 = dolfin.diff(W4f,F)               


                return PK1


	def PK1_W1(self, diff=0, *args, **kwargs):

        	a = self.parameters["material params"]["a"]
        	b = self.parameters["material params"]["b"]

		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
	    	I = Identity(d)
	    	F = I + grad(u)
		F = dolfin.variable(F)
		Cmat = F.T*F
		J = det(F)
		I1 = tr(Cmat)

        	if diff == 0:
        	    try:
        	        if float(a) > dolfin.DOLFIN_EPS:
        	            if float(b) > dolfin.DOLFIN_EPS:
        	                W1 = a / (2.0 * b) * (dolfin.exp(b * (I1 - 3)) - 1.0)
        	            else:
        	                W1 =  a / 2.0 * (I1 - 3)
        	        else:
        	            return 0.0
        	    except Exception:
        	        W1 =  a / (2.0 * b) * (dolfin.exp(b * (I1 - 3)) - 1)
        	elif diff == 1:
        	    W1 =  a / 2.0 * dolfin.exp(b * (I1 - 3))
        	elif diff == 2:
        	    W1 =  a * b / 2.0 * dolfin.exp(b * (I1 - 3))


                PK1 = dolfin.diff(W1,F)               

                return PK1



	def subplus(self, x):
		r"""
		Ramp function
		.. math::
		   \max\{x,0\}
		"""
		
		return dolfin.conditional(dolfin.ge(x, 0.0), x, 0.0)

	def heaviside(self, x):
	    	r"""
	    	Heaviside function
	    	.. math::
	    	   \frac{\mathrm{d}}{\mathrm{d}x} \max\{x,0\}
	    	"""
	
	    	return dolfin.conditional(dolfin.ge(x, 0.0), 1.0, 0.0)
	
	def get_dimesion(self, u):
	
	    	# TODO : Check argument
	    	try:
	    	    if DOLFIN_VERSION_MAJOR > 1.6:
	    	        from ufl.domain import find_geometric_dimension
	
	    	        dim = find_geometric_dimension(u)
	    	    else:
	    	        dim = u.geometric_dimension()
	
	    	except Exception as ex:
	
	    	    try:
	    	        dim = len(u)
	    	    except Exception as ex2:
	    	        logger.warning(ex)
	    	        logger.warning(ex2)
	    	        # Assume dimension is 3
	    	        logger.warning("Assume dimension is 3")
	    	        dim = 3
	
	    	return dim

