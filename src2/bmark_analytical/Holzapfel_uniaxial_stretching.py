from mshr import *
import math as math
import numpy as np
from matplotlib import pylab as plt
from scipy.optimize import fsolve


def Holzapfel_analytical_uniaxial(lbda, matparams={}):

	default_parameters = {
        	    	      "a": 0.059,
        	    	      "b": 0.023,
        	    	      "a_f": 18.472,
        	    	      "b_f": 16.026,
        	    	      "a_s": 2.481,
        	    	      "b_s": 11.120,
        	    	      "a_fs": 0.216,
        	    	      "b_fs": 11.436, 
			      "kappa" : 0.00,
				};

	default_parameters.update(matparams)

	a = default_parameters["a"]
	b = default_parameters["b"]
	a_f = default_parameters["a_f"]
	b_f = default_parameters["b_f"]
	a_s = default_parameters["a_s"]
	b_s = default_parameters["b_s"]
	a_fs = default_parameters["a_fs"]
	b_fs = default_parameters["b_fs"]
	kappa = default_parameters["kappa"]

	def I1(lbda1, lbda2):
		return lbda1**2 + lbda2**2 + 1/((lbda1*lbda2)**2.0)

	def I4(lbda):
		return lbda**2

	def I4f(I1_f, I4_f):
		return kappa*I1_f +(1- (3*kappa))*I4_f

	lbda1 = lbda

	
	def func(lbda2):

		return a*math.exp(b*(I1(lbda1, lbda2) - 3))*lbda2**2 + \
		       2*a_s*max((I4f(I1(lbda1, lbda2),I4(lbda2)) - 1),0)*math.exp(b_s*(max(I4f(I1(lbda1, lbda2),I4(lbda2)) - 1,0))**2)*lbda2**2 - \
		       a*math.exp(b*(I1(lbda1, lbda2) - 3))*1.0/((lbda1*lbda2)**2)


	lbda2 = fsolve(func, 1)		
	lbda3 = 1.0/(lbda2*lbda1)
	print "lbda 1 = ", lbda1,  "lbda 2 = ", lbda2, "lbda 3 = ", lbda3

	p = a*math.exp(b*(I1(lbda1, lbda2) - 3))*1.0/((lbda1*lbda2)**2)

	#sigma11 = a*math.exp(b*(I1(lbda1, lbda2) - 3))*lbda1**2 + 2*a_f*(max(I4f(I1(lbda1, lbda2),I4(lbda1)) - 1,0))*math.exp(b_f*(max(I4f(I1(lbda1, lbda2),I4(lbda1)) - 1,0))**2.0)*lbda1**2 - p

	sigma11 = a*math.exp(b*(I1(lbda1, lbda2) - 3))*lbda1**2 + 2*(1-3*kappa)*a_f*(max(I4f(I1(lbda1, lbda2),I4(lbda1)) - 1,0))*math.exp(b_f*(max(I4f(I1(lbda1, lbda2),I4(lbda1)) - 1,0))**2.0)*lbda1**2 - p + kappa*a_f*(max(I4f(I1(lbda1, lbda2),I4(lbda1)) - 1,0))*math.exp(b_f*(max(I4f(I1(lbda1, lbda2),I4(lbda1)) - 1,0))**2.0)*lbda1**2

	PK11 = sigma11/lbda1

	return PK11
