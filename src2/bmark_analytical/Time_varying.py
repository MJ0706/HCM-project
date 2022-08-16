#from mshr import *
#import fenicstools as ft
import math as math
#import numpy as np
#from forms2 import Forms 
#from scipy.optimize import fsolve


def Time_varying(lbda,t_a, matparams={}):

        default_parameters =   {"tau" : 30,
				"t_trans" : 320,
				"B" : 4.75, #4.75 no strech dependency seen, 
				"t0" : 275, 
                		"deg" : 4,
                		"l0" : 1.6, #1.58
                		"lr": 1.85,
        			"Tmax" : 135.7e3,
				"Ca0" : 4.35,
				"Ca0max" : 4.35
				}

	default_parameters.update(matparams)

	B = default_parameters["B"]
	l0 = default_parameters["l0"]
	Ca0max = default_parameters["Ca0max"]
	t0 = default_parameters["t0"]
	t_trans = default_parameters["t_trans"]
	tr = default_parameters["tau"]
	Ca0 = default_parameters["Ca0"]
	lr = default_parameters["lr"]

	try:
		Tmax = float(default_parameters["Tmax"])
	except TypeError:
		Tmax = default_parameters["Tmax"].val



	ls = lbda*lr
	ls_l0 = max(ls - l0, 0.0002)
	denom = math.sqrt(math.exp(B*(ls_l0)) - 1)

	ECa = Ca0max/denom

	if(t_a < t_trans):
		Ct = 0.5*(1 - math.cos(math.pi*t_a/t0))
	else:
		Ct = 0.5*(1 - math.cos(math.pi*t_trans/t0))*math.exp(-1.0*(t_a - t_trans)/tr)


	Pact = (Tmax*Ct*Ca0**2.0)/(Ca0**2.0 + ECa**2.0)

	return Pact

