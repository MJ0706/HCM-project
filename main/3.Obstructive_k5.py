import sys, pdb
from dolfin import * 
#make sure the following path is correct and all necessary folders and codes are within the same path
sys.path.append("/mnt/home/mojumder/github/HCM")
sys.path.append("/mnt/home/mojumder/github/HCM/src2")

from src2.sim_protocols.run_BiV_ClosedLoop import run_BiV_ClosedLoop as run_BiV_ClosedLoop
from src2.postprocessing.postprocessdata2_3 import postprocessdata as postprocessdata
from finalcode.ExtractStress import ExtractStress as ExtractStress
from finalcode.extract_work import ExtractWork as ExtractWork

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 

IODetails = {"casename" : "case_P3", #Change casename for specific patient
             "directory_me" : './../fine3/', #check the folder location is correct
             "directory_ep" : './../fine3/', 
             "outputfolder" : './with_dispersion/P3/k5/simulation/', #check output folder 
             "folderName" : '',
             "caseID" : 'LVelectromechanics',
             "isLV" : True}


#Change contractility mentioned in manuscript
contRactility = 1.75*95.0e3

#Change dispersion parameter mentioned in manuscript 
dispersion = 0.22 

#Change the passive material parameters as mentioned in manuscript
GuccioneParams = {"ParamsSpecified" : True, 
		 
		 "Passive model": {"Name": "HolzapfelOgden"},
     		 "Passive params": {"a": Constant(0.8*5*50.0E00),
        	    		    "b": Constant(2.0*2.0),
        	    		    "a_f": Constant(0.1*15.0E01),
        	    		    "b_f": Constant(5.893*2.5*1.5),
        	    		    "a_s": Constant(0.492E03),
        	    		    "b_s": Constant(3.393),
        	    		    "a_fs": Constant(0.070E03),
        	    		    "b_fs": Constant(3.929),
				    "kappa" : Constant(dispersion),
				   },
	         "Active model": {"Name": "Time-varying"},

		 #Change active params mentioned in HCM manuscript
     	         "Active params": {"tau" : 35, "t_trans" : 420, "B" : 4.75,  "t0" : 350,  "l0" : 1.55, \
				   "Tmax" : Constant(contRactility), "Ca0" : 4.35, "Ca0max" : 4.35, "lr" : 1.85, "kap":Constant(dispersion)},
                 "HomogenousActivation": True,
		 "deg" : 4, 
                 "Kappa": 1e5,
                 "incompressible" : True,
                 }
#Change the circulatory model parameters as mentioned in HCM manuscript
Circparam = {     "Ees_la": 10.0,
        	  "A_la": 2.67*1.5,
        	  "B_la": 0.019*1.1,
        	  "V0_la": 10.0,
        	  "Tmax_la": 150.0, 
        	  "tau_la": 25.0,
		  "tdelay_la": 140.0, 
		  "Csa": 0.0032*3.0*0.5,
		  "Cad": 0.033*0.15,
    		  "Csv" : 0.28*0.05,
    		  "Vsa0" : 360.0*0.85,
		  "Vad0" : 40.0*4.0, 
    		  "Vsv0" : 4600.0,
    		  "Rav" : 500.0*9.0*7.0,
    		  "Rsv" : 100.0,
    		  "Rsa" : 18000.0*2.5,
		  "Rad": 106000.0*1.5,
    		  "Rmv" : 200.0*5.0,
    		  "V_sv" : 4782.46989842,
		  "V_LV" : 61.6290063534,
		  "V_sa" : 368.256108224,
		  "V_ad" : 144.931051643,
		  "V_LA" : 111.902574594,
		  "Q_sv" : -0.147990026153,
		  "Q_av" : 0.0,
		  "Q_sa" : 0.00251965531624,
		  "Q_ad" : 0.0336677556665,
		  "Q_mv" : 0.0124479080844,
		  "V_LV_target": 114.0,
    		  "stop_iter" : 6
		  }


SimDetails = {    
                  "diaplacementInfo_ref": False,
                  "HeartBeatLength": 1180.0,
                  "dt": 1.0,
                  "writeStep": 2.0,
                  "GiccioneParams" : GuccioneParams, 
                  "nLoadSteps": 15, 
                  "DTI_EP": False, 
                  "DTI_ME": False, 
                  "d_iso": 1.5*0.005, 
                  "d_ani_factor": 4.0, 
                  "ploc": [[1.4, 1.4, -3.0, 2.0, 1]],
                  "pacing_timing": [[4.0, 20.0]],
		  "Isclosed": True,
		  "closedloopparam": Circparam,
		  "Pao": 80,
                  "Ppu": 10,
                  "Pmv": 10, 
                  "Ptr": 8,
                  "Ischemia": False,
             	  "isLV" : True,
                  "topid" : 4,
                  "LVendoid" : 2,
                  "RVendoid" : 0,
                  "epiid" : 1,
		  "abs_tol" : 1e-9,
		  "rel_tol" : 5e-7,
		  "isunloading": True,
		  "unloadparam": {"EDP":14.0,"maxit":10, "restol":1e-3, "drestol":1e-4, "volinc":1, "V_LV_target": 114.0, "LVangle":[70,-70]}
                 }

# Run Simulation
run_BiV_ClosedLoop(IODet=IODetails, SimDet=SimDetails)
# Postprocessing

#Execute these next 4 steps one at a time using single processor

#postprocessdata(IODet=IODetails, SimDet=SimDetails)
#ExtractStress(IODet=IODetails, SimDet=SimDetails, param = 'Ell', deno = 1.0)
#ExtractStress(IODet=IODetails, SimDet=SimDetails, param = 'fstressPk2', deno=1000.0)
#ExtractWork(IODet=IODetails, SimDet=SimDetails, param1 = 'fstressPk2', param2 = 'Eff', param3 = 'f')
#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
