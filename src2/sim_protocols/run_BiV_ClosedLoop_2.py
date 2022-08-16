import sys, shutil, pdb, math
import os as os
import numpy as np
from mpi4py import MPI as pyMPI 

from dolfin import * 
from fenicstools import *

import vtk_py
import vtk

from ..utils.oops_objects_MRC2 import printout
from ..utils.oops_objects_MRC2 import biventricle_mesh as biv_mechanics_mesh
from ..utils.oops_objects_MRC2 import lv_mesh as lv_mechanics_mesh

from ..utils.oops_objects_MRC2 import State_Variables
from ..utils.oops_objects_MRC2 import update_mesh
from ..utils.oops_objects_MRC2 import exportfiles

from ..utils.mesh_scale_create_fiberFiles import create_EDFibers

from ..ep.EPmodel import EPmodel
from ..mechanics.MEmodel3 import MEmodel

def run_BiV_ClosedLoop(IODet, SimDet): 

    deg = 4
    flags = ["-O3", "-ffast-math", "-march=native"]
    parameters["form_compiler"]["representation"]="quadrature"
    parameters["form_compiler"]["quadrature_degree"]=deg

    casename = IODet["casename"]
    directory_me = IODet["directory_me"]
    directory_ep = IODet["directory_ep"]
    outputfolder = IODet["outputfolder"]
    folderName = IODet["folderName"] + IODet["caseID"] + '/' 
    isLV = SimDet["isLV"]
    
    delTat = SimDet["dt"]

    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    # Read EP data from HDF5 Files
    mesh_ep = Mesh() 
    comm_common = mesh_ep.mpi_comm()

    meshfilename_ep = directory_ep + casename + "_refine.hdf5"
    f = HDF5File(comm_common, meshfilename_ep, 'r')
    f.read(mesh_ep, casename, False)

    File(outputfolder + folderName + "mesh_ep.pvd") << mesh_ep 

    facetboundaries_ep = MeshFunction("size_t", mesh_ep, 2)
    f.read(facetboundaries_ep, casename+"/"+"facetboundaries")

    matid_ep = CellFunction('size_t', mesh_ep)
    AHAid_ep = CellFunction('size_t', mesh_ep)
    if(f.has_dataset(casename+"/"+"matid")):
    	f.read(matid_ep, casename+"/"+"matid")
    else:
	matid_ep.set_all(0)

    if(f.has_dataset(casename+"/"+"AHAid")):
    	f.read(AHAid_ep, casename+"/"+"AHAid")
    else:
	AHAid_ep.set_all(0)

    deg_ep = 4

    Quadelem_ep = FiniteElement("Quadrature", mesh_ep.ufl_cell(), degree=deg_ep, quad_scheme="default")
    Quadelem_ep._quad_scheme = 'default'
    Quad_ep = FunctionSpace(mesh_ep, Quadelem_ep)

    VQuadelem_ep = VectorElement("Quadrature", 
                                  mesh_ep.ufl_cell(), 
                                  degree=deg_ep, 
                                  quad_scheme="default")
    VQuadelem_ep._quad_scheme = 'default'
    
    fiberFS_ep = FunctionSpace(mesh_ep, VQuadelem_ep)

    f0_ep = Function(fiberFS_ep)
    s0_ep = Function(fiberFS_ep)
    n0_ep = Function(fiberFS_ep)

    
    if SimDet["DTI_EP"] is True : 
        f.read(f0_ep, casename+"/"+"eF_DTI")
        f.read(s0_ep, casename+"/"+"eS_DTI")
        f.read(n0_ep, casename+"/"+"eN_DTI")
    else: 
        f.read(f0_ep, casename+"/"+"eF")
        f.read(s0_ep, casename+"/"+"eS")
        f.read(n0_ep, casename+"/"+"eN") 
    
    f.close()

    comm_ep = mesh_ep.mpi_comm()


    # Define state variables
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    state_obj = State_Variables(comm_ep, SimDet) 
    state_obj.dt.dt = delTat
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 

    EPparams = {"EPmesh": mesh_ep,\
                "deg": 4,\
		"matid": matid_ep,\
		"facetboundaries": facetboundaries_ep,\
		"f0": f0_ep,\
		"s0": s0_ep,\
		"n0": n0_ep,\
		"state_obj": state_obj,\
                "d_iso": SimDet["d_iso"],\
                "d_ani_factor": SimDet["d_ani_factor"],\
	        "AHAid": AHAid_ep,\
                "matid": matid_ep}

    if("ploc" in SimDet.keys()):
	EPparams.update({"ploc": SimDet["ploc"]})
    if("Ischemia" in SimDet.keys()):
	EPparams.update({"Ischemia": SimDet["Ischemia"]})
    if("pacing_timing" in SimDet.keys()):
	EPparams.update({"pacing_timing": SimDet["pacing_timing"]})

    # Define EP model and solver
    EPmodel_ = EPmodel(EPparams)
    EpiBCid_ep = EPmodel_.MarkStimulus()
	
    solver_FHN = EPmodel_.Solver()
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    # Mechanics Mesh 

    mesh_me = Mesh() 
    mesh_me_params = {"directory" : directory_me, 
                   "casename" : casename, 
                   "fibre_quad_degree" : 4, 
                   "outputfolder" : outputfolder,
                   "foldername" : folderName,
		   "state_obj": state_obj,
                   "common_communicator": comm_common,
                   "MEmesh": mesh_me,
		   "isLV": isLV}

    MEmodel_ = MEmodel(mesh_me_params, SimDet)
    solver_elas = MEmodel_.Solver()
    comm_me = MEmodel_.mesh_me.mpi_comm()

    # Set up export class
    export = exportfiles(comm_me, comm_ep, IODet, SimDet)
    export.exportVTKobj("facetboundaries_ep.pvd", facetboundaries_ep) 
    export.exportVTKobj("EpiBCid_ep.pvd", EpiBCid_ep) 

    F_ED = Function(MEmodel_.TF)

    if("AHA_segments" in SimDet.keys()):
    	AHA_segments = SimDet["AHA_segments"]
    else:
	AHA_segments = [0]

    # Get Unloaded volumes
    V_LV_unload = MEmodel_.GetLVV()
    V_RV_unload = MEmodel_.GetRVV()

    nloadstep = SimDet["nLoadSteps"]

    # Unloading LV to get new reference geometry
    isunloading = False
    if("isunloading" in SimDet.keys()):
    	if (SimDet["isunloading"] is True)  : 
		isunloading = True

    if(isunloading):

        printout("Start UnLoading", comm_me)
	#V_LV_target = MEmodel_.GetLVV()

	if("unloadparam" in SimDet.keys()):
		unloadparam = SimDet["unloadparam"]
	else:
		unloadparam = {};
			
    	nloadstep_, volinc_, V_LV_target = MEmodel_.unloading(unloadparam)
	printout("Target EDV = " + str(V_LV_target), comm_me)

        printout("Finish UnLoading and Reloading", comm_me)

 	export.writePV(MEmodel_, 0);
    	export.hdf.write(MEmodel_.mesh_me, "ME/mesh")
    	export.hdf.write(EPmodel_.mesh_ep, "EP/mesh")

	MEmodel_.LVCavityvol.vol = MEmodel_.GetLVV()
	V_LV_unload = MEmodel_.GetLVV()
	File(outputfolder + folderName +"/loadingMesh/mesh_me_"+ ".pvd") << MEmodel_.mesh_me
	File(outputfolder + folderName +"/loadingMesh/mesh_ep_"+ ".pvd") << EPmodel_.mesh_ep
	for it in np.arange(0,nloadstep):

        	#MEmodel_.LVCavityvol.vol += (V_LV_target - V_LV_unload)/nloadstep
		if("V_LV" in SimDet["closedloopparam"].keys()):
    	    		V_LV_target = SimDet["closedloopparam"]["V_LV"]
    	    		MEmodel_.LVCavityvol.vol += (V_LV_target - V_LV_unload)/nloadstep
		else:
    	    		MEmodel_.LVCavityvol.vol += 2.0
    		solver_elas.solvenonlinear() 
		printout("Pressure = " +  str(MEmodel_.GetLVP()*0.0075) +  " Vol = " + str(MEmodel_.GetLVV()), comm_me)

    		export.writePV(MEmodel_, 0);
    		export.hdf.write(MEmodel_.GetDisplacement(), "ME/u_loading", it)

    		F_ED.vector()[:] = project(MEmodel_.GetFmat(), MEmodel_.TF, solver_type='mumps').vector().array()[:] 

    # No unloading
    else:
    	export.writePV(MEmodel_, 0);
    	export.hdf.write(MEmodel_.mesh_me, "ME/mesh")
    	export.hdf.write(EPmodel_.mesh_ep, "EP/mesh")
	File(outputfolder + folderName +"/loadingMesh/mesh_me_"+ ".pvd") << MEmodel_.mesh_me
	File(outputfolder + folderName +"/loadingMesh/mesh_ep_"+ ".pvd") << EPmodel_.mesh_ep
	Fiberfilenameep = outputfolder + folderName +"/fiberdir/fiber_ep_"
	Fiberfilenameme = outputfolder + folderName +"/fiberdir/fiber_me_"
	#MEmodel_.GetFiberDirection(f0_ep, Fiberfilenameep)
	MEmodel_.GetFiberDirection(MEmodel_.f0_me, Fiberfilenameme)
    	for lmbda_value in range(0, nloadstep):

    	    if("V_LV" in SimDet["closedloopparam"].keys()):
    	    	V_LV_target = SimDet["closedloopparam"]["V_LV"]
    	    	MEmodel_.LVCavityvol.vol += (V_LV_target - V_LV_unload)/nloadstep
    	    else:
    	    	MEmodel_.LVCavityvol.vol += 2.0

    	    if("V_RV" in SimDet["closedloopparam"].keys()):
    	    	V_RV_target = SimDet["closedloopparam"]["V_RV"]
    	    	MEmodel_.RVCavityvol.vol += (V_RV_target - V_RV_unload)/nloadstep
    	    else:
    	    	MEmodel_.RVCavityvol.vol += 2.0


    	    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
    	    try:
    	    	solver_elas.solvenonlinear() 
    	    except:
    		export.hdf.close()
    	    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

	    
	    if(isLV):
    	    	printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.LVCavityvol.vol), comm_me)
    	    	export.printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.LVCavityvol.vol))
	    else:
    	    	printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.LVCavityvol.vol) + \
								      " RVV = " + str(MEmodel_.RVCavityvol.vol), comm_me)
    	    	export.printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.LVCavityvol.vol)+ \
									     " RVV = " + str(MEmodel_.LVCavityvol.vol))

    	                        	
    	    export.writePV(MEmodel_, 0);
    	    export.hdf.write(MEmodel_.GetDisplacement(), "ME/u_loading", lmbda_value)

    	    F_ED.vector()[:] = project(MEmodel_.GetFmat(), MEmodel_.TF, solver_type='mumps').vector().array()[:] 

    if("isunloadingonly" in SimDet.keys()):
    	if (SimDet["isunloadingonly"] is True)  : 
    		export.hdf.close()
		exit()

 

    MEmodel_.F_ED.vector()[:] = F_ED.vector().array()[:]
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
    # Declare communicator based on mpi4py 
    comm_me_ = comm_me.tompi4py()
    eCC, eRR, eLL, deformedMesh, deformedBoundary = MEmodel_.GetDeformedBasis({})


    fStrain = MEmodel_.GetFiberstrain(F_ED)
    fStrain_uL = MEmodel_.GetFiberstrainUL()
    #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

    # Closed-loop phase
    stop_iter = SimDet["closedloopparam"]["stop_iter"]

    # Systemic circulation
    Csa = SimDet["closedloopparam"]["Csa"]
    Cad = SimDet["closedloopparam"]["Cad"]
    Csv = SimDet["closedloopparam"]["Csv"]
    Vsa0 = SimDet["closedloopparam"]["Vsa0"]
    Vad0 = SimDet["closedloopparam"]["Vad0"] 
    Vsv0 = SimDet["closedloopparam"]["Vsv0"]
    Rav =  SimDet["closedloopparam"]["Rav"]
    Rad = SimDet["closedloopparam"]["Rad"] 
    Rsv = SimDet["closedloopparam"]["Rsv"]
    Rsa = SimDet["closedloopparam"]["Rsa"]
    Rmv = SimDet["closedloopparam"]["Rmv"]
    V_sv = SimDet["closedloopparam"]["V_sv"]
    V_sa = SimDet["closedloopparam"]["V_sa"]
    V_ad = SimDet["closedloopparam"]["V_ad"]
    V_LA = SimDet["closedloopparam"]["V_LA"]

    # Pulmonary circulation
    if(not isLV):
    	Cpa = SimDet["closedloopparam"]["Cpa"]
    	Cpv = SimDet["closedloopparam"]["Cpv"]
    	Vpa0 = SimDet["closedloopparam"]["Vpa0"]
    	Vpv0 = SimDet["closedloopparam"]["Vpv0"]
    	Rpv = SimDet["closedloopparam"]["Rpv"]
    	Rtv =  SimDet["closedloopparam"]["Rtv"]
    	Rpa = SimDet["closedloopparam"]["Rpa"]
    	Rpvv = SimDet["closedloopparam"]["Rpvv"]
    	V_pv = SimDet["closedloopparam"]["V_pv"]
    	V_pa = SimDet["closedloopparam"]["V_pa"]
    	V_RA = SimDet["closedloopparam"]["V_RA"]


    isrestart = 0
    prev_cycle = 0
    cnt = 0

    Qav = 0 
    Qmv = 0 
    Qsa = 0 
    Qad = 0
    Qsv = 0 
    Qtv = 0 
    Qpa = 0 
    Qpv = 0 
    Qpvv = 0 
    Qlvad = 0
    Qlara = 0

    if("Q_av" in SimDet["closedloopparam"].keys()):
    	Qav = SimDet["closedloopparam"]["Q_av"] 
    if("Q_mv" in SimDet["closedloopparam"].keys()):
    	Qmv = SimDet["closedloopparam"]["Q_mv"] 
    if("Q_sa" in SimDet["closedloopparam"].keys()):
    	Qsa = SimDet["closedloopparam"]["Q_sa"] 
    if("Q_ad" in SimDet["closedloopparam"].keys()):
    	Qad = SimDet["closedloopparam"]["Q_ad"] 
    if("Q_sv" in SimDet["closedloopparam"].keys()):
    	Qsv = SimDet["closedloopparam"]["Q_sv"] 
    if("Q_tv" in SimDet["closedloopparam"].keys()):
    	Qtv = SimDet["closedloopparam"]["Q_tv"] 
    if("Q_pa" in SimDet["closedloopparam"].keys()):
    	Qpa = SimDet["closedloopparam"]["Q_pa"] 
    if("Q_pv" in SimDet["closedloopparam"].keys()):
    	Qpv = SimDet["closedloopparam"]["Q_pv"] 
    if("Q_pvv" in SimDet["closedloopparam"].keys()):
    	Qpvv = SimDet["closedloopparam"]["Q_pvv"] 
    if("Q_lvad" in SimDet["closedloopparam"].keys()):
    	Qlvad = SimDet["closedloopparam"]["Q_lvad"] 
    if("Q_lara" in SimDet["closedloopparam"].keys()):
    	Qlara = SimDet["closedloopparam"]["Q_lara"] 



    #Parameters for LVAD #############################################
    LVADrpm = 0
    LVADscale = 0
    if("Q_lvad_rpm" in SimDet["closedloopparam"].keys()):
    	LVADrpm = SimDet["closedloopparam"]["Q_lvad_rpm"] 
    if("Q_lvad_scale" in SimDet["closedloopparam"].keys()):
    	LVADscale = SimDet["closedloopparam"]["Q_lvad_scale"] 

    Qlad = 0
    Qlcx = 0

    #Parameters for Shunt #############################################
    Shuntscale =  0.0
    Rsh = 1e9
    if("Shunt_scale" in SimDet["closedloopparam"].keys()):
    	Shuntscale = SimDet["closedloopparam"]["Shunt_scale"] 
    if("Rsh" in SimDet["closedloopparam"].keys()):
    	Rsh = SimDet["closedloopparam"]["Rsh"]


    potential_me = Function(FunctionSpace(MEmodel_.mesh_me,'CG',1))
    writecnt = 0

    File(outputfolder + folderName +"/loadingMesh/mesh_me_2"+ ".pvd") << MEmodel_.mesh_me
    
    while(1):

        if(state_obj.cycle > stop_iter):
            	break;


	#Time varying elastance function fro LA and RA ##################
	def et(t, Tmax, tau):
		if (t <= 1.5*Tmax):
			out = 0.5*(math.sin((math.pi/Tmax)*t - math.pi/2) + 1);
      		else:
			out = 0.5*math.exp((-t + (1.5*Tmax))/tau);

		return out        
	#################################################################
	P_LV = MEmodel_.GetLVP()
        V_LV = MEmodel_.GetLVV()

	if(not isLV):
		P_RV = MEmodel_.GetRVP()
        	V_RV = MEmodel_.GetRVV()

       
        state_obj.tstep = state_obj.tstep + state_obj.dt.dt
        state_obj.cycle = math.floor(state_obj.tstep/state_obj.BCL)
        state_obj.t = state_obj.tstep - state_obj.cycle*state_obj.BCL

	# Update deformation at F_ED
	if(state_obj.cycle > prev_cycle):
        	F_ED.vector()[:] = project(MEmodel_.GetFmat(), MEmodel_.TF).vector().array()[:] 

	prev_cycle = state_obj.cycle

        MEmodel_.t_a.vector()[:] = state_obj.t

	Psa = 1.0/Csa*(V_sa - Vsa0);
	Pad = 1.0/Cad*(V_ad - Vad0);  ###SMS add: for distal aortic compliance
    	Psv = 1.0/Csv*(V_sv - Vsv0);
    	PLV = P_LV;

	if(not isLV):
		Ppv = 1.0/Cpv*(V_pv - Vpv0);
		Ppa = 1.0/Cpa*(V_pa - Vpa0);
    		PRV = P_RV;
    
        printout("Cycle number = "+str(state_obj.cycle)+" cell time = "+str(state_obj.t)+" tstep = "+str(state_obj.tstep)+" dt = "+str(state_obj.dt.dt), comm_me)
        export.printout("Cycle number = "+str(state_obj.cycle)+" cell time = "+str(state_obj.t)+" tstep = "+str(state_obj.tstep)+" dt = "+str(state_obj.dt.dt))

	
	#### For Calculating P_LA ######################################## 
        Ees_la = SimDet["closedloopparam"]["Ees_la"]
        A_la = SimDet["closedloopparam"]["A_la"];
        B_la = SimDet["closedloopparam"]["B_la"];
        V0_la = SimDet["closedloopparam"]["V0_la"];
        Tmax_la = SimDet["closedloopparam"]["Tmax_la"];
        tau_la = SimDet["closedloopparam"]["tau_la"];
        tdelay_la = SimDet["closedloopparam"]["tdelay_la"];
	
	if (state_obj.t < SimDet["HeartBeatLength"] - tdelay_la):
		t_la = state_obj.t + tdelay_la;
        else: 
		t_la = state_obj.t - state_obj.BCL + tdelay_la;

        printout("t_LA = "+str(t_la)+" t_delay_LA = "+str(tdelay_la)+" BCL = "+str(state_obj.BCL)+" t = "+str(state_obj.t), comm_me)
        export.printout("t_LA = "+str(t_la)+" t_delay_LA = "+str(tdelay_la)+" BCL = "+str(state_obj.BCL)+" t = "+str(state_obj.t))

        PLA = et(t_la,Tmax_la,tau_la)*Ees_la*(V_LA - V0_la) + (1 - et(t_la,Tmax_la,tau_la))*A_la*(math.exp(B_la*(V_LA - V0_la)) - 1);

	if(not isLV):
        	Ees_ra = SimDet["closedloopparam"]["Ees_ra"];
        	A_ra = SimDet["closedloopparam"]["A_ra"];
        	B_ra = SimDet["closedloopparam"]["B_ra"];
        	V0_ra = SimDet["closedloopparam"]["V0_ra"];
        	Tmax_ra = SimDet["closedloopparam"]["Tmax_ra"];
        	tau_ra = SimDet["closedloopparam"]["tau_ra"];
        	tdelay_ra = SimDet["closedloopparam"]["tdelay_ra"];

        	PRA = et(t_la,Tmax_ra,tau_ra)*Ees_ra*(V_RA - V0_ra) + (1 - et(t_la,Tmax_ra,tau_ra))*A_ra*(math.exp(B_ra*(V_RA - V0_ra)) - 1);

	        printout("P_pa = "+str(Ppa*0.0075), comm_me)
	        export.printout("P_pa = "+str(Ppa))
	        printout("P_RV = " +str(PRV*0.0075), comm_me)
	        export.printout("P_RV = " +str(PRV ))
	        printout("P_pv = "+str(Ppv*0.0075), comm_me)
	        export.printout("P_pv = "+str(Ppv))
	        printout("P_RA = " +str(PRA*0.0075), comm_me)
	        export.printout("P_RA = " +str(PRA ))


	##################################################################################################################################

        printout("P_sv = "+str(Psv*0.0075), comm_me)
        export.printout("P_sv = "+str(Psv))
        printout("P_LV = " +str(PLV*0.0075), comm_me)
        export.printout("P_LV = " +str(PLV ))
        printout("P_sa = "+str(Psa*0.0075), comm_me)
        export.printout("P_sa = "+str(Psa))
        printout("P_ad = "+str(Pad*0.0075), comm_me)
        export.printout("P_ad = "+str(Pad))
        printout("P_LA = " +str(PLA*0.0075), comm_me)
        export.printout("P_LA = " +str(PLA ))

	#### conditions for Valves#######################################
    	if(PLV <= Psa):
    	     Qav = 0.0;
    	else:
    	     Qav = 1.0/Rav*(PLV - Psa);
    	

    	if(PLV >= PLA):
    	    Qmv = 0.0;
    	else: 
    	    Qmv = 1.0/Rmv*(PLA - PLV);


	if(not isLV):

		if(PRV <= Ppa):  #Pulmonary valve
		     Qpvv = 0.0;
		else:
		     Qpvv = 1.0/Rpvv*(PRV - Ppa);
		
		if(PRV >= PRA): #Tricuspid valve
		     Qtv = 0.0;
		else:
		     Qtv = 1.0/Rtv*(PRA - PRV);


	
	if(isLV):
    		Qsa = 1.0/Rsa*(Psa - Pad);
		Qad = 1.0/Rad*(Pad - Psv); 
		Qsv = 1.0/Rsv*(Psv - PLA);
	else:
 		Qsa = 1.0/Rsa*(Psa - Pad);
		Qad = 1.0/Rad*(Pad - Psv); 
		Qsv = 1.0/Rsv*(Psv - PRA);
	        Qpa = 1.0/Rpa*(Ppa - Ppv);
		Qpv = 1.0/Rpv*(Ppv - PLA);

		# Shunting
		Qlara = Shuntscale*1.0/Rsh*(PLA - PRA); 

	V_LV_prev = V_LV
	V_LA_prev = V_LA
	V_sa_prev = V_sa
	V_ad_prev = V_ad
	V_sv_prev = V_sv
	P_LV_prev = P_LV

	if(not isLV):
		V_RV_prev = V_RV
		V_RA_prev = V_RA
		V_pa_prev = V_pa
		V_pv_prev = V_pv
		P_RV_prev = P_RV

        imp = project( MEmodel_.GetIMP(), FunctionSpace(MEmodel_.mesh_me,'DG',1), form_compiler_parameters={"representation":"uflacs"})
        imp.rename("imp","imp")

        imp2 = project( MEmodel_.GetIMP2(), FunctionSpace(MEmodel_.mesh_me,'DG',1), form_compiler_parameters={"representation":"uflacs"})
        imp2.rename("imp2","imp2")

        if("probepts" in SimDet.keys()):
                x = np.array(SimDet["probepts"])
                probesIMP = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))
                probesIMP(imp)

                probesIMP2 = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))
                probesIMP2(imp2)

                probesIMP3 = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'CG',1))
                probesIMP3(MEmodel_.GetP())


      		# broadcast from proc 0 to other processes
		rank = comm_me_.Get_rank()
        	a = probesIMP3.array()  ## probe will only send to rank =0
        	if(not rank == 0):
			a = np.empty(len(x))

        	comm_me_.Bcast(a, root=0)

                
        printout("Q_lara = " + str(Qlara ), comm_me)
        export.printout("Q_lara = " + str(Qlara))
        printout("Q_lvad = " + str(Qlvad ), comm_me)
        export.printout("Q_lvad = " + str(Qlvad ))
        printout("Q_sv = " + str(Qsv ), comm_me)
        export.printout("Q_sv = " + str(Qsv ))
        printout("Q_av = " + str(Qav ), comm_me)
        export.printout("Q_av = " + str(Qav ))
        printout("Q_sa = "+ str(Qsa), comm_me)
        export.printout("Q_sa = "+ str(Qsa ))
        printout("Q_ad = "+ str(Qsa), comm_me)
        export.printout("Q_ad = "+ str(Qad ))
        printout("Q_mv = " + str(Qmv ), comm_me)
        export.printout("Q_mv = " + str(Qmv ))
	if(not isLV):
        	printout("Q_tv = " + str(Qtv ), comm_me)
        	export.printout("Q_tv = " + str(Qtv ))
        	printout("Q_pa = " + str(Qpa ), comm_me)
        	export.printout("Q_pa = " + str(Qpa ))
        	printout("Q_pv = " + str(Qpv), comm_me)
        	export.printout("Q_pv = " + str(Qpv))
        	printout("Q_pvv = " + str(Qpvv ), comm_me)
        	export.printout("Q_pvv = " + str(Qpvv ))



        printout("Q_LAD = " + str(Qlad), comm_me)
        printout("Q_LCX = " + str(Qlcx), comm_me)
	Qlad = 0
	Qlcx = 0

	# Include LVAD
	H = (Psa - PLV)*0.0075 #Pump head in mmHg
	Qlvad = LVADscale*(-0.0255*(H)+ LVADrpm*0.2 - 0.7)/60 #Flow rate of LVAD

	if(isLV):
		V_LV = V_LV + state_obj.dt.dt*(Qmv - Qav - Qlvad);
		V_sa = V_sa + state_obj.dt.dt*(Qav - Qsa - Qlad - Qlcx + Qlvad);
		V_ad = V_ad + state_obj.dt.dt*(Qsa - Qad);
    		V_sv = V_sv + state_obj.dt.dt*(Qad + Qlad + Qlcx - Qsv);
		V_LA = V_LA + state_obj.dt.dt*(Qsv - Qmv);
	else:
		V_LV = V_LV + state_obj.dt.dt*(Qmv - Qav - Qlvad);
		V_sa = V_sa + state_obj.dt.dt*(Qav - Qsa - Qlad - Qlcx + Qlvad);
		V_ad = V_ad + state_obj.dt.dt*(Qsa - Qad);
    		V_sv = V_sv + state_obj.dt.dt*(Qad + Qlad + Qlcx - Qsv);
		V_RA = V_RA + state_obj.dt.dt*(Qsv - Qtv + Qlara);
	        V_RV = V_RV + state_obj.dt.dt*(Qtv - Qpvv);
		V_pa = V_pa + state_obj.dt.dt*(Qpvv - Qpa);
		V_pv = V_pv + state_obj.dt.dt*(Qpa - Qpv);
		V_LA = V_LA + state_obj.dt.dt*(Qpv - Qmv - Qlara);

	MEmodel_.LVCavityvol.vol = V_LV
	if(not isLV):
		MEmodel_.RVCavityvol.vol = V_RV
		

        printout("V_sv = " + str(V_sv), comm_me)
        export.printout("V_sv = " + str(V_sv))
        printout("V_LV = "  + str(V_LV), comm_me)
        export.printout("V_LV = " + str(V_LV))
        printout("V_sa = " + str(V_sa), comm_me)
        export.printout("V_sa = " + str(V_sa))
        printout("V_ad = " + str(V_ad), comm_me)
        export.printout("V_ad = " + str(V_ad))
        printout("V_LA = "  + str(V_LA ),  comm_me)
        export.printout("V_LA = "  + str(V_LA ))
	if(not isLV):
        	printout("V_pv = " + str(V_pv), comm_me)
        	export.printout("V_pv = " + str(V_pv))
        	printout("V_RV = "  + str(V_RV), comm_me)
        	export.printout("V_RV = "  + str(V_RV))
        	printout("V_pa = " + str(V_pa), comm_me)
        	export.printout("V_pa = " + str(V_pa))
        	printout("V_RA = "  + str(V_RA ),  comm_me)
        	export.printout("V_RA = "  + str(V_RA ))
       
	printout("Solving elasticity", comm_me)
	export.printout("Solving elasticity")
	try:

		solver_elas.solvenonlinear()
		isrestart = 0
		state_obj.dt.dt	= delTat

	except RuntimeError:
		export.printout("Restart time step ********************************************* ")
		V_LV = V_LV_prev
		V_LA = V_LA_prev
		V_sa = V_sa_prev
		V_ad = V_ad_prev
		V_sv = V_sv_prev
		P_LV = P_LV_prev
		if(not isLV):
			V_RV = V_RV_prev
			V_RA = V_RA_prev
			V_pa = V_pa_prev
			V_pv = V_pv_prev
			P_RV = P_RV_prev


		state_obj.tstep = state_obj.tstep - state_obj.dt.dt
		state_obj.dt.dt += 1.0 #Increase time step 
		#state_obj.dt.dt += state_obj.dt.dt/2.0 #Increase time step 
		MEmodel_.Reset()
		EPmodel_.Reset()
		isrestart = 1
		if(state_obj.dt.dt < 1e-3):
			export.printout("Smallest time step reached")
			exit(1);
		continue;

	if(isrestart == 0):
		export.writePV(MEmodel_, state_obj.tstep)
		if(not isLV):
			export.writeP(MEmodel_, [Psv, PLV, Psa, PLA, Ppv, PRV, Ppa, PRA], state_obj.tstep)
			export.writeQ(MEmodel_, [Qav, Qmv, Qsa, Qsv, Qpvv, Qtv, Qpa, Qpv, Qlad, Qlcx, Qlvad], state_obj.tstep)

		else:
			export.writeP(MEmodel_, [Psv, PLV, Psa, PLA], state_obj.tstep)
			export.writeQ(MEmodel_, [Qav, Qmv, Qsa, Qsv, Qlad, Qlcx, Qlvad], state_obj.tstep)


	# Reset phi and r in EP at end of diastole
	if state_obj.t < state_obj.dt.dt:
		EPmodel_.reset();

	printout("Solving FHN", comm_me)
        solver_FHN.solvenonlinear()
	if(isrestart == 0):
	 	MEmodel_.UpdateVar() # For damping
	 	EPmodel_.UpdateVar()

	# Interpolate phi to mechanics mesh
        potential_ref = EPmodel_.interpolate_potential_ep2me_phi(V_me = Function(FunctionSpace(MEmodel_.mesh_me,'CG',1)))
        potential_ref.rename("v_ref", "v_ref")

        potential_me.vector()[:] = potential_ref.vector().array()[:]


        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -
        if(MPI.rank(comm_ep) == 0):
            print 'UPdating isActiveField and tInitiationField'
        MEmodel_.activeforms.update_activationTime(potential_n = potential_me, comm = comm_me)

        F_n = MEmodel_.GetFmat() 
	fstress_DG = project(MEmodel_.Getfstress(), FunctionSpace(MEmodel_.mesh_me,'DG', 0) , form_compiler_parameters={"representation":"uflacs"})
	fstress_DG.rename("fstress", "fstress")
    	if("probepts" in SimDet.keys()):
		probesfstress = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))
		probesfstress(fstress_DG)

	Eul_fiber_BiV_DG = project(fStrain_uL, FunctionSpace(MEmodel_.mesh_me,'DG', 0) , form_compiler_parameters={"representation":"uflacs"})
	Eul_fiber_BiV_DG.rename("Eff", "Eff")
    	if("probepts" in SimDet.keys()):
		probesEul_fiber = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))
		probesEul_fiber(Eul_fiber_BiV_DG)

		probesE_circ_BiV = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))
		probesE_long_BiV = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))
		probesE_radi_BiV = Probes(x.flatten(), FunctionSpace(MEmodel_.mesh_me,'DG',1))


	# postprocess and write 
        #  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - -

	## ----------------- Compute Natural Strain -----------------------------------------------------------------------------
        E_circ_BiV, E_circ_BiV_ = MEmodel_.GetFiberNaturalStrain(F_ED, eCC, AHA_segments)
        E_long_BiV, E_long_BiV_ = MEmodel_.GetFiberNaturalStrain(F_ED, eLL, AHA_segments)
        E_radi_BiV, E_radi_BiV_ = MEmodel_.GetFiberNaturalStrain(F_ED, eRR, AHA_segments)
	## --------------------------------------------------------------------------------------------------------------------



	E_circ_BiV_DG = project(E_circ_BiV_, FunctionSpace(MEmodel_.mesh_me,'DG', 0) , form_compiler_parameters={"representation":"uflacs"})
	E_circ_BiV_DG.rename("Ecc", "Ecc")
    	if("probepts" in SimDet.keys()):
		probesE_circ_BiV(E_circ_BiV_DG)

    	E_long_BiV_DG = project(E_long_BiV_, FunctionSpace(MEmodel_.mesh_me,'DG', 0) , form_compiler_parameters={"representation":"uflacs"})
	E_long_BiV_DG.rename("Ell", "Ell")
    	if("probepts" in SimDet.keys()):
		probesE_long_BiV(E_long_BiV_DG)

	E_radi_BiV_DG = project(E_radi_BiV_, FunctionSpace(MEmodel_.mesh_me,'DG', 0) , form_compiler_parameters={"representation":"uflacs"})
	E_radi_BiV_DG.rename("Err", "Err")
    	if("probepts" in SimDet.keys()):
		probesE_radi_BiV(E_radi_BiV_DG)


        if(cnt % SimDet["writeStep"] == 0.0):

	    export.writetpt(MEmodel_, state_obj.tstep);
	    export.hdf.write(MEmodel_.GetDisplacement(), "ME/u", writecnt)
	    export.hdf.write(potential_ref, "ME/potential_ref", writecnt)
	    export.hdf.write(E_circ_BiV_DG, "ME/Ecc", writecnt)
	    export.hdf.write(E_long_BiV_DG, "ME/Ell", writecnt)
	    export.hdf.write(E_radi_BiV_DG, "ME/Err", writecnt)
	    export.hdf.write(Eul_fiber_BiV_DG, "ME/Eff", writecnt)
	    export.hdf.write(fstress_DG, "ME/fstress", writecnt)
	    export.hdf.write(imp, "ME/imp", writecnt)
	    export.hdf.write(imp2, "ME/imp2",  writecnt)
	    export.hdf.write(MEmodel_.GetP(), "ME/imp_constraint", writecnt)

            export.hdf.write(EPmodel_.getphivar(), "EP/phi", writecnt)
            export.hdf.write(EPmodel_.getrvar(), "EP/r", writecnt)
            export.hdf.write(potential_ref, "EP/potential_ref", writecnt)

	    writecnt += 1


    	if("probepts" in SimDet.keys()):
		fIMP = probesIMP.array()
		fIMP2 = probesIMP2.array()
		fIMP3 = probesIMP3.array()
		fStress = probesfstress.array()
		fStrain_vals = probesEul_fiber.array()
		E_circ_BiV = probesE_circ_BiV.array()
		E_long_BiV = probesE_long_BiV.array()
		E_radi_BiV = probesE_radi_BiV.array()

		export.writeIMP(MEmodel_, state_obj.tstep, fIMP)
		export.writeIMP2(MEmodel_, state_obj.tstep, fIMP2)
		export.writeIMP3(MEmodel_, state_obj.tstep, fIMP3)
		export.writefStress(MEmodel_, state_obj.tstep, fStress)
		export.writefStrain(MEmodel_, state_obj.tstep, fStrain_vals)
		export.writeCStrain(MEmodel_, state_obj.tstep, E_circ_BiV)
		export.writeLStrain(MEmodel_, state_obj.tstep, E_long_BiV)
		export.writeRStrain(MEmodel_, state_obj.tstep, E_radi_BiV)

        cnt += 1
	


#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
if __name__ == "__main__":

    print 'Testing...'
    run_BiV_TimedGuccione(IODet=IODetails, SimDet=SimDetails)

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
