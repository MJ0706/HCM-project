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
    printout("Pressure = " +  str(MEmodel_.GetLVP()*0.0075) +  " Vol = " + str(MEmodel_.GetLVV()), comm_me)
    
    # Unloading LV to get new reference geometry
    isunloading = False
    if("isunloading" in SimDet.keys()):
    	if (SimDet["isunloading"] is True)  : 
		isunloading = True

    fdata = dolfin.File(outputfolder + "/deformation_gradient/F.pvd") 

    if(isunloading):

        printout("Start UnLoading", comm_me)
	#V_LV_target = MEmodel_.GetLVV()

	if("unloadparam" in SimDet.keys()):
		unloadparam = SimDet["unloadparam"]
	else:
		unloadparam = {};
			
    	nloadstep_, volinc_, V_LV_target, solver_elas = MEmodel_.unloading(unloadparam)
	printout("Target EDV = " + str(V_LV_target), comm_me)

        printout("Finish UnLoading and Reloading", comm_me)


 	export.writePV(MEmodel_, 0);
    	export.hdf.write(MEmodel_.mesh_me, "ME/mesh")
    	export.hdf.write(EPmodel_.mesh_ep, "EP/mesh")

	MEmodel_.LVCavityvol.vol = MEmodel_.GetLVV()
	V_LV_unload = MEmodel_.GetLVV()

	for it in np.arange(0,nloadstep):
		
		if("V_LV_target" in SimDet["closedloopparam"].keys()):
    	    		V_LV_target = SimDet["closedloopparam"]["V_LV_target"]
    	    		MEmodel_.LVCavityvol.vol += (V_LV_target - V_LV_unload)/nloadstep
		else:
    	    		MEmodel_.LVCavityvol.vol += 2.0

        	#MEmodel_.LVCavityvol.vol += (V_LV_target - V_LV_unload)/nloadstep
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

    	for lmbda_value in range(0, nloadstep):

    	    if("V_LV" in SimDet["closedloopparam"].keys()):
    	    	V_LV_target = SimDet["closedloopparam"]["V_LV_target"]
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
    	    	printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.GetLVV()) + \
 								      " LVP = " + str(MEmodel_.GetLVP()*0.0075)  \
                                                                      , comm_me)
    	    	export.printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.GetLVV()))
	    else:
    	    	printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.GetLVV()) + \
                                                                      " LVP = " + str(MEmodel_.GetLVP()*0.0075) + \
								      " RVV = " + str(MEmodel_.GetRVV()) + \
                                                                      " RVP = " + str(MEmodel_.GetRVP()*0.0075)  \
 								      , comm_me)
    	    	export.printout("Loading phase step = " + str(lmbda_value) + " LVV = " + str(MEmodel_.GetLVV())+ \
									     " RVV = " + str(MEmodel_.GetRVV()))

    	                        	
    	    export.writePV(MEmodel_, 0);
    	    export.hdf.write(MEmodel_.GetDisplacement(), "ME/u_loading", lmbda_value)

    	    F_ED.vector()[:] = project(MEmodel_.GetFmat(), MEmodel_.TF, solver_type='mumps').vector().array()[:] 
	    #export.hdf.write(F_ED, "ME/Fed", lmbda_value)

	    print F_ED.vector().array()
            #fdata << project(F_ED, dolfin.TensorFunctionSpace(MEmodel_.mesh_me, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
	    

    #export.hdf.close()

    if("isunloadingonly" in SimDet.keys()):
    	if (SimDet["isunloadingonly"] is True)  : 
    		export.hdf.close()
		exit()




#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
if __name__ == "__main__":

    print 'Testing...'
    run_BiV_TimedGuccione(IODet=IODetails, SimDet=SimDetails)

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
