import os
from vtk.util import numpy_support
from postprocessdatalibBiV2 import *

def postprocessdata(IODet, SimDet):

    directory = IODet["outputfolder"] + '/' 
    casename =  IODet["caseID"] 
    BCL = SimDet["HeartBeatLength"]
    cycle = SimDet["closedloopparam"]["stop_iter"]

    for ncycle in range(cycle-1,cycle):

	filename = directory+casename+"/"+"BiV_PV.txt"
	tptt, LVP, LVV, RVP, RVV, Qmv = extract_PV(filename, BCL, ncycle, False)

	filename = directory+casename+"/"+"BiV_Q.txt"
 	tptt, Qav, Qmv, Qsa, Qsv, Qpvv, Qtv, Qpa, Qpv, Qlvad = extract_Q(filename, BCL, ncycle)

	filename = directory+casename+"/"+"BiV_P.txt"
 	tptt, Psv, PLV, Psa, PLA, Ppv, PRV, Ppa, PRA = extract_P(filename, BCL, ncycle)

	LVESP, LVESV = extractESP(LVP, LVV)
	RVESP, RVESV = extractESP(RVP, RVV)

	LVEDP, LVEDV = extractEDP(LVP, LVV)
	RVEDP, RVEDV = extractEDP(RVP, RVV)

	SBP = max(Psa)*0.0075
	DBP = min(Psa)*0.0075

	mPAP = np.mean(Ppa)*0.0075

	print filename
	print "LVEF = ", (LVEDV - LVESV)/LVEDV*100, " LVEDV = ", LVEDV, " LVESV = ", LVESV, " LVEDP = ", LVEDP
	print "RVEF = ", (RVEDV - RVESV)/RVEDV*100, " RVEDV = ", RVEDV, " RVESV = ", RVESV, " RVEDP = ", RVEDP
	print "SBP = ", SBP, " DBP = ", DBP, " mPAP = ", mPAP
	print "Peak LV pressure = ", max(LVP)

	homo_directory = directory+casename+"/"

	tpt_array = readtpt(homo_directory + "tpt.txt")
	ind = np.where((tpt_array>(ncycle)*BCL)*(tpt_array<(ncycle+1)*BCL))
	tpt = tpt_array[ind]

	## Extract Surfaces and Set Material region
	isparallel = False
	if(os.path.exists(homo_directory + "facetboundaries_ep000000.pvtu")):
		isparallel = True
	LVendo, RVendo, Epi = GetSurfaces(homo_directory + "", "facetboundaries_ep", "f", isparallel)
	matid  = setMaterialRegion(homo_directory,  LVendo, RVendo, Epi)
	vtk_py.writeXMLUGrid(matid, "matid"+".vtu")

	## Get Point cloud for probing
	ptcloud, radialpos, matid_arr, vtkradialpos = getpointclouds(homo_directory , LVendo, RVendo, Epi, matid)
	vtk_py.writeXMLPData(vtkradialpos, casename+".vtp")

	# Get transmural variation of IMP
	index = find_nearest(tpt, tptt[ np.argmax(PLV)]) # Find ID correspond to peak LV pressure
	imp = probeqty(homo_directory,  'ME/imp_constraint', ptcloud, ind, index)
	imp = imp*0.0075
	

	## Get transmural variation of WD
	Sff = probetimeseries(homo_directory, "ME/fstress", ptcloud, ind, "DG", 0)
	Eff = probetimeseries(homo_directory, "ME/Eff", ptcloud, ind, "DG", 0)
	WD = np.array([-1.0*np.trapz(Sff[:,i]*0.0075, Eff[:,i]) for i in range(0,len(Sff[1,:]))])

	# Convert to vtp flie
	for i in range(0,len(Sff[:,1])):
		pdata = vtk.vtkPolyData()
		pdata.DeepCopy(vtkradialpos)
		Sff_VTK_data = numpy_support.numpy_to_vtk(num_array=0.0075*Sff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
		Sff_VTK_data.SetName("fstress_")
		pdata.GetPointData().AddArray(Sff_VTK_data)
		Eff_VTK_data = numpy_support.numpy_to_vtk(num_array=Eff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
		Eff_VTK_data.SetName("Eff_")
		pdata.GetPointData().AddArray(Eff_VTK_data)
		WD_VTK_data = numpy_support.numpy_to_vtk(num_array=WD.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
		WD_VTK_data.SetName("WD_")
		pdata.GetPointData().AddArray(WD_VTK_data)


	## Get Ecc
	Ecc = probetimeseries(homo_directory, "ME/Ecc", ptcloud, ind, "DG", 0)
	#Ecc = probetimeseries(homo_directory + "active", "Ecc", "Ecc", ptcloud, isparallel, ind)
	peakEcc = np.max(np.abs(np.mean(Ecc, axis=1)*100))
	print "Peak Ecc = ", peakEcc
	
	# Get Ell
	Ell = probetimeseries(homo_directory, "ME/Ell", ptcloud, ind, "DG", 0)
	#Ell = probetimeseries(homo_directory + "active", "Ell", "Ell", ptcloud, isparallel, ind)
	peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))
	print "Peak Ell = ", peakEll

        np.savez(directory+casename+"/"+casename+".npz", \
		 tptt    = tptt,\
		 LVP     = LVP,\
		 LVV     = LVV,\
		 RVP     = RVP,\
		 RVV     = RVV,\
		 Qav     = Qav,\
		 Qmv     = Qmv,\
		 Qsa     = Qsa,\
		 Qsv     = Qsv,\
		 Qpvv    = Qpvv,\
		 Qtv     = Qtv,\
		 Qpa     = Qpa,\
		 Qpv     = Qpv,\
		 Qlvad   = Qlvad,\
		 Psv     = 0.0075*Psv,\
		 PLV     = 0.0075*PLV,\
		 Psa     = 0.0075*Psa,\
		 PLA     = 0.0075*PLA,\
		 Ppv     = 0.0075*Ppv,\
		 PRV     = 0.0075*PRV,\
		 Ppa     = 0.0075*Ppa,\
		 PRA     = 0.0075*PRA,\
		 LVESP = LVESP, LVESV = LVESV, LVEDP = LVEDP, LVEDV = LVEDV, SBP = SBP, DBP = DBP,\
		 RVESP = RVESP, RVESV = RVESV, RVEDP = RVEDP, RVEDV = RVEDV, mPAP = mPAP,\
		 imp          =imp,\
		 radialpos    =radialpos,\
		 Eff	      =Eff, \
		 Sff	      =Sff, \
                 WD           =WD,\
		 Ecc          =Ecc,\
		 Ell          =Ell,\
		 BCL	      =BCL,\
                 tpt          =tpt,\
		 ncycle       =ncycle,\
	         matid_arr    =matid_arr
		)


