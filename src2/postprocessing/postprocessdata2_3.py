from postprocessdatalib2 import *

def postprocessdata(IODet, SimDet):

    directory = IODet["outputfolder"] + '/' 
    casename =  IODet["caseID"] 
    BCL = SimDet["HeartBeatLength"]
    cycle = SimDet["closedloopparam"]["stop_iter"]


    for ncycle in range(cycle-1,cycle):
    
         filename = directory+casename+"/"+"BiV_PV.txt"
         homo_tptt, homo_LVP, homo_LVV, homo_Qmv = extract_PV(filename, BCL, ncycle)

         filename = directory+casename+"/"+"BiV_Q.txt"
         homo_tptt, homo_Qao, homo_Qmv, homo_Qper, homo_Qla, homo_Qlad, homo_Qlcx, Q_lvad = extract_Q(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_P.txt"
         homo_tptt, homo_Pven, homo_LVPP, homo_Part, homo_PLA = extract_P(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_IMP_InC.txt"
         homo_tpt_IMP, homo_IMP = extract_probe(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_fiberStrain.txt"
         homo_tpt_Eff, homo_Eff = extract_probe(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_fiberStress.txt"
         homo_tpt_Sff, homo_Sff = extract_probe(filename, BCL, ncycle)
    
         ESP, ESV = extractESP(homo_LVP, homo_LVV)
    
         EDP, EDV = extractEDP(homo_LVP, homo_LVV)
    
         SBP = max(homo_Part)*0.0075
         DBP = min(homo_Part)*0.0075
    
         print filename
         print "EF = ", (max(homo_LVV) - min(homo_LVV))/max(homo_LVV), " EDV = ", max(homo_LVV), " ESV = ", min(homo_LVV),\
         		" EDP = ", EDP, " SBP = ", SBP, " DBP = ", DBP
    
         print "Peak LV pressure = ", max(homo_LVP)
    
         homo_directory = directory+casename+"/"
    
         tpt_array = readtpt(homo_directory + "tpt.txt")
         ind = np.where((tpt_array>(ncycle)*BCL)*(tpt_array<(ncycle+1)*BCL))
         tpt = tpt_array[ind]
         
         ## Get Point cloud for probing
         ptcloud, radialpos, vtkpos, longpos = getpointclouds(homo_directory, clipoffset=5e-1, npts=10000)
         vtk_py.writeXMLPData(vtkpos, casename+".vtp")
    
         # Get transmural variation of IMP
         index = find_nearest(tpt, homo_tptt[ np.argmax(homo_LVPP)]) # Find ID correspond to peak LV pressure
         imp = probeqty(homo_directory,  'ME/imp_constraint', ptcloud, ind, index)
         imp = imp*0.0075
         #print imp
	 print len(imp)
	 #print np.shape(imp)

	 #Get imp
	 IMP_1 = probetimeseries(homo_directory, "ME/imp_constraint", ptcloud, ind, "CG", 1)
         IMP_1 = np.array(IMP_1*0.0075)
	 peakIMP = np.max(np.abs(np.mean(IMP_1)))
	 print "Peak IMP = ", peakIMP
	 print len(IMP_1)
	 print IMP_1
	 
         # Convert to vtp flie
	 Tff_t = probetimeseries(homo_directory, "ME/fstress", ptcloud, ind, "DG", 0)
         Sff_t = probetimeseries(homo_directory, "ME/fstressPk2", ptcloud, ind, "DG", 0)
         Sff_t_2 = probetimeseries(homo_directory, "ME/fstressPk2_2", ptcloud, ind, "DG", 0)
         Sff_t_3 = probetimeseries(homo_directory, "ME/fstressPk2_3", ptcloud, ind, "DG", 0)
         Sff_t_4 = probetimeseries(homo_directory, "ME/fstressPk2_4", ptcloud, ind, "DG", 0)
         Sss_t = probetimeseries(homo_directory, "ME/sstressPk2", ptcloud, ind, "DG", 0)
         Snn_t = probetimeseries(homo_directory, "ME/nstressPk2", ptcloud, ind, "DG", 0)

	 Tff_ = probetimeseries(homo_directory, "ME/Sff", ptcloud, ind, "DG", 0)
         Sff_ = probetimeseries(homo_directory, "ME/Pk2ff", ptcloud, ind, "DG", 0)
	 Pff_ = probetimeseries(homo_directory, "ME/fstressact", ptcloud, ind, "DG", 0)
    
         ## Get transmural variation of WD
	 Lff = probetimeseries(homo_directory, "ME/lff", ptcloud, ind, "DG", 0)
         
         Eff = probetimeseries(homo_directory, "ME/Eff", ptcloud, ind, "DG", 0)
         Eff_2 = probetimeseries(homo_directory, "ME/Eff2", ptcloud, ind, "DG", 0)
         Eff_3 = probetimeseries(homo_directory, "ME/Eff3", ptcloud, ind, "DG", 0)
         Eff_4 = probetimeseries(homo_directory, "ME/Eff4", ptcloud, ind, "DG", 0)
         Ess = probetimeseries(homo_directory, "ME/Ess", ptcloud, ind, "DG", 0)
         Enn = probetimeseries(homo_directory, "ME/Enn", ptcloud, ind, "DG", 0)

         WD = np.array([-1.0*np.trapz(Sff_t[:,i]*0.0075, Eff[:,i]) for i in range(0,len(Sff_t[1,:]))])
         WD_2 = np.array([-1.0*np.trapz(Sff_t_2[:,i]*0.0075, Eff_2[:,i]) for i in range(0,len(Sff_t_2[1,:]))])
         WD_3 = np.array([-1.0*np.trapz(Sff_t_3[:,i]*0.0075, Eff_3[:,i]) for i in range(0,len(Sff_t_3[1,:]))])
         WD_4 = np.array([-1.0*np.trapz(Sff_t_4[:,i]*0.0075, Eff_4[:,i]) for i in range(0,len(Sff_t_4[1,:]))])
         WS = np.array([-1.0*np.trapz(Sss_t[:,i]*0.0075, Ess[:,i]) for i in range(0,len(Sss_t[1,:]))])
         WN = np.array([-1.0*np.trapz(Snn_t[:,i]*0.0075, Enn[:,i]) for i in range(0,len(Snn_t[1,:]))])

	 #Sffactive = probetimeseries(homo_directory, "ME/fstressact", ptcloud, ind, "DG", 0)
    
    
         ## Get Ecc
         Ecc = probetimeseries(homo_directory, "ME/Ecc", ptcloud, ind, "DG", 0)
         peakEcc = np.max(np.abs(np.mean(Ecc, axis=1)*100))
         print "Peak Ecc = ", peakEcc
         print len(Ecc)
	 print Ecc
	 
         # Get Ell
         Ell = probetimeseries(homo_directory, "ME/Ell", ptcloud, ind, "DG", 0)
         peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))
         print "Peak Ell = ", peakEll


	 # Get lff
         Lff = probetimeseries(homo_directory, "ME/lff", ptcloud, ind, "DG", 0)
         peakLff = np.max(np.abs(np.mean(Lff, axis=1)*100))
         print "Peak lff = ", peakLff

	 # Get lss
         Lss = probetimeseries(homo_directory, "ME/lss", ptcloud, ind, "DG", 0)
         peakLss = np.max(np.abs(np.mean(Lss, axis=1)*100))
         print "Peak lss = ", peakLss

	 # Get lnn
         Lnn = probetimeseries(homo_directory, "ME/lnn", ptcloud, ind, "DG", 0)
         peakLnn = np.max(np.abs(np.mean(Lnn, axis=1)*100))
         print "Peak lnn = ", peakLnn

	 # Get Sff
         #Sff_ = probetimeseries(homo_directory, "ME/Sff", ptcloud, ind, "DG", 0)
         peakSff_ = np.max(np.abs(np.mean(Sff_, axis=1)*100))
         print "Peak Sff_ = ", peakSff_

	 # Get Sss
         Sss = probetimeseries(homo_directory, "ME/Sss", ptcloud, ind, "DG", 0)
         peakSss = np.max(np.abs(np.mean(Sss, axis=1)*100))
         print "Peak Sss = ", peakSss

	 # Get Snn
         Snn = probetimeseries(homo_directory, "ME/Snn", ptcloud, ind, "DG", 0)
         peakSnn = np.max(np.abs(np.mean(Snn, axis=1)*100))
         print "Peak Snn = ", peakSnn

	 
         cnt = 0
         for i in np.arange(0,len(Sff_t[:,1]),2):
         	pdata = vtk.vtkPolyData()
         	pdata.DeepCopy(vtkpos)
         	Sff_VTK_data = numpy_support.numpy_to_vtk(num_array=0.0075*Sff_t[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	Sff_VTK_data.SetName("fstress_")
         	pdata.GetPointData().AddArray(Sff_VTK_data)
         	Eff_VTK_data = numpy_support.numpy_to_vtk(num_array=Eff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	Eff_VTK_data.SetName("Eff_")
         	pdata.GetPointData().AddArray(Eff_VTK_data)
         	WD_VTK_data = numpy_support.numpy_to_vtk(num_array=WD.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	WD_VTK_data.SetName("WD_")
         	pdata.GetPointData().AddArray(WD_VTK_data)
		Ell_VTK_data = numpy_support.numpy_to_vtk(num_array=Ell[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
		Ell_VTK_data.SetName("Ell_")
         	pdata.GetPointData().AddArray(Ell_VTK_data)
		Ecc_VTK_data = numpy_support.numpy_to_vtk(num_array=Ecc[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
		Ecc_VTK_data.SetName("Ecc_")
         	pdata.GetPointData().AddArray(Ecc_VTK_data)

         	vtk_py.writeXMLPData(pdata, directory+casename+"/"+casename+"data"+str(cnt)+".vtp")
                cnt += 1
     
    
         np.savez(directory+casename+"/"+casename+".npz", \
         	 homo_tptt    = homo_tptt,\
         	 homo_LVP     = homo_LVP,\
         	 homo_LVV     = homo_LVV,\
         	 homo_Qmv     = homo_Qmv,\
         	 homo_Qao     = homo_Qao,\
         	 homo_Qper    = homo_Qper,\
         	 homo_Qla     = homo_Qla,\
         	 homo_Qlad    = homo_Qlad,\
         	 homo_Pven    = 0.0075*homo_Pven,\
         	 homo_LVPP    = 0.0075*homo_LVPP,\
         	 homo_Part    = 0.0075*homo_Part,\
         	 homo_PLA     = 0.0075*homo_PLA,\
                 homo_tpt_IMP = 0.0075*homo_tpt_IMP,\
         	 homo_IMP     = homo_IMP,\
         	 homo_tpt_Eff = homo_tpt_Eff,\
         	 homo_Eff     = homo_Eff,\
                 homo_tpt_Sff = homo_tpt_Sff,\
         	 homo_Sff     = homo_Sff,\
         	 ESP = ESP, ESV = ESV, EDP = EDP, EDV = EDV, SBP = SBP, DBP = DBP,\
         	 #Qtotal       = Qtotal,\
         	 imp          =imp,\
         	 IMP_1          =IMP_1,\
         	 radialpos    =radialpos,\
         	 longpos      =longpos,\
         	 Eff	      =Eff, \
         	 Eff_2	      =Eff_2, \
         	 Eff_3	      =Eff_3, \
         	 Eff_4	      =Eff_4, \
         	 Ess	      =Ess, \
         	 Enn	      =Enn, \
         	 Sff_t	      =Sff_t, \
         	 Sff_t_2      =Sff_t_2, \
         	 Sff_t_3      =Sff_t_3, \
         	 Sff_t_4      =Sff_t_4, \
         	 Sss_t	      =Sss_t, \
         	 Snn_t	      =Snn_t, \
         	 Tff_t	      =Tff_t, \
         	 Sff_         =Sff_, \
         	 Pff_         =Pff_, \
         	 Tff_         =Tff_, \
                 WD           =WD,\
                 WD_2         =WD_2,\
                 WD_3         =WD_3,\
                 WD_4         =WD_4,\
                 WS           =WS,\
                 WN           =WN,\
         	 Ecc          =Ecc,\
         	 Ell          =Ell,\
         	 Lff          =Lff,\
         	 Lss          =Lss,\
         	 Lnn          =Lnn,\
         	 Sss_          =Sss,\
         	 Snn_          =Snn,\
         	 BCL	      =BCL,\
                 tpt          =tpt,\
         	 ncycle       =ncycle
         	)
    
	
	
