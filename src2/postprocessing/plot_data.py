from postprocessdatalib2 import *
from matplotlib import pylab as plt
from matplotlib import rc
import numpy as np

# Plot effects of contractility cases

rc('text', usetex=False)#True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})

BCLs = [
       	1000,\
       	1000,\
       ]
cycle = 5


casenames  = [
	      #'HCM2_7',\
	      'HCM2_7',\
	      #'HCM2_7_testa',\
	     ]

labels  = [
	      #'HCM2_7',\
	      'HCM2_7_test',\
	      'HCM2_7_testa',\
	     ]



colors=["k",
        "r",
        "b"]

PVoutfilename = "comparePV.png"
LVPoutfilename = "compareLVP.png"
LVVoutfilename = "compareLVV.png"
Qmvoutfilename = "compareQmv.png"
Vmvoutfilename = "compareVmv.png"
Qladoutfilename = "compareQlad.png"
IMPoutfilename = "compareIMP.png"
Qdisfilename = "compareQdist.png"
IMPdisfilename = "compareIMPdist.png"
WDoutfilename = "compareWD.png"
WDdistfilename = "compareWDdist.png"
RelQdisfilename = "comparerelQdist.png"

caseno = 0
MV_area = 2.5 # Assuming MV_area = 4cm^2 and Vmax = Q/(A/2) based on poiseulle
#Qtotal_arr = []

for (casename, label, color) in zip(casenames, labels, colors):

	#radialpos    = np.load(casename+".npz")['radialpos']
	#Sff 	     = np.load(casename+".npz")['Sff']
	#WD           = np.load(casename+".npz")['WD']
	#Eff 	     = np.load(casename+".npz")['Eff']
	#Ecc 	     = np.load(casename+".npz")['Ecc']
	#Ell 	     = np.load(casename+".npz")['Ell']
	#imp 	     = np.load(casename+".npz")['imp']
	##Qtotal       = np.load(casename+".npz")['Qtotal']
	homo_tptt    = np.load(casename+".npz")['homo_tptt']
	homo_LVP     = np.load(casename+".npz")['homo_LVP']
	homo_LVV     = np.load(casename+".npz")['homo_LVV']
	homo_Qmv     = np.load(casename+".npz")['homo_Qmv']
	homo_Qao     = np.load(casename+".npz")['homo_Qao']
	homo_Qper    = np.load(casename+".npz")['homo_Qper']
	homo_Qla     = np.load(casename+".npz")['homo_Qla']
	homo_Qlad    = np.load(casename+".npz")['homo_Qlad']
	homo_Pven    = np.load(casename+".npz")['homo_Pven']
	homo_LVPP    = np.load(casename+".npz")['homo_LVPP']
	homo_Part    = np.load(casename+".npz")['homo_Part']
	homo_PLA     = np.load(casename+".npz")['homo_PLA']
	BCL          = np.load(casename+".npz")['BCL']
        tpt          = np.load(casename+".npz")['tpt']
        ncycle       = np.load(casename+".npz")['ncycle']


	# Plot PV Loop
	plt.figure(1)
	plt.plot(homo_LVV, homo_LVP, label=label, color=color)

	# Plot Pressure waveform
	plt.figure(2)
	plt.plot(homo_tptt - homo_tptt[0], homo_LVPP, label=label + r"$P_{LV}$", color=color, linestyle='-', linewidth=2)
	plt.plot(homo_tptt - homo_tptt[0], homo_Pven, label=label + r"$P_{ven}$", color=color, linestyle='--', linewidth=2)
	plt.plot(homo_tptt - homo_tptt[0], homo_PLA, label=label + r"$P_{LA}$", color=color, linestyle='-.', linewidth=2)
	plt.plot(homo_tptt - homo_tptt[0], homo_Part, label=label + r"$P_{a,p}$", color=color, linestyle=':', linewidth=2)

	# Plot Volume waveform
	plt.figure(3)
	plt.plot(homo_tptt - homo_tptt[0], homo_LVV , label=label, color=color)

	# Plot Mitral valve velocity
	plt.figure(4)
	plt.plot(homo_tptt - homo_tptt[0], homo_Qmv*1000 , label=label, color=color)

	# Plot Mitral valve velocity normalized with MV area
	plt.figure(5)
	plt.plot(homo_tptt - homo_tptt[0], homo_Qmv/MV_area*1000 , label=label, color=color)
	'''
	# Get transmural variation of IMP
	plt.figure(9)
	nbins = 8
	n, __ = np.histogram(radialpos, bins=nbins)
	sy, __ = np.histogram(radialpos, bins=nbins, weights=imp)
	sy2, __ = np.histogram(radialpos, bins=nbins, weights=imp*imp)
	mean_imp = sy / n
	std_imp = np.sqrt(sy2/n - mean_imp*mean_imp)

	plt.plot(radialpos, imp, 'x', alpha=0.1, color=color)
	plt.errorbar((__[1:] + __[:-1])/2, mean_imp, yerr=std_imp, color=color, label=label)

	# Get transmural variation of IMP/LVP
	plt.figure(10)
	sy, _ = np.histogram(radialpos, bins=nbins, weights=imp/max(homo_LVPP))
	sy2, _ = np.histogram(radialpos, bins=nbins, weights=imp*imp/(max(homo_LVPP)*max(homo_LVPP)))
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)

	plt.plot(radialpos, imp/max(homo_LVPP), 'x', alpha=0.1, color=color)
	plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, color=color, label=label)
	plt.ylim(-0.25,1.25)
	plt.ylabel(r'IMP/LVP', fontsize=16)
	plt.xlabel(r"Radial position $R$")

	#if(isexpt):
#		plt.plot(IMP_Expt[:,0], IMP_Expt[:,1], 'o', color=color, label="Expt")

	if(len(casenames)<=1):
		plt.twinx().errorbar((__[1:] + __[:-1])/2, mean_imp, yerr=std_imp, color=color, label=label)
		plt.ylim(-0.25*max(homo_LVPP),1.25*max(homo_LVPP))
		plt.ylabel(r'IMP (mmHg)', fontsize=16)


	# Get transmural variation of WD
	plt.figure(11)
	n, _ = np.histogram(radialpos, bins=nbins)
	sy, _ = np.histogram(radialpos, bins=nbins, weights=WD)
	sy2, _ = np.histogram(radialpos, bins=nbins, weights=WD*WD)
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)

	plt.plot(radialpos, WD, 'x', alpha=0.1, color=color)
	plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, color=color, label=label)

	# Get endo midwall epi value for IMP, WD
	plt.figure(12)
	n, _ = np.histogram(radialpos, bins=3)
	sy, _ = np.histogram(radialpos, bins=3, weights=WD)
	sy2, _ = np.histogram(radialpos, bins=3, weights=WD*WD)
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)
	print "Work done at endo = ", mean[0] , " +/- ", std[0]
	print "Work done at midwall = ", mean[1] , " +/- ", std[1]
	print "Work done at epi = ", mean[2] , " +/- ", std[2]
	
	
	n, _ = np.histogram(radialpos, bins=3)
	sy, _ = np.histogram(radialpos, bins=3, weights=imp)
	sy2, _ = np.histogram(radialpos, bins=3, weights=imp*imp)
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)
	print "IMP at endo = ", mean[0] , " +/- ", std[0]
	print "IMP at midwall = ", mean[1] , " +/- ", std[1]
	print "IMP at epi = ", mean[2] , " +/- ", std[2]	'''

	# Plot Wf vs Q ratio
#	n, _ = np.histogram(radialpos, bins=4)
#	sy, _ = np.histogram(radialpos, bins=4, weights=WD)
#	sy2, _ = np.histogram(radialpos, bins=4, weights=WD*WD)
#	mean = sy / n
#	std = np.sqrt(sy2/n - mean*mean)
#
#	cnt = 0
#	loc_array = []
#	for m_ in mean:
#		if(cnt == 0):
#			loc = " Endo"
#			loc_array.append(loc)
#		elif(cnt == 3):
#			loc = " Epi"
#			loc_array.append(loc)
#		else:
#			loc = " Midwall " + str(cnt+1)
#			loc_array.append(loc)
#
#		cnt += 1
#
#	width = 0.8/len(casenames)
#	xloc =  np.arange(len(mean))+caseno*width
#	plt.bar(xloc, mean/Qtotal_arr[caseno], yerr=std, color=color, label=label,alpha=0.5, width=width)
#	if(caseno == 0):
#		plt.xticks(np.arange(len(mean))+0.5*len(casenames)*width - 0.5*width,  loc_array)
#
	'''	
	# Plot Sff vs Eff loop
	plt.figure(13)
	nbins = 3
	n, _ = np.histogram(radialpos, bins=nbins)
	Sff_bin_array = []
	Eff_bin_array = []
	Sff_err_bin_array = []
	Eff_err_bin_array = []
	
	for p in range(0, len(Sff[:,1])):
		Sff_av, _ = np.histogram(radialpos, bins=nbins, weights=Sff[p,:])
		Eff_av, _ = np.histogram(radialpos, bins=nbins, weights=Eff[p,:])
	
		Sff_err, _ = np.histogram(radialpos, bins=nbins, weights=Sff[p,:]*Sff[p,:])
		Sff_err = np.sqrt(Sff_err/n - Sff_av*Sff_av/(n*n))
	
		Eff_err, _ = np.histogram(radialpos, bins=nbins, weights=Eff[p,:]*Eff[p,:])
		Eff_err = np.sqrt(Eff_err/n - Eff_av*Eff_av/(n*n))
	
		Sff_bin_array.append(Sff_av/n)
		Eff_bin_array.append(Eff_av/n)
	
		Sff_err_bin_array.append(Sff_err)
		Eff_err_bin_array.append(Eff_err)
	
	
	Sff_bin_array = np.array(Sff_bin_array)
	Eff_bin_array = np.array(Eff_bin_array)
	
	Sff_err_bin_array = np.array(Sff_err_bin_array)
	Eff_err_bin_array = np.array(Eff_err_bin_array)
	
	for p in range(0, nbins):
		if(p == 0):
			location = "Endo"
			linestyle= "-"
		elif(p == nbins-1):
			location = "Epi"
			linestyle = "--"
		else:
			location = "Midwall"
			linestyle = "-."
	
		plt.plot(Eff_bin_array[:,p], Sff_bin_array[:,p]*0.0075, label=label + location, linestyle=linestyle, color=color)

	# Get Ecc
	plt.figure(14)
	peakEcc = np.max(np.abs(np.mean(Ecc, axis=1)*100))
	plt.plot(tpt - (ncycle-1)*BCL, np.mean(Ecc, axis=1)*100, linestyle='-', color=color, \
		 label=label+"("+str(int(peakEcc))+"%)",\
		)
	print "Peak Ecc = ", peakEcc
	
	# Get Ell
	plt.figure(15)
	peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))
	plt.plot(tpt - (ncycle-1)*BCL, np.mean(Ell, axis=1)*100, linestyle='-', color=color, \
		label=label+"("+str(int(peakEll))+"%)",\
		)
	print "Peak Ell = ", peakEll

	'''




   	caseno += 1
 
	

	

plt.figure(1)
plt.legend()
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.title(r'LV PV loop')
plt.xlabel(r"LV volume (ml)", fontsize=16)
plt.ylabel(r'LV pressure (mmHg)', fontsize=16)
plt.xticks(fontsize=14, rotation=0)
plt.yticks(fontsize=14, rotation=0)
plt.savefig(PVoutfilename)

plt.figure(2)
#plt.legend()
plt.title(r"LV pressure waveform")
plt.xlabel(r"Time (ms)")
plt.ylabel(r"LV pressure (mmHg)")
plt.savefig(LVPoutfilename)

plt.figure(3)
plt.legend()
plt.title("LV volume waveform")
plt.xlabel("Time (ms)")
plt.ylabel("LV volume (ml)")
plt.savefig(LVVoutfilename)

plt.figure(4)
plt.legend()
plt.title("Mitral flow rate")
plt.xlabel("Time (ms)")
plt.ylabel("Flow rate (ml/s)")
plt.savefig(Qmvoutfilename)

plt.figure(5)
plt.legend()
plt.title("Mitral flow velocity")
plt.xlabel("Time (ms)")
plt.ylabel("Velocity (cm/s)")
plt.savefig(Vmvoutfilename)

plt.figure(6)
plt.legend()
plt.title("Coronary flow")
plt.xlabel("Time (ms)")
plt.ylabel("Flow rate (ml/min)")
plt.savefig(Qladoutfilename)

plt.figure(7)
plt.legend()
plt.title("Q distribution")
plt.xlabel("Location")
plt.ylabel("Q")
plt.legend()
plt.savefig(Qdisfilename)

plt.figure(8)
plt.legend(loc='best')
plt.title("Relative Q distribution")
plt.xlabel("Location")
plt.tight_layout()
plt.savefig(RelQdisfilename)

plt.figure(9)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.legend()
plt.xlabel(r"Radial position $R$")
plt.ylabel(r"IMP (mmHg)")
plt.xlim(0,1)
plt.ylim(-10, 1.6*max(homo_LVPP))
plt.tight_layout()
plt.savefig("transmuralvariation_of_imp.png")

plt.figure(10)
plt.legend()
plt.xlim(0,1)
plt.tight_layout()
plt.savefig("transmuralvariation_of_imp_LVP_ratio.png")

plt.figure(11)
plt.legend()
plt.xlabel(r"Radial position")
plt.ylabel(r"$W_{f}$ (mmHg)")
plt.ylim(0,150)
plt.tight_layout()
plt.savefig("WD_distribution.png")

plt.figure(12)
plt.xlabel(r"Radial position")
plt.ylabel(r"$W_{f}/Q$ (mmHg/ml)")
plt.legend()
plt.savefig("QWf_ratio.png")

plt.figure(13)
plt.legend()
plt.xlabel(r"$E_{ff}$")
plt.ylabel(r"$S_{ff}$ (mmHg)")
plt.tight_layout()
plt.savefig("Stress_vs_strain.png")

plt.figure(14)
plt.legend()
#if(isexpt):
#	plt.plot(Hoit2_Ecc_Expt[:,0]*BCL, Hoit2_Ecc_Expt[:,1], 'o', label="Hoit 2", color='k')
#	plt.plot(Gorcsan_Ecc_Expt[:,0]*BCL, Gorcsan_Ecc_Expt[:,1], '*', label="Gorcsan", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{cc}$")
plt.legend()
plt.tight_layout()
plt.savefig("Ecc.png")

plt.figure(15)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{ll}$")
plt.legend()
plt.tight_layout()
plt.savefig("Ell.png")



