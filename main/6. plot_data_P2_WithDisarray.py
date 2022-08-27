import sys, pdb
sys.path.append("/mnt/home/mojumder/github/HCM")
import os
import csv
from matplotlib import pylab as plt
from matplotlib import rc
import numpy as np

# Plot effects of contractility cases

#rc('text', usetex=False)#True)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})
plt.rcParams["font.family"] = "times new roman"
plt.rcParams.update({'font.size': 12})
BCLs = [
       	910,\
	910,\
       	910,\
	910,\
       	910,\
	910,\
       ]
cycle = 10
outdir = './with_disarray_nonobs_new_01_2/'
os.mkdir(outdir)
casenames  = [

	      './with_dispersion/P2/new/k0/simulation_01/LVelectromechanics/LVelectromechanics',\
	      './with_dispersion/P2/new/k1/simulation_01/LVelectromechanics/LVelectromechanics',\
	      './with_dispersion/P2/new/k2/simulation_01/LVelectromechanics/LVelectromechanics',\
	      './with_dispersion/P2/new/k3/simulation_01/LVelectromechanics/LVelectromechanics',\
	      './with_dispersion/P2/new/k4/simulation_01/LVelectromechanics/LVelectromechanics',\

	     ]

labels  = [
	      r'$\kappa$=0',\
	      r'$\kappa$=0.07',\
	      r'$\kappa$=0.1',\
	      r'$\kappa$=0.14',\
	      r'$\kappa$=0.18',\
	      r'$\kappa$=0.22',\
	      #'ef= 80',\
	      #'ef= 85 endo,-75 epi',\
	     ]

colors=["k",
        "b",
        "r",
	"m",
	"g",
	"y",
	]

PVoutfilename = outdir +"comparePV.png"
LVPoutfilename = outdir +"compareLVP.png"
LVVoutfilename = outdir +"compareLVV.png"
Qmvoutfilename = outdir +"compareQmv.png"
Vmvoutfilename = outdir +"compareVmv.png"
Qladoutfilename =outdir + "compareQlad.png"
IMPoutfilename = outdir +"compareIMP.png"
Qdisfilename = outdir +"compareQdist.png"
IMPdisfilename = outdir +"compareIMPdist.png"
WDoutfilename = outdir +"compareWD.png"
WDdistfilename = outdir +"compareWDdist.png"
RelQdisfilename = outdir +"comparerelQdist.png"

caseno = 0
MV_area = 2.5 # Assuming MV_area = 4cm^2 and Vmax = Q/(A/2) based on poiseulle
#Qtotal_arr = []
isexpt = True


filename3 = "./clinicaldata/Patient_2_data.txt"
reader3 = csv.reader(open(filename3), delimiter=" ")
tpt_array3 = []
LVP_array3 = []
LVV_array3 = []


for row in reader3:
	tpt_array3.append(float(row[0]))
	LVP_array3.append(float(row[2]))
	LVV_array3.append(float(row[1]))

tpt_expt3 = np.array(tpt_array3)
LVP_expt3 = np.array(LVP_array3)
LVV_expt3 = np.array(LVV_array3)
	

for (casename, label, color) in zip(casenames, labels, colors):

	radialpos    = np.load(casename+".npz")['radialpos']
	longpos      = np.load(casename+".npz")['longpos']
	Sff_t 	     = np.load(casename+".npz")['Sff_t']
	Sff_t_2      = np.load(casename+".npz")['Sff_t_2']
	Sff_t_3      = np.load(casename+".npz")['Sff_t_3']
	Sff_t_4      = np.load(casename+".npz")['Sff_t_4']
	Sss_t 	     = np.load(casename+".npz")['Sss_t']
	Snn_t 	     = np.load(casename+".npz")['Snn_t']
	Tff_t 	     = np.load(casename+".npz")['Tff_t']
	Pff_         = np.load(casename+".npz")['Pff_']
	Sff_         = np.load(casename+".npz")['Sff_']
	Tff_         = np.load(casename+".npz")['Tff_']
	WD           = np.load(casename+".npz")['WD']
	WD_2         = np.load(casename+".npz")['WD_2']
	WD_3         = np.load(casename+".npz")['WD_3']
	WD_4         = np.load(casename+".npz")['WD_4']
	WS           = np.load(casename+".npz")['WS']
	WN           = np.load(casename+".npz")['WN']
	Eff 	     = np.load(casename+".npz")['Eff']
	Eff_2 	     = np.load(casename+".npz")['Eff_2']
	Eff_3 	     = np.load(casename+".npz")['Eff_3']
	Eff_4 	     = np.load(casename+".npz")['Eff_4']
	Ess 	     = np.load(casename+".npz")['Ess']
	Enn 	     = np.load(casename+".npz")['Enn']
	Ecc 	     = np.load(casename+".npz")['Ecc']
	Ell 	     = np.load(casename+".npz")['Ell']

	Lff          = np.load(casename+".npz")['Lff']
	Lss          = np.load(casename+".npz")['Lss']
	Lnn          = np.load(casename+".npz")['Lnn']

	Sss_         = np.load(casename+".npz")['Sss_']
	Snn_         = np.load(casename+".npz")['Snn_']
	imp 	     = np.load(casename+".npz")['imp']
	#Qtotal       = np.load(casename+".npz")['Qtotal']
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


	fdatadev = open(outdir+label+"_average_data.txt", "w")

	# Plot PV Loop
	plt.figure(1)
	plt.plot(homo_LVV, homo_LVP, label=label, color=color)

	# Plot Pressure waveform
	plt.figure(2)
	plt.plot(homo_tptt - homo_tptt[0], homo_LVPP, label=label + r"$P_{LV}$", color=color, linestyle='-', linewidth=2)
	#plt.plot(homo_tptt - homo_tptt[0], homo_Pven, label=label + r"$P_{ven}$", color=color, linestyle='--', linewidth=2)
	##plt.plot(homo_tptt - homo_tptt[0], homo_PLA, label=label + r"$P_{LA}$", color=color, linestyle='-.', linewidth=2)
	#plt.plot(homo_tptt - homo_tptt[0], homo_Part, label=label + r"$P_{a,p}$", color=color, linestyle=':', linewidth=2)

	# Plot Volume waveform
	plt.figure(3)
	plt.plot(homo_tptt - homo_tptt[0], homo_LVV , label=label, color=color)



	# Get transmural variation of IMP
	plt.figure(4)
	nbins = 6
	n, __ = np.histogram(radialpos, bins=nbins)
	sy,_ = np.histogram(radialpos, bins=nbins, weights=imp)
	sy2,_ = np.histogram(radialpos, bins=nbins, weights=imp*imp)
	mean_imp = sy / n
	std_imp = np.sqrt(sy2/n - mean_imp*mean_imp)

	#plt.plot(radialpos, imp, 'x', alpha=0.1, color=color)
	#plt.errorbar((__[1:] + __[:-1])/2, mean_imp, yerr=std_imp, color=color, label=label)
	plt.plot((__[:-1])/__[len(__)-2], mean_imp, color=color, label=label)
	#print 'Histogram data'
	print n 
	print sy
	#print max(radialpos)
	#print mean_imp 
	print __
	print _
	
	# Get transmural variation of IMP/LVP
	plt.figure(5)
	sy, _ = np.histogram(radialpos, bins=nbins, weights=imp/max(homo_LVPP))
	sy2, _ = np.histogram(radialpos, bins=nbins, weights=imp*imp/(max(homo_LVPP)*max(homo_LVPP)))
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)

	#plt.plot(radialpos, imp/max(homo_LVPP), 'x', alpha=0.1, color=color)
	#plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, color=color, label=label)
	#plt.errorbar((_[1:] + _[:-1])/2,mean, color=color, label=label)
	plt.plot((__[:-1])/__[len(__)-2], mean, color=color, label=label)
	#plt.ylim(-0.25,1.25)
	plt.ylabel(r"IMP/LVP", fontsize=16)
	plt.xlabel(r"Radial position $R$")

	#if(isexpt):
#		plt.plot(IMP_Expt[:,0], IMP_Expt[:,1], 'o', color=color, label="Expt")

	if(len(casenames)<=1):
		#plt.twinx().errorbar((__[1:] + __[:-1])/2, mean_imp, yerr=std_imp, color=color, label=label)
		#plt.ylim(-0.25*max(homo_LVPP),1.25*max(homo_LVPP))
		plt.ylabel(r'IMP (mmHg)', fontsize=16)


	# Get transmural variation of WD

	plt.figure(6)
	n, _ = np.histogram(radialpos, bins=nbins)
	sy, _ = np.histogram(radialpos, bins=nbins, weights=WD_2)
	sys, _ = np.histogram(radialpos, bins=nbins, weights=WS)
	syn, _ = np.histogram(radialpos, bins=nbins, weights=WN)
	sy2, _ = np.histogram(radialpos, bins=nbins, weights=WD_2*WD_2)
	syt = sy + sys + syn
	mean = sy / n
	means = sys / n
	meann = syn / n
	meant = syt/n
	std = np.sqrt(sy2/n - mean*mean)

	#plt.plot(radialpos, WD, 'x', alpha=0.1, color=color)
	#plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, color=color, label=label)
	plt.plot((__[:-1])/__[len(__)-2], mean, color=color, label=label)

	plt.figure(7)
	plt.plot((__[:-1])/__[len(__)-2], means, color=color, label=label)

	plt.figure(8)
	plt.plot((__[:-1])/__[len(__)-2], meann,color=color, label=label)
	#plt.plot((__[:-1])/__[len(__)-2], meant, linestyle=':',color=color, label=label+"_t")

	print >> fdatadev, '\n\  work_density_fiber', mean[0], mean[-1], np.mean(mean)
	print >> fdatadev, '\n\  work_density_sheet', means[0], means[-1], np.mean(means)
	print >> fdatadev, '\n\  work_density_normal', meann[0], meann[-1], np.mean(meann)

	# Plot Sff vs Eff loop
	plt.figure(9)
	plt.plot(np.mean(Eff_2, axis=1), np.mean(Sff_t_2, axis=1)*0.0075, linestyle='-', color=color, \
		label=label+'_f')


	# Plot Sff vs Eff loop
	
	nbins = 8
	n, _ = np.histogram(radialpos, bins=nbins)
	Sff_bin_array = []
	Eff_bin_array = []
	Sff_err_bin_array = []
	Eff_err_bin_array = []
	
	for p in range(0, len(Sff_t_2[:,1])):
		Sff_av, _ = np.histogram(radialpos, bins=nbins, weights=Sff_t_2[p,:])
		Eff_av, _ = np.histogram(radialpos, bins=nbins, weights=Eff_2[p,:])
	
		Sff_err, _ = np.histogram(radialpos, bins=nbins, weights=Sff_t_2[p,:]*Sff_t_2[p,:])
		Sff_err = np.sqrt(Sff_err/n - Sff_av*Sff_av/(n*n))
	
		Eff_err, _ = np.histogram(radialpos, bins=nbins, weights=Eff_2[p,:]*Eff_2[p,:])
		Eff_err = np.sqrt(Eff_err/n - Eff_av*Eff_av/(n*n))
	
		Sff_bin_array.append(Sff_av/n)
		Eff_bin_array.append(Eff_av/n)
	
		Sff_err_bin_array.append(Sff_err)
		Eff_err_bin_array.append(Eff_err)
	
	
	Sff_bin_array = np.array(Sff_bin_array)
	Eff_bin_array = np.array(Eff_bin_array)
	
	Sff_err_bin_array = np.array(Sff_err_bin_array)
	Eff_err_bin_array = np.array(Eff_err_bin_array)
	plt.figure(10)
	for p in range(0, nbins):
		if(p == 0):
			location = "Endo"
			linestyle= "-"
			plt.plot(Eff_bin_array[:,p], Sff_bin_array[:,p]*0.0075, label=label + location, linestyle=linestyle, color=color)
		elif(p == nbins-1):
			location = "Epi"
			linestyle = "--"
		else:
			location = "Midwall"
			linestyle = "-."
	
		#plt.plot(Eff_bin_array[:,p], Sff_bin_array[:,p]*0.0075, label=label + location, linestyle=linestyle, color=color)

	plt.figure(11)
	for p in range(0, nbins):
		if(p == 0):
			location = "Endo"
			linestyle= "-"
		elif(p == nbins-1):
			location = "Epi"
			linestyle = "--"
			plt.plot(Eff_bin_array[:,p], Sff_bin_array[:,p]*0.0075, label=label + location, linestyle=linestyle, color=color)
		else:
			location = "Midwall"
			linestyle = "-."
	
		#plt.plot(Eff_bin_array[:,p], Sff_bin_array[:,p]*0.0075, label=label + location, linestyle=linestyle, color=color)

	plt.figure(12)
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
	
		#plt.plot(Eff_bin_array[:,p], Sff_bin_array[:,p]*0.0075, label=label + location, linestyle=linestyle, color=color)


	# Get Ecc
	threshold = 0.5

	threshold_rad = 0.8
	plt.figure(14)
	peakEcc = np.max(np.abs(np.mean(Ecc, axis=1)*100))
	plt.plot(tpt - (ncycle)*BCL, np.mean(Ecc, axis=1)*100, linestyle='-', color=color, \
		 label=label+"("+str(int(peakEcc))+"%)",\
		)
	print "Peak Ecc = ", peakEcc



	
	# Get Ell
	plt.figure(15)
	peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))
	plt.plot(tpt - (ncycle)*BCL, np.mean(Ell, axis=1)*100, linestyle='-', color=color, \
		label=label+"("+str(int(peakEll))+"%)",\
		)
	print "Peak Ell = ", peakEll
	print >> fdatadev, '\n\  Peak longitudinal strain', peakEll


	# Get Eff
	plt.figure(16)
	peakEff = np.max(np.abs(np.mean(Eff_2, axis=1)*100))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Eff, axis=1)*100, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakEff))+"%)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Eff_2, axis=1)*100, linestyle='-', color=color, \
		label=label)
	print "Peak Eff = ", peakEff


	# Get Sffactive
	plt.figure(17)
	peakSffact = np.max(np.abs(np.mean(Sff_t_2, axis=1)/1000))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Sffactive, axis=1)/1000, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakSffact))+"kPa)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Sff_t_2, axis=1)/1000, linestyle='-', color=color, \
		label=label)
	print "Peak Sffactive = ", peakSffact
	print >> fdatadev, '\n\  Peak Pk2 total_fiber', peakSffact

   	# Get Lff
	plt.figure(18)
	peakLff = np.max((np.mean(Lff, axis=1)))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Lff, axis=1), linestyle='-', color=color, \
	#	label=label+"("+str(int(peakLff))+")",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Lff, axis=1), linestyle='-', color=color, \
		label=label)
	print "Peak Lff = ", peakLff

	# Get Lss
	plt.figure(19)
	peakLss = np.max((np.mean(Lss, axis=1)))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Lss, axis=1), linestyle='-', color=color, \
	#	label=label+"("+str(int(peakLss))+")",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Lss, axis=1), linestyle='-', color=color, \
		label=label)
	print "Peak Lss = ", peakLss

	# Get Lnn
	plt.figure(20)
	peakLnn = np.max((np.mean(Lnn, axis=1)))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Lnn, axis=1), linestyle='-', color=color, \
	#	label=label+"("+str(int(peakLnn))+")",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Lnn, axis=1), linestyle='-', color=color, \
		label=label)
	print "Peak Lnn = ", peakLnn


	# Get Cauchy sheet stress
	plt.figure(22)
	peakSss_ = np.max(np.abs(np.mean(Sss_, axis=1)/1000))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Sss_, axis=1)/1000, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakSss_))+"kPa)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Sss_, axis=1)/1000, linestyle='-', color=color, \
		label=label)
	print "Peak Sss = ", peakSss_

	# Get Cauchy normal stress
	plt.figure(23)
	peakSnn_ = np.max(np.abs(np.mean(Snn_, axis=1)/1000))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Snn_, axis=1)/1000, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakSnn_))+"kPa)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Snn_, axis=1)/1000, linestyle='-', color=color, \
		label=label)
	print "Peak Snn = ", peakSnn_



	# Get total fiber stress
	plt.figure(24)
	peakSff = np.max(np.abs(np.mean(Sff_t, axis=1)/1000))
	peakTff = np.max(np.abs(np.mean(Tff_t, axis=1)/1000))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Sff_, axis=1)/1000, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakSff_))+"kPa)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Sff_t, axis=1)/1000, linestyle='-', color=color, \
		label=label+"_PK2")
	plt.plot(tpt - (ncycle)*BCL, np.mean(Tff_t, axis=1)/1000, linestyle='--', color=color, \
		label=label+"_Cauchy")
	print "Peak Sff = ", peakSff


	# Get Ess
	plt.figure(25)
	peakEss = np.max(np.abs(np.mean(Ess, axis=1)*100))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Eff, axis=1)*100, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakEff))+"%)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Ess, axis=1)*100, linestyle='-', color=color, \
		label=label)
	print "Peak Ess = ", peakEss


	# Get Enn
	plt.figure(26)
	peakEnn = np.max(np.abs(np.mean(Enn, axis=1)*100))
	#plt.plot(tpt - (ncycle)*BCL, np.mean(Eff, axis=1)*100, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakEff))+"%)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Enn, axis=1)*100, linestyle='-', color=color, \
		label=label)
	print "Peak Enn = ", peakEnn


	# Get total sheet stress
	plt.figure(27)
	peakSss = np.max(np.abs(np.mean(Sss_t, axis=1)/1000))

	#plt.plot(tpt - (ncycle)*BCL, np.mean(Sff_, axis=1)/1000, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakSff_))+"kPa)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Sss_t, axis=1)/1000, linestyle='-', color=color, \
		label=label+"_PK2")

	print "Peak Sss = ", peakSss

	print >> fdatadev, '\n\  Peak Pk2 total_sheet', peakSss
	# Get total normal stress
	plt.figure(28)
	peakSnn = np.max(np.abs(np.mean(Snn_t, axis=1)/1000))

	#plt.plot(tpt - (ncycle)*BCL, np.mean(Sff_, axis=1)/1000, linestyle='-', color=color, \
	#	label=label+"("+str(int(peakSff_))+"kPa)",\
	#	)
	plt.plot(tpt - (ncycle)*BCL, np.mean(Snn_t, axis=1)/1000, linestyle='-', color=color, \
		label=label+"_PK2")

	print "Peak Snn = ", peakSnn

	print >> fdatadev, '\n\  Peak Pk2 total_normal', peakSnn


   	caseno += 1
 
	

	

plt.figure(1)
#plt.legend()
#if(isexpt):
	#plt.plot(LVV_expt, LVP_expt, '*', label="Exp P1", color='k')
	#plt.plot(LVV_expt2, LVP_expt2, '*', label="Exp P2", color='b')
	#plt.plot(LVV_expt3, LVP_expt3, '*', label="Exp P3", color='r')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.title(r'LV PV loop')
plt.xlabel(r"LV volume (ml)", fontsize=16)
plt.ylabel(r'LV pressure (mmHg)', fontsize=16)
plt.xticks(fontsize=14, rotation=0)
plt.yticks(fontsize=14, rotation=0)
plt.savefig( PVoutfilename)

plt.figure(2)
#plt.legend()
#if(isexpt):
	#plt.plot(tpt_expt, LVP_expt, '*', label="Exp P1", color='k')
	#plt.plot(tpt_expt2, LVP_expt2, '*', label="Exp P2", color='b')
	#plt.plot(tpt_expt3, LVP_expt3, '*', label="Exp P3", color='r')
plt.title(r"LV pressure waveform")
plt.xlabel(r"Time (ms)")
plt.ylabel(r"LV pressure (mmHg)")
plt.savefig(LVPoutfilename)

plt.figure(3)
plt.legend()
if(isexpt):
	#plt.plot(tpt_expt, LVV_expt, '*', label="Exp P1", color='k')
	#plt.plot(tpt_expt2, LVV_expt2, '*', label="Exp P2", color='b')
	plt.plot(tpt_expt3, LVV_expt3, '*', label="Exp P3", color='r')
plt.title("LV volume waveform")
plt.xlabel("Time (ms)")
plt.ylabel("LV volume (ml)")
plt.savefig(LVVoutfilename)

plt.figure(4)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.legend()
plt.xlabel(r"Radial position $R$")
plt.ylabel(r"IMP (mmHg)")
plt.xlim(0,1)
#plt.ylim(-10, 1.6*max(homo_LVPP))
plt.tight_layout()
plt.savefig(outdir+"transmuralvariation_of_imp.png")

plt.figure(5)
plt.legend()
plt.xlim(0,1)
plt.ylabel(r"IMP/LVP", fontsize=16)
plt.xlabel(r"Radial position $R$")
plt.tight_layout()
plt.savefig(outdir+"transmuralvariation_of_imp_LVP_ratio.png")

plt.figure(6)
plt.legend()
plt.xlabel(r"Radial position")
plt.ylabel(r"$W_{f}$ (mmHg)")
plt.xlim(0,1)
#plt.ylim(0,150)
plt.tight_layout()
plt.savefig(outdir+"WD_distribution_fiber.png")

plt.figure(7)
plt.legend()
plt.xlabel(r"Radial position")
plt.ylabel(r"$W_{s}$ (mmHg)")
plt.xlim(0,1)
#plt.ylim(0,150)
plt.tight_layout()
plt.savefig(outdir+"WD_distribution_sheet.png")

plt.figure(8)
plt.legend()
plt.xlabel(r"Radial position")
plt.ylabel(r"$W_{n}$ (mmHg)")
plt.xlim(0,1)
#plt.ylim(0,150)
plt.tight_layout()
plt.savefig(outdir+"WD_distribution_normal.png")

plt.figure(9)
plt.legend()
plt.xlabel(r"$E_{ff}$")
plt.ylabel(r"$S_{ff}$ (mmHg)")
plt.tight_layout()
plt.savefig(outdir+"Stress_vs_strain_mean_fiber.png")

plt.figure(10)
plt.legend()
plt.xlabel(r"$E_{ff}$")
plt.ylabel(r"$S_{ff}$ (mmHg)")
plt.tight_layout()
plt.savefig(outdir+"Stress_vs_strain_endo.png")

plt.figure(11)
plt.legend()
plt.xlabel(r"$E_{ff}$")
plt.ylabel(r"$S_{ff}$ (mmHg)")
plt.tight_layout()
plt.savefig(outdir+"Stress_vs_strain_epi.png")


plt.figure(12)
plt.legend()
plt.xlabel(r"$E_{ff}$")
plt.ylabel(r"$S_{ff}$ (mmHg)")
plt.tight_layout()
plt.savefig(outdir+"Stress_vs_strain_mid.png")

plt.figure(14)
plt.legend()
#if(isexpt):
#	plt.plot(Hoit2_Ecc_Expt[:,0]*BCL, Hoit2_Ecc_Expt[:,1], 'o', label="Hoit 2", color='k')
#	plt.plot(Gorcsan_Ecc_Expt[:,0]*BCL, Gorcsan_Ecc_Expt[:,1], '*', label="Gorcsan", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{cc}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Ecc.png")

plt.figure(15)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{ll}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Ell.png")


plt.figure(16)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{ff}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Eff.png")


plt.figure(17)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$S_{ff-active}$ (kPa)")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Sff_total.png")


plt.figure(18)
#plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$Stretch$")
#plt.legend()
plt.tight_layout()
plt.savefig(outdir + "fiber_stretch.png")

plt.figure(19)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$l_{ss}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "sheet_stretch.png")


plt.figure(20)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$l_{nn}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "sheetnormal_stretch.png")

plt.figure(21)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$Stress(kPa)$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "fiber_stress_cauchy_active.png")

plt.figure(22)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$T_{ss}$ (kPa)")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "sheet_stress_cauchy.png")


plt.figure(23)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$T_{nn}$ (kPa)")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "sheetnormal_stress_cauchy.png")


plt.figure(24)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$Stress(kPa)$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "total_stress.png")

plt.figure(25)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{ss}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Ess.png")

plt.figure(26)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{nn}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Enn.png")

plt.figure(27)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$Stress(kPa)$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "total_sheet_stress.png")

plt.figure(28)
plt.legend()
#if(isexpt):
#	plt.plot(Gorcsan_Ell_Expt[:,0]*BCL, Gorcsan_Ell_Expt[:,1], '*', label="Gorcsan", color='k')
#	plt.plot(Smiseth_Ell_Expt[:,0]*BCL, Smiseth_Ell_Expt[:,1], 'o', label="Smiseth", color='k')

plt.xlabel(r"Time (ms)")
plt.ylabel(r"$Stress(kPa)$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "total_normal_stress.png")

