import sys, pdb
sys.path.append("/mnt/home/mojumder/HPCC_HCM/github/HCM")
import os
import csv
#from heArt.src.postprocessing.postprocessdatalib2 import *
#from postprocessdatalib2 import *
from matplotlib import pylab as plt
from matplotlib import rc
import numpy as np

# Plot effects of contractility cases

#rc('text', usetex=False)#True)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})
plt.rcParams["font.family"] = "times new roman"
plt.rcParams.update({'font.size': 14})
BCLs = [
       	1000,\
	910,\
       	1180,\
       ]
cycle = 10
outdir = './Simulation_without_disarray/'

os.mkdir(outdir)
casenames  = [

	      './with_dispersion/P1/k0/simulation/LVelectromechanics/LVelectromechanics',\
	      './with_dispersion/P2/k0/simulation/LVelectromechanics/LVelectromechanics',\
	      './with_dispersion/P3/k0/simulation/LVelectromechanics/LVelectromechanics',\

	     ]

labels  = [
	      'Control',\
	      'Non-Obstructive',\
	      'Obstructive',\

	     ]



colors=["k",
        "b",
        "r",
	"m",
	"g",]

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
filename = "./clinicaldata/Patient_1_data.txt"
reader = csv.reader(open(filename), delimiter=" ")
tpt_array = []
LVP_array = []
LVV_array = []

filename2 = "./clinicaldata/Patient_2_data.txt"
reader2 = csv.reader(open(filename2), delimiter=" ")
tpt_array2 = []
LVP_array2 = []
LVV_array2 = []

filename3 = "./clinicaldata/Patient_3_data.txt"
reader3 = csv.reader(open(filename3), delimiter=" ")
tpt_array3 = []
LVP_array3 = []
LVV_array3 = []

for row in reader:
	tpt_array.append(float(row[0]))
	LVP_array.append(float(row[2]))
	LVV_array.append(float(row[1]))

	
tpt_expt = np.array(tpt_array)
LVP_expt = np.array(LVP_array)
LVV_expt = np.array(LVV_array)

for row in reader2:
	tpt_array2.append(float(row[0]))
	LVP_array2.append(float(row[2]))
	LVV_array2.append(float(row[1]))

tpt_expt2 = np.array(tpt_array2)
LVP_expt2 = np.array(LVP_array2)
LVV_expt2 = np.array(LVV_array2)

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


	
	# Plot Sff vs Eff loop
	plt.figure(4)
	plt.plot(np.mean(Eff_2, axis=1), np.mean(Sff_t_2, axis=1)*0.0075, linestyle='-', color=color, \
		label=label)




	
	# Get Ell
	plt.figure(5)
	peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))
	plt.plot(tpt - (ncycle)*BCL, np.mean(Ell, axis=1)*100, linestyle='-', color=color, \
		label=label+"("+str(int(peakEll))+"%)",\
		)
	print "Peak Ell = ", peakEll


   	caseno += 1
 
	

	

plt.figure(1)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.title(r'LV PV loop', fontsize = 20)
plt.xlabel(r"LV volume (ml)", fontsize=16)
plt.ylabel(r'LV pressure (mmHg)', fontsize=16)
plt.tight_layout()
plt.savefig( PVoutfilename)

plt.figure(2)
plt.title(r"LV pressure waveform", fontsize = 20)
plt.xlabel(r"Time (ms)", fontsize = 16)
plt.ylabel(r"LV pressure (mmHg)", fontsize = 16)
plt.tight_layout()
plt.savefig(LVPoutfilename)

plt.figure(3)
plt.legend()
if(isexpt):
	plt.plot(tpt_expt, LVV_expt, '*', label="Exp P1", color='k')
	plt.plot(tpt_expt2, LVV_expt2, '*', label="Exp P2", color='b')
	plt.plot(tpt_expt3, LVV_expt3, '*', label="Exp P3", color='r')
plt.title("LV volume waveform", fontsize = 20)
plt.xlabel("Time (ms)", fontsize = 16)
plt.ylabel("LV volume (ml)", fontsize = 16)
plt.tight_layout()
plt.ylim(0, 1.05*max(LVV_expt3))
plt.savefig(LVVoutfilename)



plt.figure(4)
plt.legend()
plt.xlabel(r"$E_{ff}$")
plt.ylabel(r"$S_{ff}$ (mmHg)")
plt.tight_layout()
plt.savefig(outdir+"Stress_vs_strain_mean_fiber.png")



plt.figure(5)
plt.legend()
plt.xlabel(r"Time (ms)")
plt.ylabel(r"$E_{ll}$")
plt.legend()
plt.tight_layout()
plt.savefig(outdir + "Ell.png")









