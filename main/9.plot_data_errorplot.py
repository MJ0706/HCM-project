import sys, pdb
sys.path.append("/home/lclee/Research/HCM")
import os
import csv
#from heArt.src.postprocessing.postprocessdatalib2 import *
#from postprocessdatalib2 import *
from matplotlib import pylab as plt
from matplotlib import rc
import numpy as np
import math as math
import csv


def GetStiffness_Fn2(alpha, beta, P0, P):

    V = 1.0/beta*math.log((P - P0)/alpha + 1)
    dPdV = alpha*beta*math.exp(beta*V)

    return dPdV

def GetEDPVR_Fn2(alpha, beta, P0, Varray):

    Parray = [P0 + alpha*(math.exp(beta*V) - 1.0) for V in Varray]

    return Parray



def GetStiffness(P, Pm, Vm):

    alpha, beta = GetAlphaBeta(Pm, Vm)
    V = (P/alpha)**(1.0/beta)
    dPdV = alpha*beta*V**(beta - 1.0)

    return dPdV

def GetUnloadedV(Pm, Vm):

    k = 0.6 - 0.006*Pm
    V0 = k*Vm;

    return V0

def GetAlphaBeta(Pm, Vm):

    An = 27.5 #27.78 #27.8+/-0.3
    Bn = 2.81  #2.76+/-0.05
    V0 = GetUnloadedV(Pm, Vm)
    V30 = V0 + (Vm - V0)/(Pm/An)**(1.0/Bn)
    alpha = 30.0/V30**((math.log(Pm/30.0)/(math.log(Vm/V30))))
    beta = math.log(Pm/30.0)/math.log(Vm/V30)

    return alpha, beta


def GetEDPVR(Pm, Vm, Varray):
    alpha, beta = GetAlphaBeta(Pm, Vm)

    V0 = GetUnloadedV(Pm, Vm)
    Parray = [alpha*V**beta  for V in Varray]

    return Parray


def readcsv(filename, ncolstart, ncolend, skip=1, delimiter="\t"):

    reader = csv.reader(open(filename), delimiter=delimiter)
    array = []
    nrow = 0
    for row in reader:
    	if(nrow >= skip):
    		try:
    			array.append([(float(row[p])) for p in range(ncolstart,ncolend+1)])
    		except ValueError:
    			break;
    	nrow += 1
    	
    
    return np.array(array)

def GetUnloadedEDPVR(unloadfilename):

    unloadPV = readcsv(unloadfilename, 0, 2, 1, delimiter=" ")
    ind = np.array(unloadPV[:,0])
    LVP = np.array(unloadPV[:,1])
    LVV = np.array(unloadPV[:,2])

    
    LVP_final = [LVP[k] for k in np.where(ind == int(np.min(ind)))[0]]
    LVV_final = [LVV[k] for k in np.where(ind == int(np.min(ind)))[0]]

    return LVP_final, LVV_final

############################################################################

# Plot effects of contractility cases

rc('text', usetex=False)#True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})


BCLs = [
       	1000,\
	910,\
       	1180,\
       ]
cycle = 10
outdir = './Simulation_without_disarray/'
#os.mkdir(outdir)
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
MV_area = 2.5 
isexpt = True
filenames = [
	     "./clinicaldata/Patient_1_data.txt",\
	     "./clinicaldata/Patient_2_data.txt",\
	     "./clinicaldata/Patient_3_data.txt",\
	     ]

filename = "Patient_1_data.txt"

	
RMSE_V = []
RMSE_P = []
MSE_V = []
MSE_P = []

DBP_model = []
SBP_model = []

Ell_model = []

cnt = 0

for (casename, label, color, filename) in zip(casenames, labels, colors, filenames):

	radialpos    = np.load(casename+".npz")['radialpos']
	longpos      = np.load(casename+".npz")['longpos']
	Sff_t 	     = np.load(casename+".npz")['Sff_t']
	Tff_t 	     = np.load(casename+".npz")['Tff_t']
	Pff_         = np.load(casename+".npz")['Pff_']
        Sff_         = np.load(casename+".npz")['Sff_']
        Tff_         = np.load(casename+".npz")['Tff_']
	WD           = np.load(casename+".npz")['WD']
	Eff 	     = np.load(casename+".npz")['Eff']
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


	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	LVP_array = []
	LVV_array = []
	for row in reader:
		tpt_array.append(float(row[0]))
		LVP_array.append(float(row[2]))
		LVV_array.append(float(row[1]))

 	print type(homo_tptt)
	tpt_expt = (np.array(tpt_array)).astype(int)
	LVP_expt = np.array(LVP_array)
	LVV_expt = np.array(LVV_array)
	
	print len(tpt_expt)
	tptt = (homo_tptt- homo_tptt[0]).astype(int)

	ind = np.zeros(len(tpt_expt)).astype(int)
	LVV = np.zeros(len(tpt_expt))
	LVP = np.zeros(len(tpt_expt))
	print np.shape(tptt), np.shape(tpt_expt)
	print ind

	for i in range(len(tpt_array)):
		ind[i] = np.where(tptt == int(tpt_expt[i]))[0]
		LVV[i] = homo_LVV[ind[i]]
		LVP[i] = homo_LVP[ind[i]]

	print ind
	LVV_error = (LVV_expt-LVV)
	LVP_error = (LVP_expt-LVP)
	
	LVV_error_sq = 0.0
	LVP_error_sq = 0.0
	LVV_expt_sq = 0.0
	LVP_expt_sq = 0.0

	for i in range (len(tpt_array)):
		LVV_error_sq += np.square(LVV_error[i])
		LVP_error_sq += np.square(LVP_error[i])
		LVV_expt_sq  += np.square(LVV_expt[i])
		LVP_expt_sq  += np.square(LVP_expt[i])


	LVV_mean_error = LVV_error_sq/LVV_expt_sq
	LVP_mean_error = LVP_error_sq/LVP_expt_sq
	

	RMSE_V.append(LVV_mean_error)

	RMSE_P.append(LVP_mean_error)


	del tpt_array[:]
	del LVV_array[:]
	del LVP_array[:]


	DBP_model.append(homo_Part[0])
	SBP_model.append(np.max(homo_Part))

	peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))

	Ell_model.append(peakEll)


print DBP_model
print SBP_model



print Ell_model


DBP_expt = [65.0, 66.0, 80.0]
SBP_expt = [126.0, 133.0, 151.0]
Ell_expt = [20.0, 19.0, 13.0]


DBP_Error = np.square(np.subtract(DBP_expt,DBP_model)/DBP_expt).mean()
SBP_Error = np.square(np.subtract(SBP_expt,SBP_model)/SBP_expt).mean()
Ell_Error = np.square(np.subtract(Ell_expt,Ell_model)/Ell_expt)
AO_Error = [DBP_Error, SBP_Error]	
print DBP_Error, SBP_Error, Ell_Error





directory1 = './with_dispersion/P1/k0/simulation/LVelectromechanics/'
unloadfilename = directory1+"BiV_PV.txt"
LVP1, LVV1 = GetUnloadedEDPVR(unloadfilename)

directory2 = './with_dispersion/P2/k0/simulation/LVelectromechanics/'
unloadfilename = directory2+"BiV_PV.txt"
LVP2, LVV2 = GetUnloadedEDPVR(unloadfilename)

directory3 = './with_dispersion/P3/k0/simulation/LVelectromechanics/'
unloadfilename = directory3+"BiV_PV.txt"
LVP3, LVV3 = GetUnloadedEDPVR(unloadfilename)

Pm1 = 8.0
Vm1 = 63.25

Pm2 = 18.0
Vm2 = 82.1 

Pm3 = 14.0
Vm3 = 114.0

V01 = GetUnloadedV(Pm1, Vm1)
P1_Varray = np.linspace(V01,Vm1, num=len(LVV1))
P1_Parray = GetEDPVR(Pm1, Vm1, P1_Varray)
err_P1 = np.concatenate([np.subtract(P1_Parray,LVP1)/P1_Parray, np.subtract(P1_Varray,LVV1)/P1_Varray])
Edpvr_P1 = np.square(err_P1).mean()

V02 = GetUnloadedV(Pm2, Vm2)
P2_Varray = np.linspace(V02,Vm2, num=len(LVV2))
P2_Parray = GetEDPVR(Pm2, Vm2, P2_Varray)
err_P2 = np.concatenate([np.subtract(P2_Parray,LVP2)/P2_Parray, np.subtract(P2_Varray,LVV2)/P2_Varray])
Edpvr_P2 = np.square(err_P2).mean()



V03 = GetUnloadedV(Pm3, Vm3)
P3_Varray = np.linspace(V03,Vm3, num=len(LVV3))
P3_Parray = GetEDPVR(Pm3, Vm3, P3_Varray)
err_P3 = np.concatenate([np.subtract(P3_Parray,LVP3)/P3_Parray, np.subtract(P3_Varray,LVV3)/P3_Varray])
Edpvr_P3 = np.square(err_P3).mean()


edpvr_err =[Edpvr_P1, Edpvr_P2, Edpvr_P3]  


plt.rcParams["font.family"] = "times new roman"
plt.rcParams.update({'font.size': 12})

labels = ['Passive \n EDPVR', 'Volume', 'Strain', 'Blood \nPressure']
x_pos = np.arange(len(labels))
CTEs = [np.linalg.norm(edpvr_err)*100, np.linalg.norm(RMSE_V)*100, np.linalg.norm(Ell_Error)*100, np.linalg.norm(AO_Error)*100]
error = [np.std(edpvr_err)*100, np.std(RMSE_V)*100, np.std(Ell_Error)*100, np.std(AO_Error)*100]

csfont = {'fontname':'Times New Roman',  'fontsize':'12'} 
fig, ax = plt.subplots()
ax.bar(x_pos, CTEs, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel('Error (%)', fontsize = 16)
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, fontsize = 16)
ax.set_title('Error plot',  fontsize = 20)
ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('bar_plot_with_error_bars_3.png')
plt.show()
plt.close()


labels2 = ['Control', 'Non \n Obstructive', 'Obstructive']








