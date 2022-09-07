from scipy.optimize import curve_fit
import os
import numpy as np
#import matplotlib as mpl
#if os.environ.get('DISPLAY','') == '':
#    print('no display found. Using non-interactive Agg backend')
#    mpl.use('Agg')
import matplotlib.pyplot as plt
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






Pm1 = 8.0
Vm1 = 63.25

Pm2 = 18.0
Vm2 = 82.1 

Pm3 = 14.0
Vm3 = 114.0

V01 = GetUnloadedV(Pm1, Vm1)
P1_Varray = np.linspace(V01,Vm1)
P1_Parray = GetEDPVR(Pm1, Vm1, P1_Varray)


V02 = GetUnloadedV(Pm2, Vm2)
P2_Varray = np.linspace(V02,Vm2)
P2_Parray = GetEDPVR(Pm2, Vm2, P2_Varray)

V03 = GetUnloadedV(Pm3, Vm3)
P3_Varray = np.linspace(V03,Vm3)
P3_Parray = GetEDPVR(Pm3, Vm3, P3_Varray)


directory1 = './with_dispersion/P1/k0/simulation/LVelectromechanics/'
unloadfilename = directory1+"BiV_PV.txt"
LVP1, LVV1 = GetUnloadedEDPVR(unloadfilename)

directory2 = './with_dispersion/P2/k0/simulation/LVelectromechanics/'
unloadfilename = directory2+"BiV_PV.txt"
LVP2, LVV2 = GetUnloadedEDPVR(unloadfilename)

directory3 = './with_dispersion/P3/k0/simulation/LVelectromechanics/'
unloadfilename = directory3+"BiV_PV.txt"
LVP3, LVV3 = GetUnloadedEDPVR(unloadfilename)

label1 = 'Control'
label2 = 'Non Obstructive'
label3 = 'Obstructive'

plt.rcParams["font.family"] = "times new roman"
plt.rcParams.update({'font.size': 12})
plt.figure(1)

plt.plot(P1_Varray, P1_Parray,'o', label="Ktotz_Control", color ='k', alpha = 0.4)
plt.plot(P2_Varray, P2_Parray, 'o', label="Ktotz_Non Obstructive", color ='b', alpha = 0.4)
plt.plot(P3_Varray, P3_Parray, 'o', label="Ktotz_Obstructive", color ='r', alpha = 0.4)

csfont = {'fontname':'Times New Roman', 'fontsize':'16'} 
plt.plot(LVV1, LVP1,  label = label1, color ='k')
plt.plot(LVV2, LVP2,  label = label2, color = 'b')
plt.plot(LVV3, LVP3,  label = label3, color = 'r')
#plt.rcParams.update({'font.size': 12})

plt.xlabel("Volume (ml)", fontsize = 16)
plt.ylabel("Pressure (mmHg)", fontsize = 16)
plt.title("Passive EDPVR", fontsize = 20)
plt.legend()
plt.tight_layout()
plt.savefig("Patient_Klotz.png")
plt.show()




    




