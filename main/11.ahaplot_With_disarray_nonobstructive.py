#########################################################

import sys
import vtk
import os
from dolfin import *
import pdb
import vtk.util.numpy_support as vtk_support
sys.path.append("/mnt/home/mojumder/HPCC_HCM/github/HCM")
import vtk_py as vtk_py
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.stats as stats
import csv


############################################################
from finalcode.mean_dev import StimulusDev2
from finalcode.AHA_segment import AHA17segment
from finalcode.bullseye_plot import LVbullseyemulti5 as LVbullseyemulti
################################################################

def aha_strain(directory1, case3, Case1):

	fdata = open( directory+Case1+'_'+case3+"_Ell.txt", "w+")

	#data_strain = [max, mean, min]
	datastrain= ['f_22','f_17', 'f_27']
	if case3 =='max':
		vtufilename = directory1+ 'Ell_'+case3
		DataName = datastrain[0]
		Stimulus = 'Ell_'+case3
	elif case3 =='mean':
		vtufilename = directory1+ 'Ell_'+case3
		DataName = datastrain[1]
		Stimulus = 'Ell_'+case3
	elif case3 =='min':
		vtufilename = directory1+ 'Ell_'+case3
		DataName = datastrain[2]
		Stimulus = 'Ell_'+case3
		
	StimulusDev2(STIMULUS=Stimulus,VTUfilename = vtufilename, fdatadev3 =fdata,  dataname = DataName, DirectorY = directory1, field_type="point")
	
	mean_val = AHA17segment(casename = directory1 + Stimulus,count = 0, stimuli = Stimulus)

	return mean_val



def aha_stress(directory1, case3, Case1):

	fdata = open( directory+Case1+'_'+case3+"_PK2.txt", "w+")

	#data_strain = [max, mean, min]
	datastrain= ['f_22','f_17', 'f_27']
	if case3 =='max':
		vtufilename = directory1+ 'fstressPk2_'+case3
		DataName = datastrain[0]
		Stimulus = 'fstressPk2_'+case3
	elif case3 =='mean':
		vtufilename = directory1+ 'fstressPk2_'+case3
		DataName = datastrain[1]
		Stimulus = 'fstressPk2_'+case3
	elif case3 =='min':
		vtufilename = directory1+ 'fstressPk2_'+case3
		DataName = datastrain[2]
		Stimulus = 'fstressPk2_'+case3
		
	StimulusDev2(STIMULUS=Stimulus,VTUfilename = vtufilename, fdatadev3 =fdata,  dataname = DataName, DirectorY = directory1, field_type="point")
	
	mean_val = AHA17segment(casename = directory1 + Stimulus,count = 0, stimuli = Stimulus)
	print mean_val
	return mean_val

def aha_stress_2(directory1, case3, Case1):

	fdata = open( directory+Case1+'_'+case3+"_PK2.txt", "w+")

	#data_strain = [max, mean, min]
	datastrain= ['f_22','f_17', 'f_27']
	if case3 =='max':
		vtufilename = directory1+ 'fstressPk2_2_'+case3
		DataName = datastrain[0]
		Stimulus = 'fstressPk2_2_'+case3
	elif case3 =='mean':
		vtufilename = directory1+ 'fstressPk2_2_'+case3
		DataName = datastrain[1]
		Stimulus = 'fstressPk2_2_'+case3
	elif case3 =='min':
		vtufilename = directory1+ 'fstressPk2_2_'+case3
		DataName = datastrain[2]
		Stimulus = 'fstressPk2_2_'+case3
		
	StimulusDev2(STIMULUS=Stimulus,VTUfilename = vtufilename, fdatadev3 =fdata,  dataname = DataName, DirectorY = directory1, field_type="point")
	
	mean_val = AHA17segment(casename = directory1 + Stimulus,count = 0, stimuli = Stimulus)
	print mean_val
	return mean_val

def aha_thickness(directory1, case3, Case1):

	fdata = open( directory+Case1+'_'+case3+"_thickness.txt", "w+")

	#data_strain = [max, mean, min]
	datastrain= ['f_22','f_17', 'f_27']
	if case3 =='max':
		vtufilename = directory1+ 'fstressPk2_2_'+case3
		DataName = datastrain[0]
		Stimulus = 'fstressPk2_2_'+case3
	elif case3 =='mean':
		vtufilename = directory1+ 'fstressPk2_2_'+case3
		DataName = datastrain[1]
		Stimulus = 'fstressPk2_2_'+case3
	elif case3 =='min':
		vtufilename = directory1+ 'fstressPk2_2_'+case3
		DataName = datastrain[2]
		Stimulus = 'fstressPk2_2_'+case3
		
	StimulusDev2(STIMULUS=Stimulus,VTUfilename = vtufilename, fdatadev3 =fdata,  dataname = DataName, DirectorY = directory1, field_type="point")
	
	mean_val = AHA17segment(casename = directory1 + Stimulus,count = 0, stimuli = "Thickness")
	print mean_val
	return mean_val


def aha_wd(directory1, case3, Case1, case4):

	fdata = open( directory+Case1+'_'+case3+"_WD.txt", "w+")

	#data_strain = [max, mean, min]
	datastrain= ['f_17']
	
	vtufilename = directory1+ 'WD_' + case4
	DataName = datastrain[0]
	Stimulus = 'WD'
	
		
	StimulusDev2(STIMULUS=Stimulus,VTUfilename = vtufilename, fdatadev3 =fdata,  dataname = DataName, DirectorY = directory1, field_type="point")
	
	mean_val = AHA17segment(casename = directory1 + Stimulus,count = 0, stimuli = Stimulus)
	print mean_val
	return mean_val



plt.rcParams["font.family"] = "times new roman"
plt.rcParams.update({'font.size': 12})

case1 = ['k0' , 'k1', 'k2','k3', 'k4'] 
n = 5; 

aharraystrain = np.zeros((n,17))
aharraystress = np.zeros((n,17))
aharrayworkdensity_f = np.zeros((n,17))
aharrayworkdensity_s = np.zeros((n,17))
aharrayworkdensity_n = np.zeros((n,17))

directory = "./aha_with_disarray_non-obstructive/"
os.mkdir(directory)

for l in range(len(case1)):
				
	#Change the directory for each HCM patients for different dispersion	
	dir1 = './with_dispersion/P2/'+ case1[l]+"/simulation/mechanics/ME/"

	
	AHA_strain = aha_strain(directory1 = dir1, case3 = 'min', Case1 = case1[l])
	AHA_stress = aha_stress(directory1 = dir1, case3 = 'max', Case1 = case1[l])
	AHA_wd_f = aha_wd(directory1 = dir1, case3 = 'mean', Case1 = case1[l], case4 = 'f')


	aharraystrain[l][:] = abs(AHA_strain)*100
	aharraystress[l][:] = (AHA_stress)
	aharrayworkdensity_f[l][:] = (AHA_wd_f)

	

LVbullseyemulti(data = aharraystrain, filename = "Long_strain", Stimulus = "Long_strain", DirectorY = directory, maxm = 30, minm=5)
LVbullseyemulti(data = aharraystress, filename = "Pk2_stress", Stimulus = "Pk2_stress", DirectorY = directory, maxm =75 , minm = 25)
LVbullseyemulti(data = aharrayworkdensity_f, filename = "Work_density_fiber", Stimulus = "Work_density", DirectorY = directory, maxm =20.0 , minm = 0.0)



# field names 
fields = ['1', '2', '3', '4','5', '6', '7', '8','9', '10', '11', '12','13', '14', '15', '16','17'] 




# name of csv file 
filename = directory+ "Aha_segmentation_Strain.csv"
    
# writing to csv file 
with open(filename, 'w') as csvfile: 
    # creating a csv writer object 
    csvwriter = csv.writer(csvfile) 
        
    # writing the fields 
    csvwriter.writerow(fields) 
        
    # writing the data rows 
    csvwriter.writerows(aharraystrain)


# name of csv file 
filename = directory+ "Aha_segmentation_PK2.csv"
    
# writing to csv file 
with open(filename, 'w') as csvfile: 
    # creating a csv writer object 
    csvwriter = csv.writer(csvfile) 
        
    # writing the fields 
    csvwriter.writerow(fields) 
        
    # writing the data rows 
    csvwriter.writerows(aharraystress)



# name of csv file 
filename = directory+ "Aha_segmentation_WD_f.csv"
    
# writing to csv file 
with open(filename, 'w') as csvfile: 
    # creating a csv writer object 
    csvwriter = csv.writer(csvfile) 
        
    # writing the fields 
    csvwriter.writerow(fields) 
        
    # writing the data rows 
    csvwriter.writerows(aharrayworkdensity_f)


