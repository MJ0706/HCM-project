#########################################################

import sys
import vtk
import os
from dolfin import *
import pdb
import vtk.util.numpy_support as vtk_support
#sys.path.append("/home/fenics/shared/")
sys.path.append("/mnt/home/mojumder/HPCC_HCM/github/HCM")
import vtk_py as vtk_py
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.stats as stats

############################################################
from finalcode.mean_dev import StimulusDev2

################################################################




def DataExtraction(directory, casename,stimuli):

	filename = directory+casename + ".vtp"
	
	pdata_reader = vtk.vtkXMLPolyDataReader()
	pdata_reader.SetFileName(filename)
	pdata_reader.Update()
	pdata = pdata_reader.GetOutput()

	array = (vtk_support.vtk_to_numpy(pdata.GetPointData().GetArray(stimuli)))

	return array	
	




	




case1 = ['P1' , 'P2', 'P3']




Stimulus = 'fstressPk2_'+'max'
DataName = 'f_22'
for l in range(len(case1)):
	directory1 = './with_dispersion/'+ case1[l]+"/k0/simulation/mechanics/ME/"

	vtufilename = directory1+ 'fstressPk2_'+'max'
	filename = "data_set_"+Stimulus
	fdatadev3 = open(directory1+case1[l]+"_thickness.txt", "w+")
	StimulusDev2(STIMULUS=Stimulus,VTUfilename = vtufilename, fdatadev3 =fdatadev3,  dataname = DataName, DirectorY = directory1, field_type="point")
	

	if l==0:
		P1_array = DataExtraction(directory1, Stimulus,stimuli= "Thickness")
	elif l==1:
		P2_array = DataExtraction(directory1, Stimulus,stimuli= "Thickness")
	else:
		P3_array = DataExtraction(directory1, Stimulus,stimuli= "Thickness")


fig, ax = plt.subplots(figsize = (8,6))	
ax.violinplot([P1_array, P2_array, P3_array][::-1], positions = [3,2,1], showmeans = True)

def set_axis_style2(ax, labels):
	ax.get_yaxis().set_tick_params(labelsize = 14, direction = 'out')
	ax.get_xaxis().set_tick_params(labelsize = 14)
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(np.arange(1, len(labels)+1))
	ax.set_xticklabels(labels, fontsize = 14)
	ax.set_xlim(0.25, len(labels)+0.75)
	ax.set_ylabel('Thickness (cm)', fontsize = 14)
set_axis_style2(ax, ['Control', 'Non-Obstructive','Obstructive'])
#plt.show()	
plt.savefig('thicknessall' +"_1.png")
#plt.show()
plt.close()


	

