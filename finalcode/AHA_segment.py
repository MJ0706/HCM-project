#########################################################

import sys
import vtk
import os
from dolfin import *
import pdb
import vtk.util.numpy_support as vtk_support
sys.path.append("/home/lclee/Research/HCM")
import vtk_py as vtk_py
import numpy as np
import math
############################################################

def clipper(domain, C, N, isinsideout):

	plane = vtk.vtkPlane()
	plane.SetOrigin(C)
	plane.SetNormal(N)

	clipper = vtk.vtkClipPolyData()
	clipper.SetClipFunction(plane)
	if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
		clipper.SetInput(domain)
	else:
		clipper.SetInputData(domain)
	clipper.SetInsideOut(isinsideout)
	clipper.Update()
	
	pData = clipper.GetOutput()

	return pData


def SixSegment(domain, count,C, stimuli):


	isinsideout = 1
	Ny1 = [0.0,1.0,0.0]
	Ny2 = [0.0,-1.0,0.0]


	#Anterior section
	domainx1 = clipper(domain, C, Ny1,isinsideout)
	#vtk_py.writeXMLPData(domainx1, "domainx1"+".vtp")

	isinsideout = 1
	N2 = [-np.cos(60), np.sin(60), 0]
	section2 = clipper(domainx1, C, N2,isinsideout)
	#vtk_py.writeXMLPData(section2, "section_"+str(6*count+2)+".vtp")
	isinsideout = 0
	section16 = clipper(domainx1, C, N2,isinsideout)
	#vtk_py.writeXMLPData(section16, "section_"+str(16)+".vtp")
	
	
	isinsideout = 0
	N11 = [-np.cos(60), -np.sin(60), 0]
	section6 = clipper(section16, C, N11,isinsideout)
	#vtk_py.writeXMLPData(section6, "section_"+str(6*count+6)+".vtp")
	isinsideout = 1
	section1 = clipper(section16, C, N11,isinsideout)
	#vtk_py.writeXMLPData(section1, "section_"+str(6*count+1)+".vtp")


	
	#Posterior section
	domainx2 = clipper(domain, C, Ny2,isinsideout)
	#vtk_py.writeXMLPData(domainx2, "domainx2"+".vtp")

	isinsideout = 0
	N3 = [np.cos(60), np.sin(60), 0]
	section3 = clipper(domainx2, C, N3,isinsideout)
	#vtk_py.writeXMLPData(section3, "section_"+str(6*count+3)+".vtp")
	isinsideout = 1
	section45 = clipper(domainx2, C, N3,isinsideout)
	#vtk_py.writeXMLPData(section45, "section_"+str(45)+".vtp")
	
	isinsideout = 0
	N4 = [np.cos(60), -np.sin(60), 0]
	section4 = clipper(section45, C, N4,isinsideout)
	#vtk_py.writeXMLPData(section4, "section_"+str(6*count+4)+".vtp")

	isinsideout = 1
	section5 = clipper(section45, C, N4,isinsideout)
	#vtk_py.writeXMLPData(section5, "section_"+str(6*count+5)+".vtp")


	thickness1 = np.average(vtk_support.vtk_to_numpy(section1.GetPointData().GetArray(stimuli)))
	thickness2 = np.average(vtk_support.vtk_to_numpy(section2.GetPointData().GetArray(stimuli)))
	thickness3 = np.average(vtk_support.vtk_to_numpy(section3.GetPointData().GetArray(stimuli)))
	thickness4 = np.average(vtk_support.vtk_to_numpy(section4.GetPointData().GetArray(stimuli)))
	thickness5 = np.average(vtk_support.vtk_to_numpy(section5.GetPointData().GetArray(stimuli)))
	thickness6 = np.average(vtk_support.vtk_to_numpy(section6.GetPointData().GetArray(stimuli)))

	slice_thickness =np.array([thickness1, thickness2, thickness3, thickness4, thickness5, thickness6])

	

	return slice_thickness


def FourSegment(domain, count,C, stimuli):
	
	isinsideout = 1

	N1 = [-np.cos(45), np.sin(45), 0]
	domainx1 = clipper(domain, C, N1,isinsideout)
	#vtk_py.writeXMLPData(domainx1, "domainx1"+".vtp")
	
	isinsideout = 0
	domainx2 = clipper(domain, C, N1,isinsideout)
	#vtk_py.writeXMLPData(domainx2, "domainx2"+".vtp")
	
	N2 = [-np.cos(45), -np.sin(45), 0]
	isinsideout = 0
	section13 = clipper(domainx1, C, N2,isinsideout)
	#vtk_py.writeXMLPData(section13, "section_"+str(6*count +1)+".vtp")
	
	isinsideout = 1
	section16 = clipper(domainx1, C, N2,isinsideout)
	#vtk_py.writeXMLPData(section16, "section_"+str(6*count +4)+".vtp")
	
	
	isinsideout = 1
	section15 = clipper(domainx2, C, N2,isinsideout)
	#vtk_py.writeXMLPData(section15, "section_"+str(6*count+3)+".vtp")
	
	isinsideout = 0
	section14 = clipper(domainx2, C, N2,isinsideout)
	#vtk_py.writeXMLPData(section14, "section_"+str(6*count+2)+".vtp")

	thickness1 = np.average(vtk_support.vtk_to_numpy(section13.GetPointData().GetArray(stimuli)))
	thickness2 = np.average(vtk_support.vtk_to_numpy(section14.GetPointData().GetArray(stimuli)))
	thickness3 = np.average(vtk_support.vtk_to_numpy(section15.GetPointData().GetArray(stimuli)))
	thickness4 = np.average(vtk_support.vtk_to_numpy(section16.GetPointData().GetArray(stimuli)))


	slice_thickness = np.array([thickness1, thickness2, thickness3, thickness4])



	return slice_thickness

def ahaSegment(domain, count, stimuli, casename):


	fdataset =  open(casename+".txt", "w")

	midx = 0.5*(domain.GetBounds()[0] + domain.GetBounds()[1])
	midy = 0.5*(domain.GetBounds()[2] + domain.GetBounds()[3])
	midz = 0.5*(domain.GetBounds()[4] + domain.GetBounds()[5])
	center = [midx, midy, midz]
	
	print(center,domain.GetBounds()[4], domain.GetBounds()[5] )

	dz = (domain.GetBounds()[5] -domain.GetBounds()[4])/4

	print(dz)
	

	isinsideout = 0
	mean = np.zeros(17)
	domain0 = domain
	for i in range(0,3):

		C = [midx, midy, domain0.GetBounds()[5]-(i+1)*dz]
		Nz = [0.0, 0.0, 1.0]
		domain1 = clipper(domain, C, Nz,isinsideout)
		#vtk_py.writeXMLPData(domain1, "domain_1"+".vtp")


		if i ==0:
			#vtk_py.writeXMLPData(domain1, "domain_1"+".vtp")
			slice1 = SixSegment(domain1, i, C, stimuli)
			mean[0:6]= slice1
			
		elif i ==1:
			#vtk_py.writeXMLPData(domain1, "domain_2"+".vtp")
			slice2 = SixSegment(domain1, i, C, stimuli)
			mean[6:12] = slice2
			
		elif i ==2:
			#vtk_py.writeXMLPData(domain1, "domain_3"+".vtp")
			slice3 = FourSegment(domain1, i, C, stimuli)
			mean[12:16] = slice3

			
			Nz =  [0.0, 0.0, -1.0]
			domain4 = clipper(domain, C, Nz,isinsideout)
			#vtk_py.writeXMLPData(domain4, "section_17.vtp")
			
			slice4 = np.average(vtk_support.vtk_to_numpy(domain4.GetPointData().GetArray(stimuli)))
			mean[-1]= slice4

			avg =  np.average(mean)
			print avg

		
		Nz = [0.0, 0.0, -1.0]
		domain = clipper(domain, C, Nz,isinsideout)
		#vtk_py.writeXMLPData(domain, "domain"+".vtp")
		

	print>>fdataset, count, slice1, slice2, slice3, slice4


	return mean


def AHA17segment(casename, count, stimuli):


	filename = casename + ".vtp"
	
	pdata_reader = vtk.vtkXMLPolyDataReader()
	pdata_reader.SetFileName(filename)
	pdata_reader.Update()
	pdata = pdata_reader.GetOutput()

	domain = pdata
	k = ahaSegment(domain, count, stimuli, casename)

	return k






























