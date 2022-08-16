from dolfin import * 
import sys, pdb
sys.path.append("/home/lclee/Research/HCM")

import vtk_py as vtk_py
import numpy as np
from numpy.linalg import norm
import os
import math
from scipy.stats import mode
import collections
import vtk
import vtk.util.numpy_support as vtk_support


def StimulusDev2(STIMULUS, VTUfilename, dataname, fdatadev3, DirectorY, field_type="point"):


    newfilename = VTUfilename + '000000.vtu'

    ugrid1 =vtk_py.readXMLUGrid(newfilename)

    print ugrid1
    
    
    dev1 = vtk_support.vtk_to_numpy(ugrid1.GetPointData().GetArray(dataname))
    print(len(dev1))
    print dev1
    
    pdata1 = vtk_py.convertUGridtoPdata(ugrid1)
    C1 = vtk_py.getcentroid(pdata1)
    ztop1 = pdata1.GetBounds()[5]
    C1 = [C1[0], C1[1], ztop1-0.5]
    clippedheart1 = vtk_py.clipheart(pdata1, C1, [0,0,1], True) 
    epi1 , endo1= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart1)


    cleanendopdata1 = vtk.vtkCleanPolyData()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
	    cleanendopdata1.SetInputData(endo1)
    else:
	    cleanendopdata1.SetInput(endo1)
    cleanendopdata1.Update()
    cleanendo1 = cleanendopdata1.GetOutput()


    cleanepipdata1 = vtk.vtkCleanPolyData()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
	    cleanepipdata1.SetInputData(epi1)
    else:
	    cleanepipdata1.SetInput(epi1)
    cleanepipdata1.Update()
    cleanepi1 = cleanepipdata1.GetOutput()




    


    if(field_type == "point"):

	pointLocator = vtk.vtkPointLocator()
	pointLocator.SetDataSet(cleanepi1)
	pointLocator.BuildLocator()

	deviation1 = vtk.vtkFloatArray()
	deviation1.SetName(STIMULUS)
	deviation1.SetNumberOfComponents(1)


	ptlocators = vtk.vtkPointLocator()
        ptlocators.SetDataSet(ugrid1)
        ptlocators.BuildLocator()

	thickness = vtk.vtkFloatArray()
	thickness.SetName("Thickness")
	thickness.SetNumberOfComponents(1)

	closest_epi_ptid = vtk.vtkIdList()

	# Create a map of points from data with duplicated vertices to one without
	ptmap = []

	for ptid in range(0, cleanendo1.GetNumberOfPoints()):
		endopt = cleanendo1.GetPoints().GetPoint(ptid)
		closest_epi_ptid = pointLocator.FindClosestPoint(endopt)
		closestepipt = cleanepi1.GetPoints().GetPoint(closest_epi_ptid)

		distance1 = math.sqrt(vtk.vtkMath().Distance2BetweenPoints(endopt, closestepipt))
		thickness.InsertNextValue(distance1)

		ptlist = vtk.vtkIdList()
		ptlocators.FindPointsWithinRadius(distance1, endopt, ptlist)

		ptmap.append(ptlist.GetId(0))

		Addition = 0.0 
		for k in range(0, ptlist.GetNumberOfIds()):
			n = ptlist.GetId(k)
			Addition = Addition + dev1[n]
			print n, dev1[n]
			
		mean_dev = Addition/ptlist.GetNumberOfIds()
		print(mean_dev)
		print('Deviation is calculated')
		
		deviation1.InsertNextValue(mean_dev)
		print >> fdatadev3,ptid,mean_dev, distance1


	cleanendo1.GetPointData().AddArray(deviation1)
	cleanendo1.GetPointData().AddArray(thickness)
	print cleanendo1
	

    vtk_py.writeXMLPData(cleanendo1, DirectorY+STIMULUS+".vtp")
    

