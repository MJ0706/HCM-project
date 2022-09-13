import dolfin as df
import sys
sys.path.append("/home/lclee/Research/HCM")
import vtk
import vtk_py as vtk_py
import glob
import numpy as np
import csv
import math
from vtk.util import numpy_support
from src2.postprocessing.postprocessdatalib2 import *

def ExtractWork(IODet, SimDet, param1, param2, param3):

	dir1 = IODet["outputfolder"] + '/' 
	casename =  IODet["caseID"] 
	BCL = SimDet["HeartBeatLength"]
	ncycle = SimDet["closedloopparam"]["stop_iter"]

	casename2 = 'mechanics'

	directory = dir1 + 'LVelectromechanics/'
	#directory = './outputs_LVelectromechanics_case_P1_60_load/LVelectromechanics/'
	homo_directory = directory
	out_directory = dir1+casename2+'/' 

	# Get indices (ind) corresponding to the time points in ncycle
	tpt_array = readtpt(homo_directory + "tpt.txt")
	ind = np.where((tpt_array>(ncycle)*BCL)*(tpt_array<(ncycle+1)*BCL))[0]

	# Read in mesh
	mesh = df.Mesh()
	hdf = df.HDF5File(mesh.mpi_comm(), homo_directory + "Data.h5", "r")
	hdf.read(mesh,"ME/mesh3",False)

	# Define field variable for reading in
	# Define field variable for reading in
	fieldvariable1 = "ME/"+param1
	fieldvariable2 = "ME/"+param2
	fieldvariable = "ME/WD_"+ param3
	# Count number of nsteps
	attr = hdf.attributes(fieldvariable1)
	nsteps = attr['count'] #Number of data associated with the field variable
	fdata = df.File(out_directory + "/u.pvd") 
	fdatamean = df.File(out_directory + "/"+fieldvariable+ ".pvd") 
	fdatamax = df.File(out_directory + "/"+fieldvariable+ "_max.pvd") 
	fdatamin = df.File(out_directory + "/"+fieldvariable+ "_min.pvd") 

	uDG2 = df.Function(df.FunctionSpace(mesh, "CG", 1))
	uDG2.vector()[:] = 0.0
	uDG2_array = uDG2.vector().array()
	uDGmax = df.Function(df.FunctionSpace(mesh, "DG", 0))
	uDGmax.vector()[:] = 0.0

	uDGmin = df.Function(df.FunctionSpace(mesh, "DG", 0))
	uDGmin.vector()[:] = 100000000.0
	maxVariableArray = uDGmax.vector().array()
	minVariableArray = uDGmin.vector().array()

	cnt = 0.0
	print uDG2.vector().array()

	sff = np.zeros((len(ind), len(uDG2.vector().array())))
	eff = np.zeros((len(ind), len(uDG2.vector().array())))
	print len(ind)
	print len(uDG2.vector().array())
	print ind

	# for each index in ind, read the displacement in CG2, project it to CG1 and output
	for p in ind:
	    u1 = df.Function(df.FunctionSpace(mesh, "CG", 1))
	    print(fieldvariable1+"/vector_%d"%p)
	    hdf.read(u1, fieldvariable1+"/vector_%d"%p)
	    
	    uCG1 = df.project(u1, df.FunctionSpace(mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
	    uCG1.rename("Sff", "Sff")

	    u2 = df.Function(df.FunctionSpace(mesh, "CG", 1))
	    print(fieldvariable2+"/vector_%d"%p)
	    hdf.read(u2, fieldvariable2+"/vector_%d"%p)
	    
	    uCG2 = df.project(u2, df.FunctionSpace(mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
	    uCG2.rename("lff", "lff")

	    sff[p-ind[0],:] = uCG1.vector().array()
	    eff[p-ind[0],:] = uCG2.vector().array()
	
	    #WD = [-1.0*np.trapz(stress[i]*0.0075, 0.5*(stretch[i])**2 -1) for i in range(0, len(stress))]
	    #print WD
	    #uDG2_array = uDG2_array + WD


	    cnt += 1

	print sff
	print len(sff[1,:])

	WD = np.array([-1.0*np.trapz(sff[:,i]*0.0075, eff[:,i]) for i in range(0,len(sff[1,:]))])
	print len(WD)

	uDG2.vector()[:] = WD*0.1323 
	print uDG2.vector().array()
	fdatamean << uDG2






