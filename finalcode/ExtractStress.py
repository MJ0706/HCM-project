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

def ExtractStress(IODet, SimDet, param, deno):

	dir1 = IODet["outputfolder"] + '/' 
	casename =  IODet["caseID"] 
	BCL = SimDet["HeartBeatLength"]
	ncycle = SimDet["closedloopparam"]["stop_iter"]

	casename2 = 'mechanics'

	directory = dir1 + 'LVelectromechanics/'

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
	fieldvariable = "ME/" +param 

	# Count number of nsteps
	attr = hdf.attributes(fieldvariable)
	nsteps = attr['count'] #Number of data associated with the field variable
	fdata = df.File(out_directory + "/u.pvd") 
	fdatamean = df.File(out_directory + "/"+fieldvariable+ "_mean.pvd") 
	fdatamax = df.File(out_directory + "/"+fieldvariable+ "_max.pvd") 
	fdatamin = df.File(out_directory + "/"+fieldvariable+ "_min.pvd") 

	uDG2 = df.Function(df.FunctionSpace(mesh, "CG", 1))
	uDG2.vector()[:] = 0.0
	uDG2_array = uDG2.vector().array()
	uDGmax = df.Function(df.FunctionSpace(mesh, "CG", 1))
	uDGmax.vector()[:] = 0.0

	uDGmin = df.Function(df.FunctionSpace(mesh, "CG", 1))
	uDGmin.vector()[:] = 100000000.0
	maxVariableArray = uDGmax.vector().array()
	minVariableArray = uDGmin.vector().array()

	cnt = 0.0
	print uDG2.vector().array()
	# for each index in ind, read the displacement in CG2, project it to CG1 and output
	for p in ind:
	    u = df.Function(df.FunctionSpace(mesh, "CG", 1))
	    print(fieldvariable+"/vector_%d"%p)
	    hdf.read(u, fieldvariable+"/vector_%d"%p)
	    
	    uCG1 = df.project(u, df.FunctionSpace(mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
	    uCG1.rename("lff", "lff")
	    
	    uDG2_array = uDG2_array + uCG1.vector().array()[:] 

	    VariableArray = uCG1.vector().array()

	    for idx, (ls, mls_old) in enumerate(zip(VariableArray, maxVariableArray)): 
		    if ls >= mls_old: 
		        maxVariableArray[idx] = ls

	    for idx, (ls, mls_old) in enumerate(zip(VariableArray, minVariableArray)): 
		    if ls <= mls_old: 
		        minVariableArray[idx] = ls

	    #print uDG2.vector().array()

	    #fdata << uCG1
	    cnt += 1


	uDG2.vector()[:] = uDG2_array/(float(cnt+1) * deno)
	print uDG2.vector().array()
	fdatamean << uDG2


	uDGmax.vector()[:] = maxVariableArray/deno
	print uDGmax.vector().array()
	fdatamax << uDGmax


	uDGmin.vector()[:] = minVariableArray/deno
	print uDGmin.vector().array()
	fdatamin << uDGmin




