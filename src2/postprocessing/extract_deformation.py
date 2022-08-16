import dolfin as df
import sys
sys.path.append("/mnt/home/fanlei1")
import vtk
import vtk_py as vtk_py
import glob
import numpy as np
import csv
import math
from vtk.util import numpy_support
from heArt.src.postprocessing.postprocessdatalibBiV2 import *

BCL = 800
ncycle = 15
#casename = 'biv_idealized2'
casename = 'HF_fine'
#directory = "./"
directory = "outputs_LVAD/"
homo_directory = directory+casename+"/"
out_directory = "./output_deformation_HF"

# Get indices (ind) corresponding to the time points in ncycle
tpt_array = readtpt(homo_directory + "tpt.txt")
ind = np.where((tpt_array>(ncycle)*BCL)*(tpt_array<(ncycle+1)*BCL))[0]

# Read in mesh
mesh = df.Mesh()
hdf = df.HDF5File(mesh.mpi_comm(), homo_directory + "Data.h5", "r")
hdf.read(mesh,"ME/mesh",False)

# Define field variable for reading in
fieldvariable = "ME/u"

# Count number of nsteps
attr = hdf.attributes(fieldvariable)
nsteps = attr['count'] #Number of data associated with the field variable
fdata = df.File(out_directory + "/u.pvd") 

# for each index in ind, read the displacement in CG2, project it to CG1 and output
for p in ind:
    u = df.Function(df.VectorFunctionSpace(mesh, "CG", 2))
    print(fieldvariable+"/vector_%d"%p)
    hdf.read(u, fieldvariable+"/vector_%d"%p)
    
    uCG1 = df.project(u, df.VectorFunctionSpace(mesh, "CG", 1))
    uCG1.rename("disp", "disp")

    fdata << uCG1

