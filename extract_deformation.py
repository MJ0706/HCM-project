import sys, pdb
from dolfin import * 
#sys.path.append("/mnt/home/mojumder/HPCC_HCM")
sys.path.append("/home/lclee/Research/HCM")
#sys.path.append("/mnt/Research")
#sys.path.append("/mnt/Research/JoyheArt")
from src2.sim_protocols.run_BiV_ClosedLoop import run_BiV_ClosedLoop as run_BiV_ClosedLoop
from src2.postprocessing.postprocessdata2_3 import postprocessdata as postprocessdata
from finalcode.ExtractStress import ExtractStress as ExtractStress
from finalcode.extract_work import ExtractWork as ExtractWork
from src2.postprocessing.postprocessdatalib2 import *

'''BCL = 1000
ncycle = 5
casename = 'case_P1'
directory = './with_dispersion/simulation/P1/k0/simulation_7/LVelectromechanics/'

BCL = 910
ncycle = 5
casename = 'case_P2'
directory = './with_dispersion/simulation/P2/k0/simulation_7/LVelectromechanics/'
'''
BCL = 1180
ncycle = 5
casename = 'case_P3'
directory = './with_dispersion/simulation/P3/k0/simulation_15/LVelectromechanics/'


homo_directory = directory
out_directory = homo_directory+ "output_deformation"

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


