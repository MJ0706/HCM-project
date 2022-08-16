import dolfin as df
import sys
sys.path.append("/mnt/home/lclee/")
import vtk
import vtk_py as vtk_py
import glob
import numpy as np
import csv
import math
from vtk.util import numpy_support

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def extract_PV(filename, BCL, ncycle, isall):

	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	LVP_array = []
	LVV_array = []
	RVP_array = []
	RVV_array = []
	Qmv_array = []
	for row in reader:
		tpt_array.append(float(row[0]))
		LVP_array.append(float(row[1]))
		LVV_array.append(float(row[2]))
		RVP_array.append(float(row[3]))
		RVV_array.append(float(row[4]))
		try:
			Qmv_array.append(float(row[6]))
		except IndexError:
			Qmv_array.append(0)
	
	
	tpt_array = np.array(tpt_array)
	LVP_array = np.array(LVP_array)
	LVV_array = np.array(LVV_array)
	RVP_array = np.array(RVP_array)
	RVV_array = np.array(RVV_array)

	Qmv_array = np.array(Qmv_array)
	
	if(isall):
		ind = np.where(np.logical_and(tpt_array >= 0*BCL, tpt_array <= (ncycle+1)*BCL))
	else:
		ind = np.where(np.logical_and(tpt_array >= ncycle*BCL, tpt_array <= (ncycle+1)*BCL))
	
	return tpt_array[ind], LVP_array[ind], LVV_array[ind], RVP_array[ind], RVV_array[ind], Qmv_array[ind]


def extract_Q(filename, BCL, ncycle):

	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	Qav_array = []
	Qmv_array = []
	Qsa_array = []
	Qsv_array = []
	Qpvv_array = []
	Qtv_array = []
	Qpa_array = []
	Qpv_array = []
	Qlvad_array = []
	for row in reader:
		tpt_array.append(float(row[0]))
		Qav_array.append(float(row[1]))
		Qmv_array.append(float(row[2]))
		Qsa_array.append(float(row[3]))
		Qsv_array.append(float(row[4]))
		Qpvv_array.append(float(row[5]))
		Qtv_array.append(float(row[6]))
		Qpa_array.append(float(row[7]))
		Qpv_array.append(float(row[8]))
		Qlvad_array.append(float(row[11]))
		
	
	tpt_array = np.array(tpt_array)
	Qav_array = np.array(Qav_array)
	Qmv_array = np.array(Qmv_array)
	Qsa_array = np.array(Qsa_array)
	Qsv_array = np.array(Qsv_array)
	Qpvv_array = np.array(Qpvv_array)
	Qtv_array = np.array(Qtv_array)
	Qpa_array = np.array(Qpa_array)
	Qpv_array = np.array(Qpv_array)
	Qlvad_array = np.array(Qlvad_array)
		
	ind = np.where(np.logical_and(tpt_array >= ncycle*BCL, tpt_array <= (ncycle+1)*BCL))
	
	return tpt_array[ind], Qav_array[ind], Qmv_array[ind], Qsa_array[ind], Qsv_array[ind], Qpvv_array[ind], Qtv_array[ind], Qpa_array[ind], Qpv_array[ind], Qlvad_array[ind]


def extract_P(filename, BCL, ncycle):

	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	Psv_array = []
	PLV_array = []
	Psa_array = []
	PLA_array = []
	Ppv_array = []
	PRV_array = []
	Ppa_array = []
	PRA_array = []
	for row in reader:
		tpt_array.append(float(row[0]))
		Psv_array.append(float(row[1]))
		PLV_array.append(float(row[2]))
		Psa_array.append(float(row[3]))
		PLA_array.append(float(row[4]))
		Ppv_array.append(float(row[5]))
		PRV_array.append(float(row[6]))
		Ppa_array.append(float(row[7]))
		PRA_array.append(float(row[8]))
	
	
	tpt_array = np.array(tpt_array)
	Psv_array = np.array(Psv_array)
	PLV_array = np.array(PLV_array)
	Psa_array = np.array(Psa_array)
	PLA_array = np.array(PLA_array)
	Ppv_array = np.array(Ppv_array)
	PRV_array = np.array(PRV_array)
	Ppa_array = np.array(Ppa_array)
	PRA_array = np.array(PRA_array)

	
	ind = np.where(np.logical_and(tpt_array >= ncycle*BCL, tpt_array <= (ncycle+1)*BCL))
	
	return tpt_array[ind], Psv_array[ind], PLV_array[ind], Psa_array[ind], PLA_array[ind], Ppv_array[ind], PRV_array[ind], Ppa_array[ind], PRA_array[ind]


def extract_probe(filename, BCL, ncycle):

	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	array_1 = []
	array_2= []
	array_3 = []
	array_4 = []
	for row in reader:
		tpt_array.append(float(row[0]))
		array_1.append(float(row[1]))
		array_2.append(float(row[2]))
		array_3.append(float(row[3]))
		array_4.append(float(row[4]))
	
	
	tpt_array = np.array(tpt_array)
	array_1 = np.array(array_1)
	array_2 = np.array(array_2)
	array_3 = np.array(array_3)
	array_4 = np.array(array_4)
	
	ind = np.where(np.logical_and(tpt_array >= ncycle*BCL, tpt_array <= (ncycle+1)*BCL))
	
	return tpt_array[ind], [array_1[ind], array_2[ind], array_3[ind], array_4[ind]]



def extractESP(LVP, LVV):

	# Get LVV associated with isovolumic relaxation
	ind = np.where(np.logical_and(LVV >= min(LVV)-0.05, LVV < min(LVV)+0.05))
	
	# Get LVP associated with isovolumic relaxation
	isoLVP = LVP[ind]
	
	# Get ind associated with ES
	ESind = ind[0][isoLVP.argmax(axis=0)]

	return LVP[ESind], LVV[ESind]

def extractEDP(LVP, LVV):

	# Get LVV associated with isovolumic contraction
	ind = np.where(np.logical_and(LVV >= max(LVV)-0.05, LVV < max(LVV)+0.05))
	
	# Get LVP associated with isovolumic contraction
	isoLVP = LVP[ind]
	
	# Get ind associated with ED
	EDind = ind[0][isoLVP.argmin(axis=0)]

	return LVP[EDind], LVV[EDind]




def findESPVR(ESP, ESV):

 	pfit = np.polyfit(ESV, ESP, 1)
	ESPVR = np.poly1d(pfit)

	return pfit, ESPVR




def readcsv(filename, ncolstart, ncolend, skip=1):

	reader = csv.reader(open(filename), delimiter=",")
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

def getradialposition(ugrid, LVendo, RVendo, Epi):

  	points = ugrid.GetPoints()

  	LVendo_ptlocator = vtk.vtkPointLocator()
  	LVendo_ptlocator.SetDataSet(LVendo)
  	LVendo_ptlocator.BuildLocator()

  	RVendo_ptlocator = vtk.vtkPointLocator()
  	RVendo_ptlocator.SetDataSet(RVendo)
  	RVendo_ptlocator.BuildLocator()

  	Epi_ptlocator = vtk.vtkPointLocator()
  	Epi_ptlocator.SetDataSet(Epi)
  	Epi_ptlocator.BuildLocator()

  	radialpos = vtk.vtkFloatArray()
  	radialpos.SetName("radial position")
  	radialpos.SetNumberOfComponents(1)

  	LVendoids = []
  	RVendoids = []
  	Epiids = []

  	for ptid in range(0, ugrid.GetNumberOfPoints()):
  	      point = np.array([points.GetPoint(ptid)[k] for k in range(0,3)])
  	      closestLVendopt = np.array(LVendo.GetPoints().GetPoint(LVendo_ptlocator.FindClosestPoint(point)))
  	      closestRVendopt = np.array(RVendo.GetPoints().GetPoint(RVendo_ptlocator.FindClosestPoint(point)))
  	      closestEpipt = np.array(Epi.GetPoints().GetPoint(Epi_ptlocator.FindClosestPoint(point)))

  	      LVdist = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point, closestLVendopt))
  	      RVdist = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point, closestRVendopt))
  	      Epidist = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point, closestEpipt))

	      if(LVdist < RVdist and Epidist < RVdist):
  	      	wallthickness = vtk.vtkMath.Norm(closestEpipt - closestLVendopt)
  	      	radialpos.InsertNextValue(LVdist/wallthickness)

  	      	if(LVdist < 0.5):
  	      		LVendoids.append(ptid)
  	      	else:
  	      		Epiids.append(ptid)

	      elif(LVdist < RVdist and RVdist < Epidist):
  	      	wallthickness = vtk.vtkMath.Norm(closestRVendopt - closestLVendopt)
  	      	radialpos.InsertNextValue(LVdist/wallthickness)

  	      	if(LVdist < 0.5):
  	      		LVendoids.append(ptid)
  	      	else:
  	      		RVendoids.append(ptid)

	      elif(RVdist < LVdist and LVdist < Epidist):
  	      	wallthickness = vtk.vtkMath.Norm(closestRVendopt - closestLVendopt)
  	      	radialpos.InsertNextValue(LVdist/wallthickness)

  	      	if(LVdist < 0.5):
  	      		LVendoids.append(ptid)
  	      	else:
  	      		RVendoids.append(ptid)


	      else:

	 	wallthickness = vtk.vtkMath.Norm(closestEpipt - closestRVendopt)
  	      	radialpos.InsertNextValue(RVdist/wallthickness)

  	      	if(RVdist < 0.5):
  	      		RVendoids.append(ptid)
  	      	else:
  	      		Epiids.append(ptid)


  	ugrid.GetPointData().AddArray(radialpos)
  	
  	return ugrid, LVendoids, RVendoids, Epiids

def getpointclouds(directory, LVendo, RVendo, Epi, matid):

	mesh = df.Mesh()
	hdf = df.HDF5File(mesh.mpi_comm(),  directory + "/" + "Data.h5", "r")
  	hdf.read(mesh,"ME/mesh",False)
	ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
	hdf.close()

	# Merge the subdomain into 1 unstructuredgrid
	merge = vtk.vtkExtractUnstructuredGrid()
	merge.SetInputData(ugrid)
	merge.MergingOn()
	merge.Update()
	ugrid = merge.GetOutput()


	# Generate point cloud
	cx = 0.5*(ugrid.GetBounds()[0] + ugrid.GetBounds()[1])
	cy = 0.5*(ugrid.GetBounds()[2] + ugrid.GetBounds()[3])
	cz = 0.5*(ugrid.GetBounds()[4] + ugrid.GetBounds()[5])
	bds = max([abs(ugrid.GetBounds()[0] - ugrid.GetBounds()[1]),\
	           abs(ugrid.GetBounds()[2] - ugrid.GetBounds()[3]),\
	           abs(ugrid.GetBounds()[4] - ugrid.GetBounds()[5])])

	ptsource = vtk.vtkPointSource()
	ptsource.SetCenter([cx,cy,cz])
	ptsource.SetRadius(bds/2.0*1.2)
	ptsource.SetNumberOfPoints(200000)
	ptsource.Update()

	selectEnclosed = vtk.vtkSelectEnclosedPoints()
	selectEnclosed.SetInputData(ptsource.GetOutput())
	selectEnclosed.SetSurfaceData(vtk_py.convertUGridtoPdata(ugrid))
	selectEnclosed.SetTolerance(1e-9)
	selectEnclosed.Update()

	thresh = vtk.vtkFloatArray()
  	thresh.SetNumberOfComponents(1);
  	thresh.InsertNextValue(0.5);
  	thresh.InsertNextValue(2.0);
  	thresh.SetName("SelectedPoints");

	selectionNode = vtk.vtkSelectionNode()
    	selectionNode.SetFieldType(1) # POINT
    	selectionNode.SetContentType(7) # INDICES
    	selectionNode.SetSelectionList(thresh) # INDICES
    	selection = vtk.vtkSelection()
    	selection.AddNode(selectionNode)

    	extractSelection = vtk.vtkExtractSelection()
    	extractSelection.SetInputData(0, selectEnclosed.GetOutput())
        extractSelection.SetInputData(1, selection)
    	extractSelection.Update()

	points = extractSelection.GetOutput().GetPoints()

	# Get radial position
	probepointpdata = vtk.vtkPolyData()
	probepointpdata.SetPoints(points)

	ugrid, LVendoids, RVendoids, epiids = getradialposition(ugrid, LVendo, RVendo, Epi)

	radpos = probe(ugrid, probepointpdata)
	radpos_array = [radpos.GetPointData().GetArray("radial position").GetValue(p) for p in range(0, radpos.GetNumberOfPoints())]

	# Get Material IDs
	celllocator = vtk.vtkCellLocator()
	celllocator.SetDataSet(matid)
	celllocator.BuildLocator()

	matid_radpos = vtk.vtkIntArray()	
	matid_radpos.SetName("region_id")	
	for p in range(0, radpos.GetNumberOfPoints()):
		cellid = celllocator.FindCell(radpos.GetPoints().GetPoint(p))
		matid_ =  matid.GetCellData().GetArray("region_id").GetValue(cellid)
		matid_radpos.InsertNextValue(matid_)	

	radpos.GetPointData().AddArray(matid_radpos)
	matid_array = [radpos.GetPointData().GetArray("region_id").GetValue(p) for p in range(0, radpos.GetNumberOfPoints())]

	vtk_py.writeXMLPData(radpos, "radialposition.vtp")
	
	return points, radpos_array, matid_array, radpos
		 
	
def probe(ugrid, pdata):

  	probeFilter = vtk.vtkProbeFilter();
  	probeFilter.SetSourceData(ugrid);
  	if(vtk.vtkVersion.GetVTKMajorVersion <= 5):
  		probeFilter.SetInput(pdata); 
  	else:
  		probeFilter.SetInputData(pdata); 

  	probeFilter.Update()

  	return probeFilter.GetOutput()

def probeqty(directory, fieldvariable, points, ind, index):

	mesh = df.Mesh()
	hdf = df.HDF5File(mesh.mpi_comm(),  directory + "/" + "Data.h5", "r")
  	hdf.read(mesh,"ME/mesh",False)
	ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
	attr = hdf.attributes(fieldvariable)
	nsteps = attr['count']

	var_space = df.FunctionSpace(mesh, "CG", 1)
	var = df.Function(var_space)

	#print nsteps
	#print ind[0][index]
    	dataset = fieldvariable+"/vector_%d"%ind[0][index]
	hdf.read(var, dataset)
	var.rename("var", "var")

	var_vtk = numpy_support.numpy_to_vtk(num_array=var.vector().array()[df.vertex_to_dof_map(var_space)].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
	var_vtk.SetName("var")
	ugrid.GetPointData().AddArray(var_vtk)

	probepointpdata = vtk.vtkPolyData()
	probepointpdata.SetPoints(points)
	npts = points.GetNumberOfPoints()

	probvar = probe(ugrid,probepointpdata).GetPointData().GetArray("var")
	
	point_fieldvararray = [probvar.GetValue(p) for p in range(0, npts)]
	
	hdf.close()

	return np.array(point_fieldvararray)



#def probetimeseries(directory, filebasename, fieldvariable, points, isparallel, ind):
def probetimeseries(directory, fieldvariable, points, ind, elemtype, deg):

	assert ((elemtype == "CG" and deg == 1) or (elemtype == "DG" and deg == 0)),"element type not supported"
	stop
	mesh = df.Mesh()
	hdf = df.HDF5File(mesh.mpi_comm(),  directory + "/" + "Data.h5", "r")
  	hdf.read(mesh,"ME/mesh",False)
	ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
	attr = hdf.attributes(fieldvariable)
	nsteps = attr['count']

	var_space = df.FunctionSpace(mesh, elemtype, deg)
	var = df.Function(var_space)
	
	probepointpdata = vtk.vtkPolyData()
	probepointpdata.SetPoints(points)
	npts = points.GetNumberOfPoints()

	point_fieldvararray = []
	cnt = 1
	for p in ind[0]:

    		dataset = fieldvariable+"/vector_%d"%p
		hdf.read(var, dataset)
		var.rename("var", "var")
		var_vtk = numpy_support.numpy_to_vtk(num_array=var.vector().array()[:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
		var_vtk.SetName("var")
		if(elemtype == "CG"):
			ugrid.GetPointData().AddArray(var_vtk)
		elif(elemtype == "DG"):
			ugrid.GetCellData().AddArray(var_vtk)

		probvar = probe(ugrid,probepointpdata).GetPointData().GetArray("var")
		
		point_fieldvararray.append([probvar.GetValue(p) for p in range(0, npts)])

		cnt += 1
	
	hdf.close()

	return np.array(point_fieldvararray)

def readtpt(filename):

	reader = csv.reader(open(filename), delimiter=" ")
	tpt_array = []
	for row in reader:
		tpt_array.append(int(float(row[0])))

	return np.array(tpt_array)


def setMaterialRegion(directory, filebasename, fieldvariable, isparallel, LVendo, RVendo, Epi):


	filenames = glob.glob(directory + "/" + filebasename + "*")
	filenames.sort()
	filenames = [filename for filename in filenames if filename[-3:] != "pvd"]

	if(filenames[0][-4:] == "pvtu" and isparallel):
		ugrid = vtk_py.readXMLPUGrid(filenames[0])
	
	elif(filenames[0][-4:] == ".vtu" and (not isparallel)):
		ugrid = vtk_py.readXMLUGrid(filenames[0])


	vtk_py.addRegionsToBiV(ugrid, LVendo, RVendo, Epi)

	return ugrid

def setMaterialRegion(directory, LVendo, RVendo, Epi):


	mesh = df.Mesh()
	hdf = df.HDF5File(mesh.mpi_comm(),  directory + "/" + "Data.h5", "r")
  	hdf.read(mesh,"ME/mesh",False)
	ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
	hdf.close()

	vtk_py.addRegionsToBiV(ugrid, LVendo, RVendo, Epi)

	return ugrid



def GetSurfaces(directory, filebasename, fieldvariable, isparallel):

	filenames = glob.glob(directory + "/" + filebasename + "*")
	filenames.sort()
	filenames = [filename for filename in filenames if filename[-3:] != "pvd"]

	if(filenames[0][-4:] == "pvtu" and isparallel):
		ugrid = vtk_py.readXMLPUGrid(filenames[0])
	
	elif(filenames[0][-4:] == ".vtu" and (not isparallel)):
		ugrid = vtk_py.readXMLUGrid(filenames[0])

	Epi = vtk_py.extractUGridBasedOnThreshold(ugrid, fieldvariable, 1)
	Epi  = vtk_py.convertUGridtoPdata(Epi)
	LVendo = vtk_py.extractUGridBasedOnThreshold(ugrid, fieldvariable, 2)
	LVendo  = vtk_py.convertUGridtoPdata(LVendo)
	RVendo = vtk_py.extractUGridBasedOnThreshold(ugrid, fieldvariable, 3)
	RVendo  = vtk_py.convertUGridtoPdata(RVendo)

	return LVendo, RVendo, Epi


def GetDisplacement(directory, variabledata, cycle):

	#assert ((elemtype == "CG" and deg == 1) or (elemtype == "DG" and deg == 0)),"element type not supported"

	mesh = df.Mesh()
	hdf = df.HDF5File(mesh.mpi_comm(),  directory + "/" + "Data.h5", "r")
  	hdf.read(mesh,"ME/mesh",False)
	ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
	#attr = hdf.attributes(fieldvariable)
	#nsteps = attr['count']
	elemtype = "CG" 
	deg = 1
	var_space = df.FunctionSpace(mesh, "CG", 1)
	var = df.Function(var_space)
	datafile = File(directory+"/" +"deformation/u_disp.pvd")
	data = variabledata

	for i in range(0,len(data[:,1])):
         	#pdata = vtk.vtkPolyData()
         	#pdata.DeepCopy(vtkradialpos)
         	#Sff_VTK_data = numpy_support.numpy_to_vtk(num_array=0.0075*Sff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	#Sff_VTK_data.SetName("fstress_")
         	#pdata.GetPointData().AddArray(Sff_VTK_data)
         	#Eff_VTK_data = numpy_support.numpy_to_vtk(num_array=Eff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	#Eff_VTK_data.SetName("Eff_")
         	#pdata.GetPointData().AddArray(Eff_VTK_data)
         	#WD_VTK_data = numpy_support.numpy_to_vtk(num_array=WD.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	#WD_VTK_data.SetName("WD_")
         	#pdata.GetPointData().AddArray(WD_VTK_data)
		var.vector()[:] = disp[i,:]
		#disp_VTK_data = numpy_support.numpy_to_vtk(num_array=disp[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	#disp_VTK_data.SetName("disp_")
         	#pdata.GetPointData().AddArray(disp_VTK_data)
		#print pdata
		datafile << var
		
