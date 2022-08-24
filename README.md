# Hypertrophic cardiomyopathy project

Simulator of mechanics in the heart based on [FEniCS](https://fenicsproject.org/) library.

<!-- TOC -->
  - [Installation and Running the code](#installation-and-running-the-code)
  - [Organization of the code](#organization-of-the-code)
  - [Simulation protocols](#simulation-protocols)

<!-- /TOC -->

### Installation and Running the code
A singularity "build" [file](./SingularitY/Singularity_fenics2017_msucompbiomechlab) is provided that will install all the libraries required to run the code.

1. Install singularity by following the instruction in [here](https://sylabs.io/guides/3.5/admin-guide/installation.html)

2. Build a singularity container using the "build" [file](./SingularitY/Singularity_fenics2017_msucompbiomechlab) with
```
sudo singularity build <container_name>.img Singularity_fenics2017_msucompbiomechlab
```

3. Once the container is built, you can launch the Singularity container by
```
 singularity run <container_name>.img
```

4. The code can be run within the singularity container. For example, for the code [createLV_refine](./ed_mesh_create/Patient_1/createLV_refine.py)  
```
python createLV_refine.py
```
or in parallel
```
mpirun.mpich -np <# processors> createLV_refine.py
```

### Organization of the code
The code is organized as follows:
- [mechanics module](./src2/mechanics)
- [simulation protocols](./src2/sim_protocols/README.md)
- [utilities](./src2/utils)
- [benchmark analytical solution](./src2/bmark_analytical)
- [postprocessing](./src2/postprocessing)

Demo python scripts are also provided to simulate
- Create end diastole mesh file 
  - [Patient_1](./ed_mesh_create/Patient_1/createLV_refine.py)
  - [Patient_2](./ed_mesh_create/Patient_2/createLV_refine.py)
  - [Patient_3](./ed_mesh_create/Patinet_3/createLV_refine.py)
- Create Hdf5 file to run simuations using [create_baselinegeo_animal.py](./ed_mesh_create/create_baselinegeo_animal.py)
  - make sure the directory is correct for specific patient at [create_baselinegeo_animal.py](./ed_mesh_create/create_baselinegeo_animal.py)
- Simulation protocol & post-processing
  - Detail of the code is explained in [Control_HCM_main.py](./main/Control_HCM_main.py)
  - Postprocessing of the code is explained last 4 lines in [Control_HCM_main.py](./main/Control_HCM_main.py)
  - Klotz plot will be plotted using [fitEDPVR_1.py](./main/fitEDPVR_1.py)
  - PV plot for without disarray case will be outlined using [ plot_data_withoutdis.py](./main/plot_data_withoutdis.py)
  - PV plot for disarray case will be outlined using [plot_data_P2.py](./main/plot_data_P2.py)
  - Error bar plot will be plotted by [plot_data_error.py](./main/plot_data_error.py)
  - Deformation can be extracted using [extract_deformation.py](./main/extract_deformation.py)
  - AHA plot can be outlined by [ahaplot.py](./main/ahaplot.py)
  - If you face any issues running with code, email (mojumder@msu.edu) 




