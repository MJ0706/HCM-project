# Hypertrophic cardiomyopathy project

### The mathematical details, data, simulation protocols and results are explained in the manuscript. If you have access to the manuscript, please, follow it first. 

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
mpirun.mpich -np <# processors> python createLV_refine.py
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
  - [Patient_1](./ed_mesh_create/Patient_1/createLV_refine.py) : Control Patient
  - [Patient_2](./ed_mesh_create/Patient_2/createLV_refine.py) : Non-obstructive HCM patient
  - [Patient_3](./ed_mesh_create/Patinet_3/createLV_refine.py) : Obstructive HCM patient
- Create Hdf5 file to run simuations using [create_baselinegeo_animal.py](./ed_mesh_create/create_baselinegeo_animal.py)
  - make sure the directory is correct for specific patient at [create_baselinegeo_animal.py](./ed_mesh_create/create_baselinegeo_animal.py)
- Simulation protocol & post-processing
  - Detail of the code is explained in patient specific code. 
    - Control patient [1.Control.py](./main/1.Control.py)
    - Nonobstructive HCM patient [2. Nonobstructive_main.py](./main/2.Nonobstructive_main.py)
    - Obstructive HCM patient [3. Obstructive_main.py](./main/3.Obstructive_main.py)
    - The simulation cases with various degree of dispersion for non-obstructive HCM patient are: 
        - For kappa = 0.07, [2.Nonobstructive_k1.py](./main/2.Nonobstructive_k1.py) 
        - For kappa = 0.1, [2.Nonobstructive_k2.py](./main/2.Nonobstructive_k2.py)        
        - For kappa = 0.14, [2.Nonobstructive_k3.py](./main/2.Nonobstructive_k3.py)
        - For kappa = 0.18, [2.Nonobstructive_k4.py](./main/2.Nonobstructive_k4.py)
    - The simulation cases with various degree of dispersion for Obstructive HCM patient are: 
        - For kappa = 0.07, [2.Obstructive_k1.py](./main/2.Obstructive_k1.py) 
        - For kappa = 0.1, [2.Obstructive_k2.py](./main/2.Obstructive_k2.py)        
        - For kappa = 0.14, [2.Obstructive_k3.py](./main/2.Obstructive_k3.py)
        - For kappa = 0.18, [2.Obstructive_k4.py](./main/2.Obstructive_k4.py)
        - For kappa = 0.22, [2.Obstructive_k5.py](./main/2.Obstructive_k5.py)
  - Postprocessing of the code is explained at last 4 lines (that are commented) in codes of [13.Postprocessing.py](./main/13.Postprocessing.py). Make sure you run each step at a time using single processor.
  - Klotz plot will be plotted using [4. KlotzPlot.py](/main/4.KlotzPlot.py). Make sure, only passive simulation results are used to plot the Klotz curve using this code. Please, check the references in the manuscript to learn more about Klotz curve. 
  - PV plot for without disarray case will be outlined using [ 5. plot_data_WithoutDisarray.py](./main/5.plot_data_WithoutDisarray.py). Make sure the input directory is correct while running this code. 
  - PV plot for disarray case will be outlined using [6. plot_data_P2_WithDisarray.py](./main/6.plot_data_P2_WithDisarray.py) for non-obstructive patient and [8. plot_data_P3_WithDisarray.py](./main/8.plot_data_P3_WithDisarray.py) for obstructive patient. Make sure the input directory is correct while running this code. 
  - Error bar plot will be plotted by [9. plot_data_errorplot.py](./main/9.plot_data_errorplot.py)
  - AHA plot for without disarray cases can be outlined by [10. ahaplot_WithoutDisarray.py](./main/10.ahaplot_WithoutDisarray.py)
  - AHA plot with disarray cases can be outlined for [Non-obstructive](./main/11.ahaplot_With_disarray_nonobstructive.py) and [Obstructive](./main/11.ahaplot_With_disarray_obstructive.py) patient for various degree of disarray.
  - Deformation can be extracted using [12. extract_deformation.py](/main/12.extract_deformation.py)
  - If you face any issues running with code, email (mojumder@msu.edu) 




