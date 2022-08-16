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

4. The code can be run within the singularity container. For example, for the code [createLV_refine](./ed_mesh_create/createLV_refine.py)  
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
- [closed loop simulation of LV](./main/Control_HCM_main.py)


### Simulation protocols
