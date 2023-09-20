# DPD-Silk-Analysis-Package
This is an updated version of the simulation generation and analysis procedure originally published in "Predicting Silk Fiber Mechanical Properties through Multiscale Simulation and Protein Design" ACS Biomaterials Science and Engineering, 2017

The following instructions document the procedure to generate, simulate, and analyze simulations that were originally presented in Lin et al, Nature Communications 2015, and first released for open use by Rim et al, ACS Biomaterials Science and Engineering. Full citations for these works can be found below. This repository is managed by Sinan Keten’s Computational Nanodynamics Laboratory at Northwestern University.

1.	Lin, S., Ryu, S., Tokareva, O. et al. Predictive modelling-based design and experiments for synthesis and spinning of bioinspired silk fibres. Nat Commun 6, 6892 (2015). https://doi.org/10.1038/ncomms7892

2.	Nae-Gyune Rim, Erin G. Roberts, Davoud Ebrahimi et al. Predicting Silk Fiber Mechanical Properties through Multiscale Simulation and Protein Design, ACS Biomaterials Science & Engineering 2017 3 (8), 1542-1556 [DOI: 10.1021/acsbiomaterials.7b00292](https://doi.org/10.1021/acsbiomaterials.7b00292)

## Notable changes that were made from the original manuscripts:
1.	**Generate_Configuration.m** was modified to include the capability of adding ‘sticky’ terminal ends, introducing a 6th bead type to the system. The file was renamed **Generate_Configuration_Sticky.m**.
2.	The executable, **network_noprint**, made available in supplementary information of Rim et al. was replaced with **network_noprint.py**. This script was written as a direct replacement after running into errors using the original **network_noprint** executable to analyze higher molecular weight chains.
3.	The LAMMPS input file **equil_shear_stretch_sticky.in** was modified from **equil_shear_stretch.in** to run the LAMMPS simulation with terminal ends.

## Software Packages Used:

- LAMMPS: Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) is a molecular dynamics program from Sandia National Laboratories. The LAMMPS software and its installation manual can be found at http://lammps.sandia.gov/. Here we assume software is installed under the Linux environment. We used the LAMMPS, November 17th, 2016 version.

- OCTAVE: Octave is a free open-source programming language primarily intended for numerical computations. It is compatible with MATLAB and can be downloaded from (http://www.gnu.org/software/octave/download.html). In this paper, we assume software is installed under the Linux environment.

- VMD: Visual Molecular Dynamics is a molecular visualization program for displaying, animating, and analyzing large biomolecular systems using 3-D graphics and built-in scripting. VMD supports computers running MacOS X, Unix, or Windows, is distributed free of charge, and includes source code. https://www.ks.uiuc.edu/Research/vmd/

## Software Setup
### LAMMPS Setup:
Download the LAMMPS tarball for the version released on November 17, 2016 from the LAMMPS website (https://download.lammps.org/tars/). Expand the archive in an appropriate location for your new project. Navigate into the directory **lammps-17Nov2016/src** and overwrite **pair_soft.cpp** using the **pair_soft_modified.cpp** filed provided in the directory **DPD-Silk-Analysis-Package/Run_Simulation**. To overwrite, be sure to rename **pair_soft_modified.cpp** to *pair_soft.cpp**. Compile LAMMPS using [make](https://docs.lammps.org/Build_make.html) instructions rather than cmake. Note that the mpi software package is necessary to compile **lmp_mpi**.

### Octave Setup:
Install octave packages per the procedure detailed in http://www.gnu.org/software/octave/download.html.

### VMD Setup:
Install VMD Version 1.9.3 found at https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD.

