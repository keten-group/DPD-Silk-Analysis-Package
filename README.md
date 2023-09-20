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

## Procedure
See File Description section for details about each file.

### Running a Sample Simulation – Use files in directory ‘Simulation_Run’
1)	Generate initial configurations of randomly distributed peptide chains in a water box in the form of ".psf" and LAMMPS “.data” using **Generate_Configuration_Sticky.m**. See file description for more details.
2)	In the LAMMPS input file, the name of the data file generated in the previous step should be modified (inside the “#file name” section).
3)	In the file **submit.sh**, change the variable _lmp_ to the path to your compiled lammps executable. Change the mpirun command as necessary to specify number of processors (-np) and name of lammps log file (-log). Change or remove the options specified in lines beginning with #SBATCH as necessary. Current #SBATCH options are compatible with [Northwestern Quest Computing Cluster](https://www.it.northwestern.edu/departments/it-services-support/research/computing/quest/). This step will generate the following “dcd” trajectory and stress data files: **equil_11111.dcd, equil_11111_unwrap.dcd, shear_1111.dcd, shear_11111_unwrap.dcd, equil_after_shear_11111.dcd, equil_after_shear_11111_unwrap.dcd, stretch_11111.dcd, stretch_11111_unwrap.dcd, all_stress_11111.txt**.
4)	Once the simulation has completed, open VMD and run the command ‘source view.vmd’ to generate a visually appealing representation of the simulation.

### Analyzing a Sample Simulation – Use files in directory ‘Simulation_Analysis’
1)	Copy the output files from the LAMMPS simulation into the directory **DPD-Silk-Analysis-Package/Simulation_Analysis** (**equil_11111.dcd, equil_11111_unwrap.dcd, shear_11111.dcd, shear_11111_unwrap.dcd, equil_after_shear_11111.dcd, equil_after_shear_11111_unwrap.dcd, stretch_11111.dcd, stretch_11111_unwrap.dcd, all_stress_11111.txt**). 
2)	Copy the psf file without water (the file with the “protein_only.psf” extension) as **ref.psf** into this directory. 
3)	Copy the LAMMPS data file as **ref.data** into this directory.
4)	Open VMD and run the tcl script **make_refpdb.tcl** using the command ‘source make_refpdb.tcl’. This will write a file called ref.pdb, which will be used in the following analysis scripts.
5)	In **individual_frame_coordinate.sh**, change each loop to specify the number of frames to be analyzed from each dcd. This is equal to the number of frames output by each lammps production run (equil, shear, stretch) divided by the -stride specified in the ./catdcd command. Change the variables _Nequil_shot, Nshear_shot, and Nstretch_shot_ in **separate_coordinate_single.m** to match the number of frames to be analyzed. Update the names of each dcd file to reflect the random seed chosen for the simulation.
6)	Run the bash script: **individual_frame_coordinates.sh**. In the first three lines (starting with “catdcd”) from the configurations saved as “dcd” file, -stride 2 means that every other frame is saved into “pdb” format as **coord.pdb** for further analysis. For instance, if there are 14 configurations generated in the equilibrium part of the example, 7 configurations will be saved in the “pdb” format and analyzed inside the first loop.
7)	In **get_connectivity.m** change _MAX_peptide_repeat_ such that it is greater than the maximum number of peptide repeats. 
8)	In **individual_frame_network_analysis.sh** change each loop to specify the number of frames to be analyzed from each dcd as was done for **individual_frame_coordinates.sh**.
9)	Run the bash script: **individual_frame_network_analysis.sh**. 
10)	In **run_all_node_bridge_multiplicity.sh**, change each loop to specify the number of frames to be analyzed from each dcd. This is equal to the number of frames output by each lammps production run (equil, shear, stretch) divided by the -stride specified in the ./catdcd command.
11)	To plot the node-bridge diagram, run the bash script **run_all_node_bridge_multiplicity.sh**.
12)	Generate ".png" files of the node-bridge diagram for each frame by running **collect_all_movie.py**. All the resulting images will be placed in a directory called **movie**.
13)	In **Connectivity_Analysis.m** change the variables _Nequil, Nshear, and Nstretch_ to reflect the number of frames being analyzed. Change the variable _timestep_ to match that specified in the LAMMPS input file. In the ‘% get equil data’ section, change the variable _timestep_ and _dumpfreq_ to match the equilibration run. In the ‘% get shear data’ section, change the _dumpfreq_ and _timestep_ to match the shear run. In the ‘% get stretch data’ section, change the _dumpfreq_ and _timestep_ to match the stretch run.
14)	To run **Connectivity_Analysis.m**, either use MATLAB or run the bash script **run_connectivity_analysis.sh**.
15)	In **run_all_network_conductance.sh** change each loop to specify the number of frames to be analyzed from each dcd. This is equal to the number of frames output by each lammps production run (equil, shear, stretch) divided by the -stride specified in the ./catdcd command.
16)	To plot the network conductance, run the bash script **run_all_network_conductance.sh**.
17)	To conduct a stress strain analysis and plot a stress strain curve, run **ss_analysis.py**.

