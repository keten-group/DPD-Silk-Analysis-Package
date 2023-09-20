# DPD Silk Analysis Package
This is an updated version of the simulation generation and analysis procedure originally published in "Predicting Silk Fiber Mechanical Properties through Multiscale Simulation and Protein Design" ACS Biomaterials Science and Engineering, 2017

The following instructions document the procedure to generate, simulate, and analyze simulations that were originally presented in Lin et al, Nature Communications 2015, and first released for open use by Rim et al, ACS Biomaterials Science and Engineering. Full citations for these works can be found below. This repository is managed by Sinan Keten’s Computational Nanodynamics Laboratory at Northwestern University.

1.	Lin, S., Ryu, S., Tokareva, O. et al. Predictive modelling-based design and experiments for synthesis and spinning of bioinspired silk fibres. Nat Commun 6, 6892 (2015). https://doi.org/10.1038/ncomms7892

2.	Nae-Gyune Rim, Erin G. Roberts, Davoud Ebrahimi et al. Predicting Silk Fiber Mechanical Properties through Multiscale Simulation and Protein Design, ACS Biomaterials Science & Engineering 2017 3 (8), 1542-1556 [DOI: 10.1021/acsbiomaterials.7b00292](https://doi.org/10.1021/acsbiomaterials.7b00292)

## Notable changes that were made from the original manuscripts:
1.	**Generate_Configuration.m** was modified to include the capability of adding ‘sticky’ terminal ends, introducing a 6th bead type to the system. The file was renamed **Generate_Configuration_Sticky.m**.
2.	The executable, **network_noprint**, made available in supplementary information of Rim et al. was replaced with **network_noprint.py**. This script was written as a direct replacement after running into errors using the original **network_noprint** executable to analyze higher molecular weight chains.
3.	The LAMMPS input file **equil_shear_stretch.in** was modified to include parameters for the terminal region sticky ends.

## Software Packages Used:

- LAMMPS: Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) is a molecular dynamics program from Sandia National Laboratories. The LAMMPS software and its installation manual can be found at http://lammps.sandia.gov/. Here we assume software is installed under the Linux environment. We used the LAMMPS, November 17th, 2016 version.

- OCTAVE: Octave is a free open-source programming language primarily intended for numerical computations. It is compatible with MATLAB and can be downloaded from (http://www.gnu.org/software/octave/download.html). In this paper, we assume software is installed under the Linux environment.

- VMD: Visual Molecular Dynamics is a molecular visualization program for displaying, animating, and analyzing large biomolecular systems using 3-D graphics and built-in scripting. VMD supports computers running MacOS X, Unix, or Windows, is distributed free of charge, and includes source code. https://www.ks.uiuc.edu/Research/vmd/

## Software Setup
### LAMMPS Setup:
Download the LAMMPS tarball for the version released on November 17, 2016 from the LAMMPS website (https://download.lammps.org/tars/). Expand the archive in an appropriate location for your new project. Navigate into the directory **lammps-17Nov2016/src** and overwrite **pair_soft.cpp** using the **pair_soft_modified.cpp** filed provided in the directory **DPD-Silk-Analysis-Package/Run_Simulation**. To overwrite, be sure to rename **pair_soft_modified.cpp** to **pair_soft.cpp**. Compile LAMMPS using [make](https://docs.lammps.org/Build_make.html) instructions rather than cmake. Note that the mpi software package is necessary to compile **lmp_mpi**.

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
4)	Open VMD and run the tcl script **make_refpdb.tcl** using the command ‘source **make_refpdb.tcl**’. This will write a file called **ref.pdb**, which will be used in the following analysis scripts.
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

## File Descriptions

- **Generate_Configuration_Sticky.m**  - Writes ‘.psf’ files and a ‘.data’ to be used as input for a LAMMPS simulation. Simulation parameters such as density, protein volume fraction, and silk motif characteristics can be modified. Modify _repeat_motif_spider_ to define motif repeating units. Modify _num_hydrophobic_spider_ to define the number of ‘a’ beads in an A block. Modify _num_hydrophilic_spider_ to define the number of ‘b’ beads in a B block. Modify _num_histidine_spider_ to define how many ‘b’ beads in the H block. Modify _num_sticky_spider_ to define how many beads to place in each terminal region. OUTPUTS: **spider_....psf, spider_..._protein_only.psf, spider_..._evap.psf, spider_....data**

- **pair_soft_modified.cpp** – Code to overwrite **pair_soft.cpp** when compiling lammps. Pair style soft is then used to implement a nonbonded harmonic potential between ‘a’ type beads to represent harmonic bonds.

- **equil_shear_stretch.in** – LAMMPS input file to run equilibration, shear, post-shear equilibration, and stretch simulations which does not include parameters for bead type ‘c’ sticky terminal regions. Data file name in “# file name” section must match the ‘.data’ file output by **Generate_Configuration_Sticky.m**. INPUTS: Requires **spider_....data**. OUTPUTS: **equil_11111.dcd, equil_11111_unwrap.dcd, shear_11111.dcd, shear_11111_unwrap.dcd, equil_after_shear_11111_unwrap.dcd, equil_after_shear_11111_unwrap.dcd, stretch_11111.dcd, stretch_11111_unwrap.dcd, all_stress_11111.txt**.

- **view.vmd** – Generates an appealing visualization for the default simulation. Run this script in VMD using the command ‘source view.vmd’.

- **submit.sh** – Bash commands necessary to run the LAMMPS simulation using mpi.

- **make_refpdb.tcl** – Writes a reference pdb file called **ref.pdb** from **ref.data** which is required by later analysis scripts.

- **individual_frame_coordinates.sh** – Bash script that runs executable **catdcd** and **separate_coordinate_single.m** to generate an individual directory for each frame of the equilibration (equil), shear, and stretch simulation runs. Each directory will have a snapshot pdb file corresponding to that frame called **coord.pdb** which is generated by running **get_coordinate.m** in each directory.

- **individual_frame_network_analysis.sh**  - Bash script that runs **network_noprint.py** and **get_connectivity.m** in each directory.

- **catdcd** – OUTPUTS: **equil.pdb, shear.pdb, stretch.pdb**

- **separate_coordinate_single.m** – OUTPUTS: **Nbead.txt, equil_evolve_*/coord.pdb, shear_evolve_*/coord.pdb, stretch_evolve_*/coord.pdb**.

- **get_coordinate.m**  - OUTPUTS: **boundary_condition.txt, Ndata.txt,  coordinate.txt, coordinates.txt**.

- **network_noprint.py** – OUTPUTS: **cluster_sizes.txt, clusters.txt, clusters_resid.txt**.

- **get_connectivity.m** – OUTPUTS: **connectivity_analysis.txt, connectivity_analysis_numeric.txt, resid_for_connecting_polymer.txt, Ncluster_of_nodes_sizes.txt, cluster_coordinate.txt, link_resid_of_each_cluster.txt, links_per_cluster.txt, connectivity_matrix.txt, cluster_connectivity.txt, cluster_angles.txt, cluster_connectivity.net**.

- **run_all_node_bridge_multiplicity.sh** - Bash script that runs the matlab script **node_bridge_diagram_multiplicity.m** in each of the directories **equil_evolve_*, shear_evolve_*, and stretch_evolve_***.

- **node_bridge_diagram_multiplicity.m** – Draws a diagram showing nodes and brdiges between them. OUTPUTS: movie.png.

- **collect_all_movie.py** – Organizes the output ‘movie.png’ files from ‘node_bridge_diagram_multiplicity.m’ into a single directory called ‘movie.’

- **Connectivity_Analysis.m** – Conducts a connectivity analysis and outputs figures. OUTPUTS: **ave_size_Crystal_equil.png, ave_size_Crystal_shear.png, ave_size_Crystal_stretch.png, num_Beta_equil.png, num_Beta_shear.png, num_Beta_stretch.png, num_Conn_equil.png, num_Conn_shear.png, num_Conn_stretch.png**.

- **run_connectivity_analysis.sh** – Contains the bash commands to run **Connectivity_Analysis.m** in MATLAB without a display.

- **Network_Conductance.m** – Plot the network conductance.

- **run_all_network_conductance.sh** – Bash script that runs the MATLAB script **Network_Conductance.m** in each of the directories **equil_evolve_*, shear_evolve_*, and stretch_evolve_***.

- **ss_analysis.py** – Plot the effective stress versus strain and output as **ss_analysis_out/effective_stress.png**.



