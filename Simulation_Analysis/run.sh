#!/bin/bash

bash run_prep_dcd_pdb_psf.sh

./catdcd -o equil.pdb -otype pdb -stype psf -s ref.psf -stride 2 equil_11111.dcd # Change dcd filename to reflect randomseed.

./catdcd -o shear.pdb -otype pdb -stype psf -s ref.psf -stride 2 shear_11111.dcd # Change dcd filename to reflect randomseed.

./catdcd -o stretch.pdb -otype pdb -stype psf -s ref.psf -stride 2 stretch_11111.dcd # Change dcd filename to reflect randomseed.

module load octave
octave separate_coordinate_single.m

for i in {0..14..1} # Change this to number of frames in equil_{randomseed}.dcd
  do
    cp DFS.m Nbead.txt network_noprint get_coordinate.m get_connectivity.m equil_evolve_$i
    cd equil_evolve_$i
    octave get_coordinate.m
    ./network_noprint
    octave get_connectivity.m
    cd ..
  done

for i in {0..28..1} # Change this to number of frames in shear_{randomseed}.dcd
  do
    cp DFS.m Nbead.txt network_noprint get_coordinate.m get_connectivity.m shear_evolve_$i
    cd shear_evolve_$i
    octave get_coordinate.m
    ./network_noprint
    octave get_connectivity.m
    cd ..
  done

for i in {0..5..1} # Change this to number of frames in stretch_{randomseed}.dcd
  do
    cp DFS.m Nbead.txt network_noprint get_coordinate.m get_connectivity.m stretch_evolve_$i
    cd stretch_evolve_$i
    octave get_coordinate.m
    ./network_noprint
    octave get_connectivity.m
    cd ..
  done

module load matlab
for i in {1..14..1} # Change this loop length to reflect the number of frames in equil_{randomseed}.dcd.
  do
    cp node_bridge_diagram_multiplicity.m equil_evolve_$i
    cd equil_evolve_$i
    matlab -nodisplay -r "node_bridge_diagram_multiplicity ; exit"
    cd ..
  done

for i in {0..28..1} # Change this loop length to reflect the number of frames in shear_{randomseed}.dcd. Remember to subtract 1 for erate 0.0005
  do
    cp node_bridge_diagram_multiplicity.m shear_evolve_$i
    cd shear_evolve_$i
    matlab -nodisplay -r "node_bridge_diagram_multiplicity ; exit" # From: https://arc.umich.edu/software/matlab/#:~:text=To%20run%20a%20MATLAB%20script,m%20from%20the%20current%20directory.&text=Note%20that%20the%20MATLAB%20script,an%20exit%20command%20in%20it.
    cd ..
  done

for i in {0..5..1} # Change this loop length to reflect the number of frames in stretch_{randomseed}.dcd.
  do
    cp node_bridge_diagram_multiplicity.m stretch_evolve_$i
    cd stretch_evolve_$i
    matlab -nodisplay -r "node_bridge_diagram_multiplicity ; exit"
    cd ..
  done

python collect_all_movie.py

module load matlab

matlab -nodisplay -r "Connectivity_Analysis ; exit"

module load matlab
for i in {1..14..1} # Change loop length to reflect the number of frames in equil_{randomseed}.dcd.
  do
    cp Network_Conductance.m equil_evolve_$i
    cd equil_evolve_$i
    matlab -nodisplay -r "Network_Conductance ; exit"
    cd ..
  done

for i in {0..28..1} # Change loop length to reflect the number of frames in shear_{randomseed}.dcd. Remember to subtract 1 for erate 0.0005
  do
    cp Network_Conductance.m shear_evolve_$i
    cd shear_evolve_$i
    matlab -nodisplay -r "Network_Conductance ; exit" # From: https://arc.umich.edu/software/matlab/#:~:text=To%20run%20a%20MATLAB%20script,m%20from%20the%20current%20directory.&text=Note%20that%20the%20MATLAB%20script,an%20exit%20command%20in%20it.
    cd ..
  done

for i in {0..5..1} # Change loop length to reflect the number of frames in stretch_{randomseed}.dcd.
  do
    cp Network_Conductance.m stretch_evolve_$i
    cd stretch_evolve_$i
    matlab -nodisplay -r "Network_Conductance ; exit"
    cd ..
  done