#!/bin/bash

./catdcd -o equil.pdb -otype pdb -stype psf -s ref.psf -stride 2 equil_11111.dcd # Change dcd filename to reflect randomseed.

./catdcd -o shear.pdb -otype pdb -stype psf -s ref.psf -stride 2 shear_11111.dcd # Change dcd filename to reflect randomseed.

./catdcd -o stretch.pdb -otype pdb -stype psf -s ref.psf -stride 2 stretch_11111.dcd # Change dcd filename to reflect randomseed.

module load octave # Load octave on a Linux machine under bash.
module load python/anaconda3.6 # Load python/anaconda3.6 on a Linux machine under bash.
octave separate_coordinate_single.m

for i in {0..14..1} # Change this to loop through the number of frames in equil_{randomseed}.dcd
  do
    cp DFS.m Nbead.txt network_noprint.py get_coordinate.m get_connectivity.m equil_evolve_$i
    cd equil_evolve_$i
    octave get_coordinate.m
    python network_noprint.py
    octave get_connectivity.m
    cd ..
  done

for i in {0..28..1} # Change this to loop through the number of frames in shear_{randomseed}.dcd
  do
    cp DFS.m Nbead.txt network_noprint.py get_coordinate.m get_connectivity.m shear_evolve_$i
    cd shear_evolve_$i
    octave get_coordinate.m
    python network_noprint.py
    octave get_connectivity.m
    cd ..
  done

for i in {0..5..1} # Change this to loop through the number of frames in stretch_{randomseed}.dcd
  do
    cp DFS.m Nbead.txt network_noprint.py get_coordinate.m get_connectivity.m stretch_evolve_$i
    cd stretch_evolve_$i
    octave get_coordinate.m
    python network_noprint.py
    octave get_connectivity.m
    cd ..
  done
