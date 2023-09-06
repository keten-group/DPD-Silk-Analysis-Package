#!/bin/bash

./catdcd -o equil.pdb -otype pdb -stype psf -s ref.psf -stride 2 equil_11111.dcd # Change dcd filename to reflect randomseed.

./catdcd -o shear.pdb -otype pdb -stype psf -s ref.psf -stride 2 shear_11111.dcd # Change dcd filename to reflect randomseed.

./catdcd -o stretch.pdb -otype pdb -stype psf -s ref.psf -stride 2 stretch_11111.dcd # Change dcd filename to reflect randomseed.

module load octave # Load octave on a Linux machine under bash.
module load matlab # Load matlab on a Linux machien under bash.

octave separate_coordinate_single.m

for i in {1..7..1} # Change this to loop through the number of frames in equil_{randomseed}.dcd divided by the -stride specified by the ./catdcd command above.
# Do not analyze the initial conformation to avoid writing 'cluster_connectivity.txt' for a very long time.
  do
    cp DFS.m Nbead.txt get_coordinate.m equil_evolve_$i
    cd equil_evolve_$i
    octave get_coordinate.m
    cd ..
  done

for i in {0..14..1} # Change this to loop through the number of frames in shear_{randomseed}.dcd divided by the -stride specified by the ./catdcd command above.
  do
    cp DFS.m Nbead.txt get_coordinate.m shear_evolve_$i
    cd shear_evolve_$i
    octave get_coordinate.m
    cd ..
  done

for i in {0..5..1} # Change this to loop through the number of frames in stretch_{randomseed}.dcd divided by the -stride specified by the ./catdcd command above.
  do
    cp DFS.m Nbead.txt get_coordinate.m stretch_evolve_$i
    cd stretch_evolve_$i
    octave get_coordinate.m
    cd ..
  done
