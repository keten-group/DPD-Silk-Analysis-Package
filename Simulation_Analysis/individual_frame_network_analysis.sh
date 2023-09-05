#!/bin/bash

module load python/anaconda3.6 # Load python/anaconda3.6 on a Linux machine under bash.
module load matlab # Load matlab on a Linux machine under bash.

for i in {1..7..1} # Change this to loop through the number of frames in equil_{randomseed}.dcd divided by the -stride specified by the ./catdcd command above.
# Do not analyze the initial conformation to avoid writing 'cluster_connectivity.txt' for a very long time.
  do
    cp network_noprint.py get_connectivity.m equil_evolve_$i
    cd equil_evolve_$i
    python network_noprint.py
    matlab -nodisplay -r "get_connectivity ; exit" # get_connectivity.m was originally written to be run using octave. Here, we run it with matlab.
    cd ..
  done

for i in {0..14..1} # Change this to loop through the number of frames in shear_{randomseed}.dcd divided by the -stride specified by the ./catdcd command above.
  do
    cp network_noprint.py get_connectivity.m shear_evolve_$i
    cd shear_evolve_$i
    python network_noprint.py
    matlab -nodisplay -r "get_connectivity ; exit" # get_connectivity.m was originally written to be run using octave. Here, we run it with matlab.
    cd ..
  done

for i in {0..5..1} # Change this to loop through the number of frames in stretch_{randomseed}.dcd divided by the -stride specified by the ./catdcd command above.
  do
    cp network_noprint.py get_connectivity.m stretch_evolve_$i
    cd stretch_evolve_$i
    python network_noprint.py
    matlab -nodisplay -r "get_connectivity ; exit" # get_connectivity.m was originally written to be run using octave. Here, we run it with matlab.
    cd ..
  done
