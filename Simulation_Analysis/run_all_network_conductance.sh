#!/bin/bash
module load matlab
for i in {1..14..1} # Change loop length to reflect the number of frames in equil_{randomseed}.dcd.
  do
    cp Network_Conductance.m equil_evolve_$i
    cd equil_evolve_$i
    matlab -nodisplay -r "Network_Conductance ; exit"
    cd ..
  done

for i in {0..28..1} # Change loop length to reflect the number of frames in shear_{randomseed}.dcd.
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
