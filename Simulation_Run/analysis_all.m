%------------ PART1: run node_bridge_diagram_multiplicity.m to get node-bridge fig for each frame ------
% this will output a 'movie.png' in each subfolder 'xxx_evolve_x'
clear;
delete_temp=1 % control if delete .m copied to each subfolder or not
config=0 % control if plot config or not

if (config==1)
	% equil
	for i=0:7
		checkname1=sprintf('./equil_evolve_%d/boundary_condition.txt', i);
		checkname2=sprintf('./equil_evolve_%d/cluster_connectivity.txt', i);
		checkname3=sprintf('./equil_evolve_%d/cluster_coordinate.txt', i);
		temp=isfile(checkname1) && isfile(checkname2) && isfile(checkname3);
	    if (temp>0)
			copyfolder0=sprintf('./equil_evolve_%d', i);
			copyfile('node_bridge_diagram_multiplicity.m', copyfolder0);
		    dirname=sprintf('./equil_evolve_%d/node_bridge_diagram_multiplicity.m', i);
			run(dirname);
          	if (delete_temp>0)
               delete(dirname);
            end
		end
	end

	% shear
	for i=0:14
		checkname1=sprintf('./shear_evolve_%d/boundary_condition.txt', i);
		checkname2=sprintf('./shear_evolve_%d/cluster_connectivity.txt', i);
		checkname3=sprintf('./shear_evolve_%d/cluster_coordinate.txt', i);
		temp=isfile(checkname1) && isfile(checkname2) && isfile(checkname3);
	    if (temp>0)
			copyfolder0=sprintf('./shear_evolve_%d', i);
			copyfile('node_bridge_diagram_multiplicity.m', copyfolder0);
		    dirname=sprintf('./shear_evolve_%d/node_bridge_diagram_multiplicity.m', i);
			run(dirname);
          	if (delete_temp>0)
               delete(dirname);
            end
		end
	end

	% stretch
	for i=0:5
		checkname1=sprintf('./stretch_evolve_%d/boundary_condition.txt', i);
		checkname2=sprintf('./stretch_evolve_%d/cluster_connectivity.txt', i);
		checkname3=sprintf('./stretch_evolve_%d/cluster_coordinate.txt', i);
		temp=isfile(checkname1) && isfile(checkname2) && isfile(checkname3);
	    if (temp>0)
			copyfolder0=sprintf('./stretch_evolve_%d', i);
			copyfile('node_bridge_diagram_multiplicity.m', copyfolder0);
		    dirname=sprintf('./stretch_evolve_%d/node_bridge_diagram_multiplicity.m', i);
			run(dirname);
          	if (delete_temp>0)
               delete(dirname);
            end
		end
	end
end


%------------PART2: Network_Conductance.m to get conductance and resistance for each frame----------
% this will output a 'conductance_resistance.txt' in each subfolder 'xxx_evolve_x' with format '%f %f'
conduct=0 % control if calc conductance or not

if (conduct==1)
	% equil
	for i=0:7
		checkname1=sprintf('./equil_evolve_%d/boundary_condition.txt', i);
		checkname2=sprintf('./equil_evolve_%d/cluster_connectivity.txt', i);
		checkname3=sprintf('./equil_evolve_%d/cluster_coordinate.txt', i);
		temp=isfile(checkname1) && isfile(checkname2) && isfile(checkname3);
	    if (temp>0)
			copyfolder0=sprintf('./equil_evolve_%d', i);
			copyfile('Network_Conductance.m', copyfolder0);
		    dirname=sprintf('./equil_evolve_%d/Network_Conductance.m', i);
			run(dirname);
          	if (delete_temp>0)
               delete(dirname);
            end
		end
    end
	% shear
	for i=0:14
		checkname1=sprintf('./shear_evolve_%d/boundary_condition.txt', i);
		checkname2=sprintf('./shear_evolve_%d/cluster_connectivity.txt', i);
		checkname3=sprintf('./shear_evolve_%d/cluster_coordinate.txt', i);
		temp=isfile(checkname1) && isfile(checkname2) && isfile(checkname3);
	    if (temp>0)
			copyfolder0=sprintf('./shear_evolve_%d', i);
			copyfile('Network_Conductance.m', copyfolder0);
		    dirname=sprintf('./shear_evolve_%d/Network_Conductance.m', i);
			run(dirname);
          	if (delete_temp>0)
               delete(dirname);
            end
		end
	end
	% stretch
	for i=0:5
		checkname1=sprintf('./stretch_evolve_%d/boundary_condition.txt', i);
		checkname2=sprintf('./stretch_evolve_%d/cluster_connectivity.txt', i);
		checkname3=sprintf('./stretch_evolve_%d/cluster_coordinate.txt', i);
		temp=isfile(checkname1) && isfile(checkname2) && isfile(checkname3);
	    if (temp>0)
			copyfolder0=sprintf('./stretch_evolve_%d', i);
			copyfile('Network_Conductance.m', copyfolder0);
		    dirname=sprintf('./stretch_evolve_%d/Network_Conductance.m', i);
			run(dirname);
          	if (delete_temp>0)
               delete(dirname);
            end
		end
	end
end

%------------PART3: Connectivity_Analysis.m to get network properties evolving as time for all frames in equil and shear----------
% this will output a 'equil.txt', 'shear.txt' and 'stretch.txt' in main folder with 3 arrays in each .txt
% output freq: dcd dump freq x2 (we output pdb at a spacing of 2)
% array1: num_of_frames, Number of Beta Crystals
% array2: num_of_frames, Average Crystal Size
% array3: num_of_frames, Number of Connections
% Note: if necessary .txt (boundary_condition.txt, cluster_connectivity.txt, cluster_coordinate.txt) missing, we write 0 

connect=1 
plot_indi = 1 % control if plot or not for each stage
if (connect>0)
    run('Connectivity_Analysis.m');
end
