anal_full_process_single_seed.sh                                                                    0000775 0001750 0001750 00000001744 12073616351 016573  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  #!/bin/bash

./catdcd -o equil.pdb -otype pdb -stype psf -s ref.psf -stride 2 equil_11111.dcd

./catdcd -o shear.pdb -otype pdb -stype psf -s ref.psf -stride 3 shear_11111.dcd

./catdcd -o stretch.pdb -otype pdb -stype psf -s ref.psf -stride 1 stretch_11111.dcd

octave separate_coordinate_single.m

#break;

for i in {0..3..1}
  do
    cp DFS.m network_noprint get_coordinate.m get_connectivity.m equil_evolve_$i
    cd equil_evolve_$i
    octave get_coordinate.m
    ./network_noprint
    octave get_connectivity.m
    cd ..
  done

for i in {0..3..1}
  do
    cp DFS.m n45etwork_noprint get_coordinate.m get_connectivity.m shear1_evolve_$i
    cd shear1_evolve_$i
    octave get_coordinate.m
    ./network_noprint
    octave get_connectivity.m
    cd ..
  done

for i in {0..5..1}
  do
    cp DFS.m network_noprint get_coordinate.m get_connectivity.m stretch1_evolve_$i
    cd stretch1_evolve_$i
    octave get_coordinate.m
    ./network_noprint
    octave get_connectivity.m
    cd ..
  done
                            DFS.m                                                                                               0000664 0001750 0001750 00000000376 12073613501 010764  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  function DFS(i,nnmax,nn)
    global n_marker nnode_in_cluster
    if (n_marker(i)==0)
	n_marker(i)=1;
	nnode_in_cluster=nnode_in_cluster+1;
        for j=1:nnmax(i)
	    if (n_marker(nn(j,i))==0)
	        DFS(nn(j,i),nnmax,nn);
	    end
	end
    end
end
                                                                                                                                                                                                                                                                  draw_ss_curve.m                                                                                     0000664 0001750 0001750 00000004202 12043065371 013211  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  totalsteps = 100000;
volume = 96000;
volume_after_evaporation = volume*0.8*0.8*0.8;
maxstrain = 1;

fid=fopen('all_stress_after_equil.txt');
nLines=0;
while (fgets(fid)~=-1),
 nLines = nLines+1;
end
fclose(fid);
fid=fopen('all_stress_after_equil.txt');
fgets(fid);fgets(fid);
strain_data = zeros(nLines-2,1);
stress_data = zeros(nLines-2,6);
for i=1:nLines-2
    line = fgets(fid);
    data = sscanf(line,'%d %f %f %f %f %f %f');
    strain_data(i)=data(1);
    stress_data(i,:)=data(2:7)';
end
fclose(fid);
strain_data = strain_data/totalsteps*maxstrain;
stress_data = stress_data/volume;
stress_in_x = stress_data(:,1)-0.5*(stress_data(:,2)+stress_data(:,3));
figure(1);
plot(strain_data,stress_in_x(:,1),'.-');
ylim([0 0.5]);
title('after equilibration');

fid=fopen('all_stress_after_shear.txt');
nLines=0;
while (fgets(fid)~=-1),
 nLines = nLines+1;
end
fclose(fid);
fid=fopen('all_stress_after_shear.txt');
fgets(fid);fgets(fid);
strain_data = zeros(nLines-2,1);
stress_data = zeros(nLines-2,6);
for i=1:nLines-2
    line = fgets(fid);
    data = sscanf(line,'%d %f %f %f %f %f %f');
    strain_data(i)=data(1);
    stress_data(i,:)=data(2:7)';
end
fclose(fid);
strain_data = strain_data/totalsteps*maxstrain;
stress_data = stress_data/volume;
stress_in_x = stress_data(:,1)-0.5*(stress_data(:,2)+stress_data(:,3));
figure(2);
plot(strain_data,stress_in_x(:,1),'.-');
ylim([0 0.5]);
title('after shear');

fid=fopen('all_stress_after_evaporation.txt');
nLines=0;
while (fgets(fid)~=-1),
 nLines = nLines+1;
end
fclose(fid);
fid=fopen('all_stress_after_evaporation.txt');
fgets(fid);fgets(fid);
strain_data = zeros(nLines-2,1);
stress_data = zeros(nLines-2,6);
for i=1:nLines-2
    line = fgets(fid);
    data = sscanf(line,'%d %f %f %f %f %f %f');
    strain_data(i)=data(1);
    stress_data(i,:)=data(2:7)';
end
fclose(fid);
strain_data = strain_data/totalsteps*maxstrain;
stress_data = stress_data/volume_after_evaporation;
stress_in_x = stress_data(:,1)-0.5*(stress_data(:,2)+stress_data(:,3));
figure(3);
plot(strain_data,stress_in_x(:,1),'.-');
ylim([0 0.5]);
title('after evaporation');
                                                                                                                                                                                                                                                                                                                                                                                              gen_node_bridge.m                                                                                   0000664 0001750 0001750 00000002636 12073613501 013443  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  load('cluster_connectivity.txt');
load('cluster_coordinate.txt');
load('boundary_condition.txt');

bh_x=boundary_condition(1)/2;
bh_y=boundary_condition(2)/2;
bh_z=boundary_condition(3)/2;
num_nodes=length(cluster_coordinate(:,1));
num_bridges=nnz(cluster_connectivity);
num_bond_types=max(max(cluster_connectivity));


fp=fopen('node_bridge.data','w');
fprintf(fp,'node_bridge_model\n');
fprintf(fp,'%d atoms\n',num_nodes);
fprintf(fp,'%d bonds\n',num_bridges);
fprintf(fp,'%d angles\n',0);
fprintf(fp,'%d dihedrals\n',0);
fprintf(fp,'%d impropers\n\n',0);
fprintf(fp,'%d atom types\n',1);
fprintf(fp,'%d bond types\n',num_bond_types);
fprintf(fp,'-%f %f xlo xhi\n',bh_x,bh_x);
fprintf(fp,'-%f %f ylo yhi\n',bh_y,bh_y);
fprintf(fp,'-%f %f zlo zhi\n',bh_z,bh_z);
fprintf(fp,'0.0 0.0 0.0 xy xz yz\n\n');
fprintf(fp,'Masses\n\n');
fprintf(fp,'%d %f\n',1,1);

fprintf(fp,'\nAtoms\n\n');
for i=1:num_nodes
  fprintf(fp,'%d %d %d %f %f %f\n',i,0,1,cluster_coordinate(i,1),cluster_coordinate(i,2),cluster_coordinate(i,3));
end

fprintf(fp,'\nBonds\n\n');

max_multiplicity=0;
k=0;
for i=1:num_nodes-1
  for j=i+1:num_nodes
    if (cluster_connectivity(i,j)>0)
      k=k+1;
      fprintf(fp,'%d %d %d %d\n',k,cluster_connectivity(i,j),i,j);
      if (cluster_connectivity(i,j)>max_multiplicity)
        max_multiplicity = cluster_connectivity(i,j);
      end
    end
  end
end

fclose(fp);

printf('maximum multiplicity=%d\n',max_multiplicity);
                                                                                                  get_connectivity.m                                                                                  0000664 0001750 0001750 00000017635 12073613501 013733  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% connectivity analysis file            %%%
%%% for PDB file with orthogonal cell     %%%
%%% 1st version, 2012-07-31 Seunghwa Ryu  %%%
%%% 2nd version, 2012-08-13 Seunghwa Ryu  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cluster analysis
MAX_peptide_repeat=12; % this number must be larger than maximum peptide repeat
disp('load cluster information');
load('cluster_sizes.txt');
load('clusters.txt');
load('clusters_resid.txt');
load('coordinates.txt');

% connectivity analysis
disp('start connectivity analysis');

cluster_size_MAX=max(cluster_sizes);
Ncluster=length(cluster_sizes);
N_largest_cluster=max(cluster_sizes);
max_resid=max(max(clusters_resid));
connectivity_matrix=zeros(max_resid,MAX_peptide_repeat);
total_connections=0;
total_connection_res=0;
res_connection_matrix=zeros(max_resid,1);
links_per_cluster=zeros(Ncluster,1);
link_resid_of_each_cluster=zeros(Ncluster,100);

for i=1:max_resid
count=0;
    for j=1:Ncluster
	temp=(clusters_resid(j,:)==i);
	resid_exist=sum(temp');
        if (resid_exist>0)
	    count=count+1;
	    connectivity_matrix(i,count)=j;
	end
    end
    connection=length(find(connectivity_matrix(i,:)>0))-1;
    if (connection>0)
	total_connection_res=total_connection_res+1;
        res_connection_matrix(i)=1;
    end
end

MAX_links_per_cluster=0;
links_per_cluster=0;
connecting_polymers=find( (res_connection_matrix'==1));

for i=1:Ncluster
    count=0;
    temp=unique(clusters_resid(i,:));
    temp_index=find(temp>0);
    temp=temp(temp_index);temp_length=length(temp);
    for j=1:temp_length
	for k=1:total_connection_res
	    if (temp(j)==connecting_polymers(k))
		count=count+1;
		link_resid_of_each_cluster(i,count)=temp(j);
	    end
	end
    end
%    links_per_cluster(i)=count;
end

disp('check link between clusters');
% link check
tempa=connectivity_matrix>0;
tempb=sum(tempa');
ind_double=find(tempb>1);
cluster_connectivity=zeros(Ncluster,Ncluster);

for i=1:length(ind_double)
%for i=1:1
    atomid_for_polymer=find(coordinates(:,5)==ind_double(i));
    nodeid_for_polymer=zeros(size(atomid_for_polymer));
    node_search_id=connectivity_matrix(ind_double(i),:);
    node_search_id=node_search_id(find(node_search_id>0));
    for j=1:length(atomid_for_polymer) 
        for k=1:length(node_search_id)
  	    temp=clusters(node_search_id(k),:)==atomid_for_polymer(j);
	    if (sum(temp)>0)
		nodeid_for_polymer(j)=node_search_id(k);
	    end
        end
    end
    for i=1:length(nodeid_for_polymer)-1
	if (nodeid_for_polymer(i)!=nodeid_for_polymer(i+1))
	    id1=nodeid_for_polymer(i);
	    id2=nodeid_for_polymer(i+1);
	    cluster_connectivity(id1,id2)=cluster_connectivity(id1,id2)+1;
	    cluster_connectivity(id2,id1)=cluster_connectivity(id2,id1)+1;
	end
    end
end

links_per_cluster=sum(cluster_connectivity);
MAX_links_per_cluster=max(links_per_cluster);
total_connections=sum(links_per_cluster)/2;
Ncluster_isol=sum(links_per_cluster==0);

nnmax=zeros(Ncluster,1);
nnmax=sum(cluster_connectivity>0);
nn=zeros(MAX_links_per_cluster,Ncluster);
Ncluster_of_nodes=0;
Ncluster_of_nodes_sizes=zeros(Ncluster,1);
for i=1:Ncluster
    count=0;
    for j=1:Ncluster
   	if (cluster_connectivity(i,j)>0)
	    count=count+1;
	    nn(count,i)=j;
	end
    end
end
%size(nn)
%size(nnmax)

global n_marker nnode_in_cluster
n_marker=zeros(Ncluster,1);
for i=1:Ncluster
    nnode_in_cluster=0;
    old_n_marker=n_marker;
    DFS(i,nnmax,nn);
    if (nnode_in_cluster>0)
        Ncluster_of_nodes=Ncluster_of_nodes+1;
        Ncluster_of_nodes_sizes(Ncluster_of_nodes)=nnode_in_cluster;
    end
end

% save data
disp('saving final analysis results');

fp=fopen('connectivity_analysis.txt','w');
fprintf(fp,'Number of clusters = %d\n',Ncluster);
fprintf(fp,'Number of isolated clusters = %d\n',Ncluster_isol);
fprintf(fp,'Max cluster size = %d\n',cluster_size_MAX);
fprintf(fp,'Mean cluster size = %f\n',mean(cluster_sizes));
fprintf(fp,'Median cluster size = %f\n',median(cluster_sizes));
fprintf(fp,'Max links per cluster =%d\n',MAX_links_per_cluster);
fprintf(fp,'Number of connections = %d\n',total_connections);
fprintf(fp,'\n');
fprintf(fp,'Number of clusters of nodes = %d\n',Ncluster_of_nodes);
fprintf(fp,'MAX number of clusters of nodes =%d\n',max(Ncluster_of_nodes_sizes));
fprintf(fp,'\n');
fprintf(fp,'Number of peptides connecting crystals = %d\n',total_connection_res);
fprintf(fp,'Total number of peptides =%d\n',max_resid);
fclose(fp);

fp=fopen('connectivity_analysis_numeric.txt','w');
fprintf(fp,'%d\n',Ncluster);
fprintf(fp,'%d\n',Ncluster_isol);
fprintf(fp,'%d\n',cluster_size_MAX);
fprintf(fp,'%f\n',mean(cluster_sizes));
fprintf(fp,'%f\n',median(cluster_sizes));
fprintf(fp,'%d\n',MAX_links_per_cluster);
fprintf(fp,'%d\n',total_connections);
fprintf(fp,'%d\n',Ncluster_of_nodes);
fprintf(fp,'%d\n',max(Ncluster_of_nodes_sizes));
fprintf(fp,'%d\n',total_connection_res);
fprintf(fp,'%d\n',max_resid);
fclose(fp);

fp=fopen('resid_for_connecting_polymer.txt','w');
fprintf(fp,'%d ',connecting_polymers);
fclose(fp);

fp=fopen('Ncluster_of_nodes_sizes.txt','w');
for i=1:Ncluster_of_nodes
    fprintf(fp,'%d \n',Ncluster_of_nodes_sizes(i));
end
fclose(fp);

fp=fopen('cluster_sizes.txt','w');
for i=1:Ncluster
    fprintf(fp,'%d\n',cluster_sizes(i));
end
fclose(fp);

fp=fopen('clusters.txt','w');
for i=1:Ncluster
    for j=1:cluster_size_MAX
	fprintf(fp,'%d ',clusters(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);

fid = fopen('coordinate.txt','r');
line = fgets(fid);
data = sscanf(line,'%f %f %f %d');
Hx=data(1);Hy=data(2);Hz=data(3);
fclose(fid);

fp=fopen('cluster_coordinate.txt','w');
for i=1:Ncluster
    temp=[ 0 0 0];
    temp0=coordinates(clusters(i,1),1:3);
    for j=1:cluster_size_MAX
	if (clusters(i,j)>0)
            rj=coordinates(clusters(i,j),1:3);
            dr=(rj-temp0);
            dr=dr-round(dr);
	    temp=temp+dr;
	end
    end
    temp=temp/cluster_sizes(i);
    temp=temp+temp0;
    temp=temp-round(temp);
    temp(1)=temp(1)*Hx;temp(2)=temp(2)*Hy;temp(3)=temp(3)*Hz;
    fprintf(fp,'%f %f %f \n',temp(1),temp(2),temp(3));
end
fclose(fp);

fp=fopen('clusters_resid.txt','w');
for i=1:Ncluster
    for j=1:cluster_size_MAX
        fprintf(fp,'%d ',clusters_resid(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);

fp=fopen('link_resid_of_each_cluster.txt','w');
for i=1:Ncluster
    for j=1:MAX_links_per_cluster
	fprintf(fp,'%d ',link_resid_of_each_cluster(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);

fp=fopen('links_per_cluster.txt','w');
for i=1:Ncluster
    fprintf(fp,'%d\n',links_per_cluster(i));
end

fp=fopen('connectivity_matrix.txt','w');
for i=1:max_resid
    for j=1:MAX_peptide_repeat
	fprintf(fp,'%d ',connectivity_matrix(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);

fp=fopen('cluster_connectivity.txt','w');
for i=1:Ncluster
    for j=1:Ncluster
	fprintf(fp,'%d ',cluster_connectivity(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);

Nedges=nnz(cluster_connectivity)/2;
edge_angles=zeros(Nedges,5);

load('cluster_coordinate.txt');
count=0;
for i=1:Ncluster
    for j=i+1:Ncluster
	if (cluster_connectivity(i,j)>0)
	    count=count+1;
	    edge_angles(count,1)=i;
	    edge_angles(count,2)=j;
	    edge_angles(count,3)=cluster_connectivity(i,j);
	    drij=cluster_coordinate(i,:)-cluster_coordinate(j,:);
	    drij=drij./[Hx Hy Hz];
	    drij=drij-round(drij);
	    drij=drij.*[Hx Hy Hz];
	    drij_n=norm(drij);
	    edge_angles(count,4)=drij_n;
	    angle_=180/pi*acos(abs( drij(1)/drij_n  ));
	    edge_angles(count,5)=angle_;
	end
    end
end

fp=fopen('cluster_angles.txt','w');
for i=1:Nedges
	fprintf(fp,'%d %d %d %f %f\n',edge_angles(i,1),edge_angles(i,2),...
	           edge_angles(i,3),edge_angles(i,4),edge_angles(i,5));
end
fclose(fp);

fp=fopen('cluster_connectivity.net','w');
fprintf(fp,'*vertices %d\n',Ncluster);
for i=1:Ncluster
    fprintf(fp,'%d \"%d\" \n',i,cluster_sizes(i));
end
fprintf(fp,'*edges\n');
for i=1:Ncluster
    for j=1:Ncluster
        if (cluster_connectivity(i,j)>0)
	    fprintf(fp,'%d %d\n',i,j);
        end   
    end
end
fclose(fp);
                                                                                                   get_coordinate.m                                                                                    0000664 0001750 0001750 00000005515 12073613501 013336  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% connectivity analysis file            %%%
%%% for PDB file with orthogonal cell     %%%
%%% 1st version, 2012-07-31 Seunghwa Ryu  %%%
%%% 2nd version, 2012-08-13 Seunghwa Ryu  %%%
%%%  works only if no. of peptides < 2000 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_recursion_depth(30000);

rcut = 0.85;
NC_MAX=300;
load('Nbead.txt');
filename ='coord.pdb';

num_a_bead = 0;

% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);
H22 = data(2);
H33 = data(3);
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

fid__ = fopen('boundary_condition.txt','w');
fprintf(fid__,'%f %f %f',H11,H22,H33);
fclose(fid__);

% counting how many a beads in the file
disp('counting number of hydrophobic beads');
for i = 1:Nbead
line = fgets(fid);
data = sscanf(line,'%*s %*s %d %*s %*d %d %f %f %f %*f %*f %*f')';
    if (data(1)==1 || data(1)==4)
	num_a_bead = num_a_bead + 1;
    end
end
fclose(fid);

% declare & read Ndata,
% 1st col : atom index
% 2nd col : atom species, 
% 3rd col : resid
% 4th, 5th, 6th col : x,y,z coordinates
disp('saving the data in memory');
Ndata = zeros(num_a_bead,6);
fid = fopen(filename,'r');line = fgets(fid);
k=0;
for i = 1: Nbead
    line = fgets(fid);
    data = sscanf(line,'%*s %*s %d %*s %*d %d %f %f %f %*f %*f %*f')';
    if (data(1) ==1 || data(1) == 4)
 	k=k+1;
	SR=h\data(3:5)';data(3:5)=SR';
	Ndata(k,2:6) = data(1:end);
	Ndata(k,1)=i;
    end
end
fclose(fid);


% saving data files
disp('saving Ndata nn files');

fp=fopen('Ndata.txt','w');
for i=1:num_a_bead
    fprintf(fp,'%d %d %d %f %f %f\n',Ndata(i,1),Ndata(i,2),Ndata(i,3),Ndata(i,4),Ndata(i,5),Ndata(i,6));
end
fclose(fp);

fp=fopen('coordinate.txt','w');
fprintf(fp,'%f %f %f %d\n',H11,H22,H33,num_a_bead);
key=0;
for i=1:num_a_bead
    if (key==0)
	fprintf(fp,'%f %f %f %d %d\n',Ndata(i,4),Ndata(i,5),Ndata(i,6),...
	    Ndata(i,1),Ndata(i,3));
    else
	if (Ndata(i,3)==999)
            fprintf(fp,'%f %f %f %d %d\n',Ndata(i,4),Ndata(i,5),Ndata(i,6),...
                Ndata(i,1),Ndata(i,3));
	else
            fprintf(fp,'%f %f %f %d %d\n',Ndata(i,4),Ndata(i,5),Ndata(i,6),...
                Ndata(i,1),Ndata(i,3)+1000);
	end
    end
    if (key==0 && Ndata(i,3)==999)
        key=1;
    end
end

fp=fopen('coordinates.txt','w');
key=0;
for i=1:num_a_bead
    if (key==0)
        fprintf(fp,'%f %f %f %d %d\n',Ndata(i,4),Ndata(i,5),Ndata(i,6),...
            Ndata(i,1),Ndata(i,3));
    else
        if (Ndata(i,3)==999)
            fprintf(fp,'%f %f %f %d %d\n',Ndata(i,4),Ndata(i,5),Ndata(i,6),...
                Ndata(i,1),Ndata(i,3));
        else
            fprintf(fp,'%f %f %f %d %d\n',Ndata(i,4),Ndata(i,5),Ndata(i,6),...
                Ndata(i,1),Ndata(i,3)+1000);
        end
    end
    if (key==0 && Ndata(i,3)==999)
        key=1;
    end
end


fclose(fp);

                                                                                                                                                                                   separate_coordinate.m                                                                               0000664 0001750 0001750 00000017477 12043562653 014405  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% separate pdb file generation          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen('ref.pdb','r');
i=0;
while feof(fid) == 0
    line = fgetl(fid);
   i=i+1;
   %do something with line
end
fclose(fid);
Nbead=i-2;

fid = fopen('Nbead.txt','w');
fprintf(fid,'%d',Nbead);
fclose(fid);

filename ='equil.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:5
    dirname=sprintf('equil_evolve_%d',i);
    mkdir(dirname);
end

% counting how many a beads in the file

for j=0:5

filename=sprintf('equil_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end

fclose(fid);

%%% shear file
filename ='shear1.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:7
    dirname=sprintf('shear1_evolve_%d',i);
    mkdir(dirname);
end

for j=0:7

filename=sprintf('shear1_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end
fclose(fid);

%%% shear file
filename ='shear2.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:7
    dirname=sprintf('shear2_evolve_%d',i);
    mkdir(dirname);
end

for j=0:7

filename=sprintf('shear2_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end
fclose(fid);
%%% shear file
filename ='shear3.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:7
    dirname=sprintf('shear3_evolve_%d',i);
    mkdir(dirname);
end

for j=0:7

filename=sprintf('shear3_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end
fclose(fid);
%%% shear file
filename ='shear4.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:7
    dirname=sprintf('shear4_evolve_%d',i);
    mkdir(dirname);
end

for j=0:7

filename=sprintf('shear4_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end
fclose(fid);
%%% shear file
filename ='shear5.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:7
    dirname=sprintf('shear5_evolve_%d',i);
    mkdir(dirname);
end

for j=0:7

filename=sprintf('shear5_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end
fclose(fid);


%%% stretch file
filename ='stretch1.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:5
    dirname=sprintf('stretch1_evolve_%d',i);
    mkdir(dirname);
end

% counting how many a beads in the file
for j=0:5

filename=sprintf('stretch1_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');

H11=H11o*(1+0.2*j);
H22=sqrt(H11o*H22o*H33o/H11);
H33=sqrt(H11o*H22o*H33o/H11);

fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end
fclose(fp);
end
fclose(fid);


filename ='stretch2.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:5
    dirname=sprintf('stretch2_evolve_%d',i);
    mkdir(dirname);
end

% counting how many a beads in the file
for j=0:5

filename=sprintf('stretch2_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');

H11=H11o*(1+0.2*j);
H22=sqrt(H11o*H22o*H33o/H11);
H33=sqrt(H11o*H22o*H33o/H11);

fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end
fclose(fp);
end
fclose(fid);


filename ='stretch3.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:5
    dirname=sprintf('stretch3_evolve_%d',i);
    mkdir(dirname);
end

% counting how many a beads in the file
for j=0:5

filename=sprintf('stretch3_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');

H11=H11o*(1+0.2*j);
H22=sqrt(H11o*H22o*H33o/H11);
H33=sqrt(H11o*H22o*H33o/H11);

fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end
fclose(fp);
end
fclose(fid);


filename ='stretch4.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:5
    dirname=sprintf('stretch4_evolve_%d',i);
    mkdir(dirname);
end

% counting how many a beads in the file
for j=0:5

filename=sprintf('stretch4_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');

H11=H11o*(1+0.2*j);
H22=sqrt(H11o*H22o*H33o/H11);
H33=sqrt(H11o*H22o*H33o/H11);

fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end
fclose(fp);
end
fclose(fid);


filename ='stretch5.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:5
    dirname=sprintf('stretch5_evolve_%d',i);
    mkdir(dirname);
end

% counting how many a beads in the file
for j=0:5

filename=sprintf('stretch5_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');

H11=H11o*(1+0.2*j);
H22=sqrt(H11o*H22o*H33o/H11);
H33=sqrt(H11o*H22o*H33o/H11);

fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end
fclose(fp);
end
fclose(fid);


                                                                                                                                                                                                 separate_coordinate_single.m                                                                        0000664 0001750 0001750 00000005023 12073616644 015731  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% separate pdb file generation          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nequil_shot=3;
Nshear_shot=3;
Nstretch_shot=5;

fid = fopen('ref.pdb','r');
i=0;
while feof(fid) == 0
    line = fgetl(fid);
   i=i+1;
   %do something with line
end
fclose(fid);
Nbead=i-2;

fid = fopen('Nbead.txt','w');
fprintf(fid,'%d',Nbead);
fclose(fid);

filename ='equil.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:Nequil_shot
    dirname=sprintf('equil_evolve_%d',i);
    mkdir(dirname);
end

% write coord.pdb for each directory
for j=0:Nequil_shot

filename=sprintf('equil_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end

fclose(fid);

%%% shear file
filename ='shear1.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:Nshear_shot
    dirname=sprintf('shear1_evolve_%d',i);
    mkdir(dirname);
end

% write coord.pdb for each directory
for j=0:Nshear_shot

filename=sprintf('shear1_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');
H11=H11o;
H22=H22o;
H33=H33o;
fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end

fclose(fp);
end
fclose(fid);


%%% stretch file
filename ='stretch1.pdb';
% read file
fid = fopen(filename,'r');
line = fgets(fid);
data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
H11 = data(1);H11o=H11;
H22 = data(2);H22o=H22;
H33 = data(3);H33o=H33;
h = [ H11 0 0 ; 0 H22 0 ; 0 0 H33 ];

first_line=line;

for i=0:Nstretch_shot
    dirname=sprintf('stretch1_evolve_%d',i);
    mkdir(dirname);
end

% write coord.pdb for each directory
for j=0:Nstretch_shot

filename=sprintf('stretch1_evolve_%d/coord.pdb',j);

fp= fopen(filename,'w');

H11=H11o*(1+0.2*j);
H22=sqrt(H11o*H22o*H33o/H11);
H33=sqrt(H11o*H22o*H33o/H11);

fprintf(fp,'CRYST1 %9.3f %9.3f %9.3f 90.00 90.00 90.00 P 1     1\n',H11,H22,H33);

for i = 1:Nbead+1
    line = fgets(fid);
    fputs(fp,line);
end
fclose(fp);
end
fclose(fid);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             network_noprint                                                                                     0000775 0001750 0001750 00000031446 12043065401 013361  0                                                                                                    ustar   shryu                           shryu                                                                                                                                                                                                                  ELF          >    �@     @       �!          @ 8 	 @         @       @ @     @ @     �      �                   8      8@     8@                                          @       @     $      $                          `     `     H      �                    @      @`     @`     �      �                   T      T@     T@     D       D              P�td         @     @     4       4              Q�td                                                  R�td         `     `     �      �             /lib64/ld-linux-x86-64.so.2          GNU                       GNU �`�t^�Sͽz��t�@Ԗ5                                                     .                      O                      e                      ]                      I                      l                      d                                             V                      C                                             4                       __gmon_start__ _Jv_RegisterClasses libm.so.6 round sqrt libc.so.6 fopen fputc fclose malloc fscanf fprintf __libc_start_main GLIBC_2.2.5                       9          ui	   ~         $          ui	   ~       �`                     `                    `                    `                    `                     `                   ( `                   0 `                   8 `        	           @ `        
           H `                   H����   �  �=  H���        �5�  �%�  @ �%�  h    ������%�  h   ������%�  h   ������%�  h   �����%�  h   �����%�  h   �����%�  h   �����%�  h   �p����%�  h   �`����%�  h	   �P���1�I��^H��H���PTI�� @ H���@ H�ǚ@ �������H��H��  H��t��H��Ð������������UH��SH���=P   uK�0` H�J  H��(` H��H��H9�s$fD  H��H�%  ��(` H�  H9�r��  H��[]�fff.�     H�=�   UH��t�    H��t]�8` ��]Ð�UH��SH��(H�}�H��  H�U�H��H�H� H����   H��  H�U�H��H�H�    H��  H��H��  H�E�    �{H�}  H��  H�E�H�H��H��H��H)�H��HE�H��H�H� H��H�H� H��u4H�K  H�E�H�H��H��H��H)�H��HE�H��H�H� H���*���H�E�H�
  H�U�H��H�H� H;E������a���H��([]�UH��SH���   ��@ ��@ �����H��p���H�� ���H��p�����@ H�Ǹ    �d���H��(���H��p�����@ H�Ǹ    �D���H��0���H��p�����@ H�Ǹ    �$���H��8���H��p�����@ H�Ǹ    ����H��8�����0�����(����� ���H�ƿ�@ �   �����H��8���H��H��� ���H��x���H��8���H��H�������H�E�H��8���H��H�������H�E�H��8���H��H������H�E�H��8���H��H������H�E�H��8���H��H������H��  H��8���H��H��H��H��H)�H��H���`���H�a  H��8���H��H��H��H��H)�H��H���6���H�?  H��8���H��H������H�  H��8���H��H������H�E�H��8���H��H�������H�E�H��8���Hi�`	  H�������H�E�H��8���Hi�`	  H������H�E��`	  ����H�E�Hǅ@���    ��   H��@���H��H��H�x���H��p�����@ H�Ǹ    �'���H��@���H��H��HU�H��p�����@ H�Ǹ    �����H��@���H��H��HU�H��p�����@ H�Ǹ    �����H��@���H��H��HU�H��p�����@ H�Ǹ    ����H��@���H��H��HU�H��p�����@ H�Ǹ    �{���H��@���H��8���H9�@�����������H��p���H���+�����@ ��@ ����H��p���Hǅ@���    �Z  HǅX���    HǅH���    �  H��@���H;�H�����  H��@���H��H�x���� H��H���H��H�x�����\��E�H��@���H��HE�� H��H���H��HE���\��E�H��@���H��HE�� H��H���H��HE���\��E��E��)����M�f(��\�f(��E��E��	����M�f(��\�f(��E��E�������M�f(��\�f(��E��� ����YE��� ����Y�f(��YM���(����YE���(����Y��YE��X���0����YE���0����Y��YE��X��E��QE�f.�zf.�t
�E�������E��@  f.E�����u>H��@���H��HE�H� ��H��H���H��HE�H� ��)��ȉ���1�)Ѓ���   H��H���H��HE�H�H��H���H�pH��p���H��H��@ H�Ǹ    ����H�+  H��@���H�H��H��H��H)�H��H�X���H��H�H��H���H�H��  H��@���H�H��H��H��H)�H��H�X���H��H�H��H���H��HE�H� H�H��X���H��H���H��8���H9�H��������:���HǅP���    �+H��p����    �    ��@ H�Ǹ    �D���H��P����   H+�X���H;�P�������u�H�+  H��@���H��H�H��X���H�H��@���H��8���H9�@������������H��p���H��������@ ��@ �����H��p���Hǅ@���    �9H��  H��@���H��H�H�H��p�����@ H�Ǹ    ����H��@���H��8���H9�@�������u�H��p���H������Hǅh���    Hǅ`���    Hǅ@���    �  H�"      HǅH���    �HH��H���H��HE�H�  H��H���H��H�H�H�H��H���H��HE�H�     H��H���H��8���H9�H�������u�H��@���H�������H��  H���i  HǅH���    �MH��H���H��HE�H�~  H��H���H��H�H�
H��H���H��HU�H�H��H)�H��H�H��H���H��8���H9�H�������u�H��`���H��HE�H�  H�HǅX���    HǅH���    �   H��H���H��HE�H� H��ugH��8���H��`���H�X���H��HE�H��H���H��H�H��8���H��`���H�X���H��HE�H��H���H��HU�H�H�H��X���H��H���H��8���H9�H��������`���H��`���H�G  H;�h���~H�7  H��h���H��@���H��8���H9�@��������������@ ��@ � ���H��p���Hǅ@���    �3H��@���H��HE�H�H��p�����@ H�Ǹ    �����H��@���H��@���H;�`�������u�H��p���H���E�����@ ��@ ����H��p���Hǅ@���    �   HǅH���    �BH��8���H��@���H�H���H��HE�H�H��p�����@ H�Ǹ    �(���H��H���H��H���H;�h�������u�H��p���H�ƿ
   �����H��@���H��@���H;�`��������g���H��p���H���w�����@ ��@ �����H��p���Hǅ@���    �   HǅH���    �BH��8���H��@���H�H���H��HE�H�H��p�����@ H�Ǹ    �Z���H��H���H��H���H;�h�������u�H��p���H�ƿ
   �	���H��@���H��@���H;�`��������g���H��p���H�������    H���   []Ð���������H�l$�L�d$�H�-s	  L�%l	  L�l$�L�t$�L�|$�H�\$�H��8L)�A��I��H��I������H��t1�@ L��L��D��A��H��H9�u�H�\$H�l$L�d$L�l$ L�t$(L�|$0H��8��    �Ð�������������UH��SH��H��  H���t�` D  H����H�H���u�H��[]Ð�H������H���        r coordinate.txt %lf %ld %lf %lf %lf %ld 
 w nn.txt %ld %ld
 %d %d
 nnmax.txt %ld
 cluster_sizes.txt clusters.txt %ld  clusters_resid.txt       333333�?;4      ����P   |���x   �����   x����   ����              zR x�  $      �����    FJw� ?;*3$"    $   D   ����   A�Cg��       $   l   �����   A�CR��      $   �   �����    Q��_@����X      �   ���                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ��������        ��������                       $              9              �@            h@     ���o    �@            �@            �@     
       �                                           �`            �                            �@            �@                   	              ���o    �@     ���o           ���o    z@                                                                                                             @`                     @     &@     6@     F@     V@     f@     v@     �@     �@     �@                     GCC: (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.gnu.build-id .gnu.hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .ctors .dtors .jcr .dynamic .got .got.plt .data .bss .comment                                                                              8@     8                                    #             T@     T                                     1             t@     t      $                              D   ���o       �@     �                                   N             �@     �      8                          V             �@     �      �                              ^   ���o       z@     z                                  k   ���o       �@     �      @                            z             �@     �                                  �             �@     �      �                           �             �@     �                                    �              @            �                             �             �@     �      �                             �             h@     h                                    �             x@     x      �                              �             @           4                              �             P@     P      �                              �             `                                         �             (`     (                                    �             8`     8                                    �             @`     @      �                           �             �`     �                                   �             �`     �      h                             �             P `     P                                     �             ` `     `       8                              �      0               `       *                                                   �       �                                                    )      P         .                 	                      X0      �                                                           8@                   T@                   t@                   �@                   �@                   �@                   z@                   �@                  	 �@                  
 �@                   �@                    @                   �@                   h@                   x@                   @                   P@                   `                   (`                   8`                   @`                   �`                   �`                   P `                   ` `                                        �@                 ��                     `             *     (`             8     8`             E      @             [     ` `            j     h `            x     p@                 ��                �      `             �      @             �     8`             �     0@             �    ��                �      `             �     @`             �      `             �     �`                  @                                  3     P `             >   ��` `             E                     Y    h@             _    �@           g                     {    p `            �                     �   0`             �                     �                     �    P `             �    � `            �                                              X `             #    x@            2    �@     �       B                     V   ��� `             [    � `            e    �@             l   ��` `             x    � `            �    �@     �      �                     �                      �    x `            �                     �    �@              call_gmon_start crtstuff.c __CTOR_LIST__ __DTOR_LIST__ __JCR_LIST__ __do_global_dtors_aux completed.6531 dtor_idx.6533 frame_dummy __CTOR_END__ __FRAME_END__ __JCR_END__ __do_global_ctors_aux network.cpp __init_array_end _DYNAMIC __init_array_start _GLOBAL_OFFSET_TABLE_ __libc_csu_fini round@@GLIBC_2.2.5 data_start _edata fclose@@GLIBC_2.2.5 _fini _Z3DFSl printf@@GLIBC_2.2.5 nbead_in_cluster fscanf@@GLIBC_2.2.5 __DTOR_END__ fputc@@GLIBC_2.2.5 __libc_start_main@@GLIBC_2.2.5 __data_start nnmax fprintf@@GLIBC_2.2.5 __gmon_start__ __dso_handle _IO_stdin_used __libc_csu_init malloc@@GLIBC_2.2.5 _end nn_atomid _start __bss_start nn_resid main fopen@@GLIBC_2.2.5 _Jv_RegisterClasses n_marker sqrt@@GLIBC_2.2.5 _init                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           