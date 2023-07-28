load('cluster_connectivity.txt');
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
