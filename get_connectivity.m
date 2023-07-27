%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
link_resid_of_each_cluster=zeros(Ncluster,500);

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
	if (nodeid_for_polymer(i)~=nodeid_for_polymer(i+1))
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
fprintf(fp,'Number of nodes = %d\n',Ncluster);
fprintf(fp,'Number of isolated nodes = %d\n',Ncluster_isol);
fprintf(fp,'Max node size = %d\n',cluster_size_MAX);
fprintf(fp,'Mean node size = %f\n',mean(cluster_sizes));
fprintf(fp,'Median node size = %f\n',median(cluster_sizes));
fprintf(fp,'Max bridges of a single node =%d\n',MAX_links_per_cluster);
fprintf(fp,'Number of bridges = %d\n',total_connections);
fprintf(fp,'\n');
fprintf(fp,'Number of clusters of nodes = %d\n',Ncluster_of_nodes);
fprintf(fp,'MAX number of clusters of nodes =%d\n',max(Ncluster_of_nodes_sizes));
fprintf(fp,'\n');
fprintf(fp,'Number of peptides bridging crystals = %d\n',total_connection_res);
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
