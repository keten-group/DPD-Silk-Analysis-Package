%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% connectivity analysis file            %%%
%%% for PDB file with orthogonal cell     %%%
%%% 1st version, 2012-07-31 Seunghwa Ryu  %%%
%%% 2nd version, 2012-08-13 Seunghwa Ryu  %%%
%%%  works only if no. of peptides < 2000 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%max_recursion_depth(30000);

rcut = 0.85;
load('../Nbead.txt');
filename ='coord.pdb';

num_a_bead = 0;

% read file
fid = fopen(filename,'r');
line = fgets(fid);
fprintf(line)

data = sscanf(line,'%*s %f %f %f %f %f %f %*c %f %f');
%line
%data
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
line = fgets(fid); % every line is an atom info

data = sscanf(line,'%*s %*s %d %*s %*d %d %f %f %f %*f %*f %*f'); % <-- problem: cannot read in values after first %d
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
fid = fopen(filename,'r');
line = fgets(fid);
k=0;

special=1; %<--

for i = 1: Nbead
    line = fgets(fid);

    if (special==1) % for HA3B 
        data = sscanf(line,'%*s %*s %d %*s %*d %d %f %f %f %*f %*f %*f');
        data = data';

    else % for H(AB)2, HAB3
        linesp = strsplit(line, {" "});
        % use {} to extract from cell array
        if (length(linesp{1,4})==1) 
            data=[linesp{1,3}; linesp{1,5}; linesp{1,6}; linesp{1,7}; linesp{1,8}];
        else
            data=[linesp{1,3}; linesp{1,4}(2:end); linesp{1,5}; linesp{1,6}; linesp{1,7}];
        end
        data=str2num(data)';
    end 

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
                Ndata(i,1),Ndata(i,3)+1000*special); %<-- fixed the ACEII char issue, no need for +1000
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
                Ndata(i,1),Ndata(i,3)+1000*special); % <-- no need for +1000
        end
    end
    if (key==0 && Ndata(i,3)==999)
        key=1;
    end
end


fclose(fp);

