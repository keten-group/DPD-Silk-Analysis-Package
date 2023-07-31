%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% separate pdb file generation          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nequil_shot=7; % Change to match number of frames in equil_{randomseed}.dcd
Nshear_shot=14; % Change to match number of frames in shear_{randomseed}.dcd
Nstretch_shot=5; % Change to match number of frames in stretch_{randomseed}.dcd

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
filename ='shear.pdb';
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
    dirname=sprintf('shear_evolve_%d',i);
    mkdir(dirname);
end

% write coord.pdb for each directory
for j=0:Nshear_shot

filename=sprintf('shear_evolve_%d/coord.pdb',j);

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
filename ='stretch.pdb';
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
    dirname=sprintf('stretch_evolve_%d',i);
    mkdir(dirname);
end

% write coord.pdb for each directory
for j=0:Nstretch_shot

filename=sprintf('stretch_evolve_%d/coord.pdb',j);

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

