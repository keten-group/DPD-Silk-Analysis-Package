%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


