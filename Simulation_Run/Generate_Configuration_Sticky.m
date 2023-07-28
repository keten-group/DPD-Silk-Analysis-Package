%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spider Silk Experiment Configuration Generation Code %%%
%%% 1st version, July 17th, 2012 Seunghwa Ryu            %%%
%%% 2nd version, July 20th, 2012 Seunghwa Ryu            %%%
%%% 3rd version, Aug 8th, 2012 Seunghwa Ryu              %%%
%%%              silkworm silk first in the data file    %%%
%%% 4th version, Oct 22th, 2012 Seunghwa Ryu             %%%
%%%              final protocal level                    %%%
%%% 5th version, Jan 08th, 2013 Seunghwa Ryu             %%%
%%%              publication level, with 2 bond styles   %%%
%%%                                                      %%%
%%%          generate  lammps data input file, psf file  %%%
%%%       still this code does not have angle & dihedral %%%
%%% 6th version, July 24th, 2023 Timothy Russell         %%%
%%%         added capability to generate 'sticky'        %%%
%%%         terminal ends in the form of bead type 'c'.  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ALL ADDITIONS TO ORIGINAL SCRIPT DICTATED WITH "XXX" %%%
%%% including additional comments for readability        %%%
%%% questions to come back to a represented with "???"   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% input variable section               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total number of bead type, change when different indexes are used for
% spider h, silkworm a, silkworm b beads

% default choice in the paper
spider_a_index=1;   spider_b_index=2;   w_index = 3; 

%%% XXX new index for the sticky terminal end, represented by the c bead
%%% XXX for now both of the sticky ends are represented by the same bead 
%%% XXX type, there is no distinction between N
spider_c_index=6;


%%% XXX histadine tag represented as a b bead
nbeadtype=6; spider_h_index = 2; silkworm_a_index=4; silkworm_b_index=5; 

% other choices
% nbeadtype=6; spider_h_index = 6; silkworm_a_index=4; silkworm_b_index=5;
% nbeadtype=3; spider_h_index = 2; silkworm_a_index=1; silkworm_b_index=2;


%%% ??? unsure if the number of bond types holds or if new bond types  
%%% ??? are required to model the sticky c beads
nbondtype=2; % total number of bond type, 2 for 'a'-'a', 1 for others
rho=3; % density of DPD bead in the simulation         
d=1.025; % initial polymer bead-to-bead distance
boxsize_x=60.0;
boxsize_y=40.0;
boxsize_z=40.0;

% spider silk parameter
num_repeat_spider=2; % number of repeat motif
num_histidine_spider=16;  % number of histidine bead in H motif, default=16
num_hydrophobic_spider=5; % number of hydrophobic bead in A (hydrophobic) motif
num_hydrophilic_spider=7; % number of hydrophilic bead in B (hydrophilic) motif

%%% XXX The number of sticky end beads in each end of the termninal motif
num_sticky_spider=0;


repeat_motif_spider='A1B1'; % repeat motif description ( A1B3, A1B1, A2B2 for example )

% silkworm silk parameter
num_repeat_silkworm=1;  % number of repeat motif
num_hydrophobic_silkworm=10; % number of hydrophobic bead in A (hydrophobic) motif
num_hydrophilic_silkworm=10;  % number of hydrophilic bead in B (hydrophilic) motif
repeat_motif_silkworm='A1B1'; % default, no change, very simple assumption that beta crystal ~ 50%

% fraction of silk polymer
fraction_of_silk_peptide=0.2; % volume percent of peptides in the solution
relative_fraction_of_spider_silk=1.0; % relative volume percent of silk peptides in the solution
relative_fraction_of_silkworm_silk=1-relative_fraction_of_spider_silk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end of input variable                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% beyond this point, processing starts %%%

% obtaining the full sequence for spider
num_bead_repeat_motif_spider=(...
    (str2num(repeat_motif_spider(2)))*(...
    (repeat_motif_spider(1)=='A')*num_hydrophobic_spider+...
    (repeat_motif_spider(1)=='B')*num_hydrophilic_spider)+...
    (str2num(repeat_motif_spider(4)))*(...
    (repeat_motif_spider(3)=='A')*num_hydrophobic_spider+...
    (repeat_motif_spider(3)=='B')*num_hydrophilic_spider)...
    );

%%% XXX the number of terminal beads of each polymer is added to the 
%%% XXX length of the polymer
num_bead_eachpolymer_spider=...
    num_histidine_spider+...
    num_sticky_spider*2+... %%% XXX 
    num_bead_repeat_motif_spider*...
    num_repeat_spider;

spider_bead_type=spider_b_index*ones(1,num_bead_eachpolymer_spider);

%%% XXX The type index which is used to begin of the location of repeating
%%% XXX motif of the spider silk now needs to include the length of the
%%% XXX sticky terminal end because there are more beads before the repeat
%%% XXX sequence of the peptide
%%% XXX the sticky terminal end will be dictated with a 'C'
typeindex=num_sticky_spider+num_histidine_spider;
if (num_histidine_spider>0 & num_sticky_spider==0)
    full_spider_sequence='H';
    spider_bead_type(1:num_histidine_spider)= spider_h_index*ones(1,num_histidine_spider);
%%% XXX if the both the number of histidine beads and sticky beads are
%%% XXX greater than 0 add the sticky end to the beginning of the sequence
%%% XXXthen the histidine tag
elseif (num_histidine_spider>0 & num_sticky_spider>0)
    full_spider_sequence='HC';
%%% XXX the initial group of spider beads will be C beads
    spider_bead_type(1:num_histidine_spider)=spider_h_index*ones(1,num_histidine_spider);
%%% XXX after the initial C beads then the histidine beads come next
    spider_bead_type(num_histidine_spider+1:num_histidine_spider+num_sticky_spider)=spider_c_index*ones(1,num_sticky_spider);
%%% XXX the final group will be the other set of C beads
    spider_bead_type(num_bead_eachpolymer_spider-num_sticky_spider+1:end)=spider_c_index*ones(1,num_sticky_spider);
elseif (num_histidine_spider==0 & num_sticky_spider>0)
%%% XXX no histidine tag so the beginning of the peptide is only C beads
    full_spider_sequence='C';
    spider_bead_type(1:num_sticky_spider)=spider_c_index*ones(1,num_sticky_spider);
%%% XXX the end of the peptide is still C beads
    spider_bead_type(num_bead_eachpolymer_spider-num_sticky_spider+1:end)=spider_c_index*ones(1,num_sticky_spider);
else
    full_spider_sequence='';
end

for i=1:num_repeat_spider
    full_spider_sequence=strcat(full_spider_sequence,repeat_motif_spider);
    if (repeat_motif_spider(1)=='A')
        spider_a_bead_num=str2num(repeat_motif_spider(2))*num_hydrophobic_spider;
        spider_b_bead_num=num_bead_repeat_motif_spider-spider_a_bead_num;
        spider_bead_type(typeindex+1:typeindex+spider_a_bead_num)=...
                spider_a_index*ones(1,spider_a_bead_num);
    else
        spider_a_bead_num=str2num(repeat_motif_spider(4))*num_hydrophobic_spider;
        spider_b_bead_num=num_bead_repeat_motif_spider-spider_a_bead_num;
        spider_bead_type(typeindex+1+spider_b_bead_num:...
                typeindex+spider_b_bead_num+spider_a_bead_num)=...
                spider_a_index*ones(1,spider_a_bead_num);
    end
    typeindex=typeindex+spider_a_bead_num+spider_b_bead_num;
end

%%% XXX if the length of the sticky beads is greater than 0 add the other 
%%% XXX sticky end to full spider sequence
if (num_sticky_spider>0)
    full_spider_sequence=strcat(full_spider_sequence,'C');
end

% obtaining full sequence for silkworm
num_bead_repeat_motif_silkworm=(...
    (str2num(repeat_motif_silkworm(2)))*(...
    (repeat_motif_silkworm(1)=='A')*num_hydrophobic_silkworm+...
    (repeat_motif_silkworm(1)=='B')*num_hydrophilic_silkworm)+...
    (str2num(repeat_motif_silkworm(4)))*(...
    (repeat_motif_silkworm(3)=='A')*num_hydrophobic_silkworm+...
    (repeat_motif_silkworm(3)=='B')*num_hydrophilic_silkworm)...
    );
num_bead_eachpolymer_silkworm=...
    num_bead_repeat_motif_silkworm*...
    num_repeat_silkworm;

silkworm_bead_type=silkworm_b_index*ones(1,num_bead_eachpolymer_silkworm);typeindex=0;
full_silkworm_sequence='';
for i=1:num_repeat_silkworm
    full_silkworm_sequence=strcat(full_silkworm_sequence,repeat_motif_silkworm);
    if (repeat_motif_silkworm(1)=='A')
        silkworm_a_bead_num=str2num(repeat_motif_silkworm(2))*num_hydrophobic_silkworm;
        silkworm_b_bead_num=num_bead_repeat_motif_silkworm-silkworm_a_bead_num;
        silkworm_bead_type(typeindex+1:typeindex+silkworm_a_bead_num)=...
                silkworm_a_index*ones(1,silkworm_a_bead_num);
    else
        silkworm_a_bead_num=str2num(repeat_motif_silkworm(4))*num_hydrophobic_silkworm;
        silkworm_b_bead_num=num_bead_repeat_motif_silkworm-silkworm_a_bead_num;
        silkworm_bead_type(typeindex+1+silkworm_b_bead_num:...
                typeindex+silkworm_b_bead_num+silkworm_a_bead_num)=...
                silkworm_a_index*ones(1,silkworm_a_bead_num);
    end
    typeindex=typeindex+num_bead_repeat_motif_silkworm;
end

% determine number of beads in the simulation. 
boxvolume=boxsize_x*boxsize_y*boxsize_z;
num_tot_beads=floor(rho*boxvolume);
num_tot_silk_beads_temp=floor(num_tot_beads*fraction_of_silk_peptide);

num_tot_spider_beads_temp=floor(num_tot_silk_beads_temp*...
					relative_fraction_of_spider_silk);
if (relative_fraction_of_spider_silk==0)
    num_tot_spider_beads_temp=0;
end

num_tot_silkworm_beads_temp=num_tot_silk_beads_temp-num_tot_spider_beads_temp;
if (relative_fraction_of_silkworm_silk==0)
    num_tot_silkworm_beads_temp=0;
end

num_spider_peptide=ceil(num_tot_spider_beads_temp/num_bead_eachpolymer_spider);
num_silkworm_peptide=ceil(num_tot_silkworm_beads_temp/num_bead_eachpolymer_silkworm);

num_tot_spider_beads=num_spider_peptide*num_bead_eachpolymer_spider;
num_tot_silkworm_beads=num_silkworm_peptide*num_bead_eachpolymer_silkworm;
num_tot_silk_beads=num_tot_spider_beads+num_tot_silkworm_beads;


%%% determine the filename for the datafile
datafilename=strcat('spider_',full_spider_sequence,'x',num2str(num_spider_peptide),...
		'_silkworm_',full_silkworm_sequence,'x',num2str(num_silkworm_peptide),...
                '.data');
psffilename=strcat('spider_',full_spider_sequence,'x',num2str(num_spider_peptide),...
                '_silkworm_',full_silkworm_sequence,'x',num2str(num_silkworm_peptide),...
                '.psf');
psffilename_protein_only=strcat('spider_',full_spider_sequence,'x',num2str(num_spider_peptide),...
                '_silkworm_',full_silkworm_sequence,'x',num2str(num_silkworm_peptide),...
                '_protein_only.psf');
psffilename_evap=strcat('spider_',full_spider_sequence,'x',num2str(num_spider_peptide),...
                '_silkworm_',full_silkworm_sequence,'x',num2str(num_silkworm_peptide),...
                '_evap.psf');
psffilename_evap_water_only=strcat('spider_',full_spider_sequence,'x',num2str(num_spider_peptide),...
                '_silkworm_',full_silkworm_sequence,'x',num2str(num_silkworm_peptide),...
                '_evap_water_only.psf');

%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating configuration %
%%%%%%%%%%%%%%%%%%%%%%%%%%

tot_h=num_histidine_spider*num_spider_peptide;
tot_spider_a=num_repeat_spider*spider_a_bead_num*num_spider_peptide;
tot_spider_b=num_repeat_spider*spider_b_bead_num*num_spider_peptide;
tot_silkworm_a=num_repeat_silkworm*silkworm_a_bead_num*num_silkworm_peptide;
tot_silkworm_b=num_repeat_silkworm*silkworm_b_bead_num*num_silkworm_peptide;
%%% XXX add the total number of c beads
tot_spider_c=2*num_sticky_spider* num_spider_peptide;
%%% XXX the number of water beads is subtracted from the number of a,b,c
%%% XXX and silkworm beads
tot_water=num_tot_beads-tot_h-tot_spider_a-tot_spider_b-tot_spider_c-tot_silkworm_a-tot_silkworm_b;

%%% XXX the total number of bonds now includes the C beads
num_tot_bonds=num_spider_peptide*(num_bead_eachpolymer_spider-1)+...
	      num_silkworm_peptide*(num_bead_eachpolymer_silkworm-1);
angles=0;
dihedrals=0;
impropers=0;
mass=ones(1,nbeadtype);

bh_x=0.5*boxsize_x;
bh_y=0.5*boxsize_y;
bh_z=0.5*boxsize_z;

%%% creating the data file
disp(sprintf('creating the data file, %s',datafilename));
% write intro
fp=fopen(datafilename,'w');
fprintf(fp,'%s\n',datafilename);
fprintf(fp,'# silk fraction=%f, relative portion of spider silk=%f\n',...
		num_tot_silk_beads/num_tot_beads,num_tot_spider_beads/num_tot_silk_beads);
fprintf(fp,'%d atoms',num_tot_beads);
fprintf(fp,' # h %d, sp-a %d, sp-b %d,',tot_h,tot_spider_a,tot_spider_b);
%%% XXX add the number of c beads to the data file information
fprintf(fp,' c %d,',tot_spider_c);
fprintf(fp,' w %d, sw-a %d, sw-b %d\n',tot_water,tot_silkworm_a,tot_silkworm_b);
fprintf(fp,'%d bonds\n',num_tot_bonds);
fprintf(fp,'%d angles\n',angles);
fprintf(fp,'%d dihedrals\n',dihedrals);
fprintf(fp,'%d impropers\n\n',impropers);
fprintf(fp,'%d atom types\n',nbeadtype);
fprintf(fp,'%d bond types\n',nbondtype);
fprintf(fp,'-%f %f xlo xhi\n',bh_x,bh_x);
fprintf(fp,'-%f %f ylo yhi\n',bh_y,bh_y);
fprintf(fp,'-%f %f zlo zhi\n',bh_z,bh_z);
fprintf(fp,'0.0 0.0 0.0 xy xz yz\n\n');
fprintf(fp,'Masses\n\n');
for i=1:nbeadtype
    fprintf(fp,'%d %f\n',i,mass(i));
end
fprintf(fp,'\n');
fprintf(fp,'Atoms\n\n');

bead_number=0;
polymer_number=0;

% write beads coordinate

for i=1:num_silkworm_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_silkworm
        bead_number=bead_number+1;
        fprintf(fp,'%d %d %d %f %f %f\n',bead_number,polymer_number,silkworm_bead_type(j),x,y,z);
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

for i=1:num_spider_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_spider
        bead_number=bead_number+1;
	fprintf(fp,'%d %d %d %f %f %f\n',bead_number,polymer_number,spider_bead_type(j),x,y,z);
	old_x=x;old_y=y;old_z=z;
	while (1)
 	    theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
   	    if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
		break;
	    else
		x=old_x;y=old_y;z=old_z;
	    end
	end
    end
end

for i=bead_number:num_tot_beads-1

    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;
    bead_number=bead_number+1;
 %%% XXX added a semi colon to prevent print out 
    fprintf(fp,'%d %d %d %f %f %f\n',bead_number,0,w_index,x,y,z);
end

% write bond information
fprintf(fp,'\nBonds\n\n');

bond_number=0;bead_number=0;

for i=1:num_silkworm_peptide
    for j=1:num_bead_eachpolymer_silkworm-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
	if (silkworm_bead_type(j)==silkworm_bead_type(j+1) && silkworm_bead_type(j)==silkworm_a_index )
	    bond_type=2;
	else
            bond_type=1;
	end
        fprintf(fp,'%d %d %d %d\n',bond_number,bond_type,bead_number,bead_number+1);
    end
    bead_number=bead_number+1;
end

for i=1:num_spider_peptide
    for j=1:num_bead_eachpolymer_spider-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        if (spider_bead_type(j)==spider_bead_type(j+1) && spider_bead_type(j)==spider_a_index )
            bond_type=2;
        else
            bond_type=1;
        end	
        fprintf(fp,'%d %d %d %d\n',bond_number,bond_type,bead_number,bead_number+1);
    end
    bead_number=bead_number+1;
end

fclose(fp);


%%% creating the psf file
disp(sprintf('creating the corresponding psf file, %s',psffilename));

% write intro
fp=fopen(psffilename,'w');
fprintf(fp,'PSF %s\n\n',psffilename);
fprintf(fp,'2 !NTITLE\n');
fprintf(fp,'REMARKS  h %d, sp-a %d, sp-b %d,',tot_h,tot_spider_a,tot_spider_b);
%%% XXX add the number of c beads to the psf file information
fprintf(fp,' c %d,',tot_spider_c);
fprintf(fp,' w %d, sw-a %d, sw-b %d\n',tot_water,tot_silkworm_a,tot_silkworm_b);
fprintf(fp,'REMARKS %d bonds\n',num_tot_bonds);
fprintf(fp,'\n');
fprintf(fp,'%d !NATOM\n',num_tot_beads);

bead_number=0;
polymer_number=0;

% write beads
for i=1:num_silkworm_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_silkworm
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,silkworm_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

for i=1:num_spider_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_spider
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,spider_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

for i=bead_number:num_tot_beads-1

    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;
    bead_number=bead_number+1;
    fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,0,w_index);
end

% write bond information
fprintf(fp,'\n %d !NBOND: bonds\n',num_tot_bonds);

bond_number=0;bead_number=0;linecounter=0;

for i=1:num_silkworm_peptide
    for j=1:num_bead_eachpolymer_silkworm-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

for i=1:num_spider_peptide
    for j=1:num_bead_eachpolymer_spider-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

fclose(fp);



%%% creating the water only psf file
disp(sprintf('creating the corresponding psf file, %s',psffilename_protein_only));

% write intro
fp=fopen(psffilename_protein_only,'w');
fprintf(fp,'PSF %s\n\n',psffilename_protein_only);
fprintf(fp,'2 !NTITLE\n');
fprintf(fp,'REMARKS  h %d, sp-a %d, sp-b %d,',tot_h,tot_spider_a,tot_spider_b);
%%% XXX add the number of c beads to the water only psf file information
fprintf(fp,' c %d,',tot_spider_c);
fprintf(fp,' w %d, sw-a %d, sw-b %d\n',tot_water,tot_silkworm_a,tot_silkworm_b);
fprintf(fp,'REMARKS %d bonds\n',num_tot_bonds);
fprintf(fp,'\n');
fprintf(fp,'%d !NATOM\n',num_tot_beads-tot_water);

bead_number=0;
polymer_number=0;

% write beads
for i=1:num_silkworm_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_silkworm
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,silkworm_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end
for i=1:num_spider_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_spider
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,spider_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

% write bond information
fprintf(fp,'\n %d !NBOND: bonds\n',num_tot_bonds);

bond_number=0;bead_number=0;linecounter=0;

for i=1:num_silkworm_peptide
    for j=1:num_bead_eachpolymer_silkworm-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

for i=1:num_spider_peptide
    for j=1:num_bead_eachpolymer_spider-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

fclose(fp);


%%% creating the evap psf file
disp(sprintf('creating the corresponding psf file, %s',psffilename_evap));

boxsize_x=boxsize_x*0.8;
boxsize_y=boxsize_y*0.8;
boxsize_z=boxsize_z*0.8;

% write intro
fp=fopen(psffilename_evap,'w');
fprintf(fp,'PSF %s\n\n',psffilename_evap);
fprintf(fp,'2 !NTITLE\n');
fprintf(fp,'REMARKS  h %d, sp-a %d, sp-b %d,',tot_h,tot_spider_a,tot_spider_b);
%%% XXX add the number of c beads to the evap psf file information
fprintf(fp,' c %d,',tot_spider_c);
fprintf(fp,' w %d, sw-a %d, sw-b %d\n',tot_water,tot_silkworm_a,tot_silkworm_b);
fprintf(fp,'REMARKS %d bonds\n',num_tot_bonds);
fprintf(fp,'\n');
fprintf(fp,'%d !NATOM\n',num_tot_beads);

bead_number=0;
polymer_number=0;

% write beads
for i=1:num_silkworm_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_silkworm
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,silkworm_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

for i=1:num_spider_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_spider
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,spider_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

for i=bead_number:num_tot_beads-1

    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;
    bead_number=bead_number+1;
    fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,0,w_index);
end

% write bond information
fprintf(fp,'\n %d !NBOND: bonds\n',num_tot_bonds);

bond_number=0;bead_number=0;linecounter=0;

for i=1:num_silkworm_peptide
    for j=1:num_bead_eachpolymer_silkworm-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

for i=1:num_spider_peptide
    for j=1:num_bead_eachpolymer_spider-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end

    end
    bead_number=bead_number+1;
end

fclose(fp);

if 0

%%% creating the evap water only psf file
disp(sprintf('creating the corresponding psf file, %s',psffilename_evap_water_only));

% write intro
fp=fopen(psffilename_evap_water_only,'w');
fprintf(fp,'PSF %s\n\n',psffilename_evap_water_only);
fprintf(fp,'2 !NTITLE\n');
fprintf(fp,'REMARKS  h %d, sp-a %d, sp-b %d,',tot_h,tot_spider_a,tot_spider_b);
fprintf(fp,' w %d, sw-a %d, sw-b %d\n',tot_water,tot_silkworm_a,tot_silkworm_b);
fprintf(fp,'REMARKS %d bonds\n',num_tot_bonds);
fprintf(fp,'\n');
fprintf(fp,'%d !NATOM\n',num_tot_beads);

bead_number=0;
polymer_number=0;

% write beads
for i=1:num_silkworm_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_silkworm
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,silkworm_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

for i=1:num_spider_peptide

    polymer_number=polymer_number+1;
    x=0.5*(2*rand-1)*boxsize_x;
    y=0.5*(2*rand-1)*boxsize_y;
    z=0.5*(2*rand-1)*boxsize_z;

    for j=1:num_bead_eachpolymer_spider
        bead_number=bead_number+1;
        fprintf(fp,'%10d %5d 1 XXX %5d CCC 1.000000 1.0000  0\n',bead_number,polymer_number,spider_bead_type(j));
        old_x=x;old_y=y;old_z=z;
        while (1)
            theta=rand*pi;
            phi=rand*2*pi;
            x=x+d*sin(theta)*cos(phi);
            y=y+d*sin(theta)*sin(phi);
            z=z+d*cos(theta);
            if ( abs(x)<bh_x && abs(y)<bh_y && abs(z)<bh_z )
                break;
            else
                x=old_x;y=old_y;z=old_z;
            end
        end
    end
end

% write bond information
fprintf(fp,'\n %d !NBOND: bonds\n',num_tot_bonds);

bond_number=0;bead_number=0;linecounter=0;

for i=1:num_silkworm_peptide
    for j=1:num_bead_eachpolymer_silkworm-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

for i=1:num_spider_peptide
    for j=1:num_bead_eachpolymer_spider-1
        bond_number=bond_number+1;
        bead_number=bead_number+1;
        fprintf(fp,'%10d %10d',bead_number,bead_number+1);
        linecounter=linecounter+1;
        if (linecounter==4)
            fprintf(fp,'\n');
            linecounter=0;
        end
    end
    bead_number=bead_number+1;
end

fclose(fp);

end
