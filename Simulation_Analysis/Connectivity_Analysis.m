close all;
% get equil data
timestep=0.03;dumpfreq=15000; % Change these to reflect parameters of equilibration run.
Nequil=7; % Change to reflect the number of frames in equil_{randomseed}.dcd.
Nshear=14; % Change to reflect the number of frames in shear_{randomseed}.dcd.
Nstretch=5; % Change to reflect the number of frames in stretch_{randomseed}.dcd.
Equil_data=zeros(Nequil,4);
for index=1:Nequil
	dirname=sprintf('equil_evolve_%d',index);
	checkname1=sprintf('%s/connectivity_analysis_numeric.txt',dirname);
	checkname2=sprintf('%s/cluster_sizes.txt',dirname);
	temp=isfile(checkname1) && isfile(checkname2);
    if (temp>0)
		load(sprintf('%s/connectivity_analysis_numeric.txt',dirname));
		load(sprintf('%s/cluster_sizes.txt',dirname));
		Equil_data(index,1)=connectivity_analysis_numeric(1);
		Equil_data(index,2)=connectivity_analysis_numeric(5);
		Equil_data(index,3)=std(cluster_sizes);
		Equil_data(index,4)=connectivity_analysis_numeric(7);
	else
		Equil_data(index+1, 1:4)=[0;0;0;0];
	end

end

plot_indi = 1; 
equiltime=(1:Nequil)*timestep*dumpfreq;

if (plot_indi>0)
	figure1=figure(1)
	plot(equiltime,Equil_data(:,1),'ro-');
	xlabel('\tau');
	ylabel('Number of Beta Crystals');
	xlim([0 max(equiltime)*1.1]);
	ylim([0 max(Equil_data(:,1))*1.1]);
	print(figure1,'-dpng','num_Beta_equil.png');
	close(figure1);

	figure2=figure(2)
	errorbar(equiltime,Equil_data(:,2),Equil_data(:,3),'bo-');
	xlabel('\tau');
	ylabel('Average Crystal Size');
	xlim([0 max(equiltime)*1.1]);
	print(figure2,'-dpng','ave_size_Crystal_equil.png');
	close(figure2);

	figure3=figure(3)
	plot(equiltime,Equil_data(:,4),'bo-');
	xlabel('\tau');
	ylabel('Number of Connections');
	xlim([0 max(equiltime)*1.1]);
	%ylim([0 max(Equil_data(:,4))*1.1]);
	ylim([0 inf]); %<--
	print(figure3,'-dpng','num_Conn_equil.png');
	close(figure3);

	end


fp=fopen('equil.txt','w');
for i=1:length(equiltime)
    fprintf(fp,'%d %f\n',equiltime(i),Equil_data(i,1));
end
for i=1:length(equiltime)
    fprintf(fp,'%d %f\n',equiltime(i),Equil_data(i,2));
end
for i=1:length(equiltime)
    fprintf(fp,'%d %f\n',equiltime(i),Equil_data(i,4));
end
fclose(fp);

%--------------------------------------
% get shear data
timestep=0.03;dumpfreq=15000; % Change these to reflect parameters of equilibration run.
Shear_data=zeros(Nshear+1,4);

for index=0:Nshear
	dirname=sprintf('shear_evolve_%d',index);
	checkname1=sprintf('%s/connectivity_analysis_numeric.txt',dirname);
	checkname2=sprintf('%s/cluster_sizes.txt',dirname);
	temp=isfile(checkname1) && isfile(checkname2);
    if (temp>0)
		load(sprintf('%s/connectivity_analysis_numeric.txt',dirname));
		load(sprintf('%s/cluster_sizes.txt',dirname));
		Shear_data(index+1,1)=connectivity_analysis_numeric(1);
		Shear_data(index+1,2)=connectivity_analysis_numeric(5);
		Shear_data(index+1,3)=std(cluster_sizes);
		Shear_data(index+1,4)=connectivity_analysis_numeric(7);
	else
		Shear_data(index+1, 1:4)=[0;0;0;0];
	end
end

sheartime=(0:Nshear)*timestep*dumpfreq;

if (plot_indi>0)

	figure4=figure(4)
	plot(sheartime,Shear_data(:,1),'ro-');
	xlabel('\tau');
	ylabel('Number of Beta Crystals');
	xlim([0 max(sheartime)*1.1]);
	%ylim([0 max(Shear_data(:,1))*1.1]);
	ylim([0 inf]); %<--
	print(figure4,'-dpng','num_Beta_shear.png');
	close(figure4);

	figure5=figure(5)
	plot(sheartime,Shear_data(:,2),'ro-');
	%errorbar(sheartime,Shear_data(:,2,1),Shear_data(:,3,1),'bo-');
	xlabel('\tau');
	ylabel('Average Crystal Size');
	xlim([0 max(sheartime)*1.1]);
	print(figure5,'-dpng','ave_size_Crystal_shear.png');
	close(figure5);

	figure6=figure(6)
	plot(sheartime,Shear_data(:,4),'ro-');
	xlabel('\tau');
	ylabel('Number of Connections');
	xlim([0 max(sheartime)*1.1]);
	%ylim([0 max(Shear_data(:,4))*1.1]);
	ylim([0 inf]); %<--
	print(figure6,'-dpng','num_Conn_shear.png');
	close(figure6);

	end

fp=fopen('shear.txt','w');
for i=1:length(sheartime)
    fprintf(fp,'%d %f\n',sheartime(i),Shear_data(i,1));
end
for i=1:length(sheartime)
    fprintf(fp,'%d %f\n',sheartime(i),Shear_data(i,2));
end
for i=1:length(sheartime)
    fprintf(fp,'%d %f\n',sheartime(i),Shear_data(i,4));
end
fclose(fp);

%-------------------------------------------
% get stretch data
timestep=0.03;dumpfreq=10000; % Change these to reflect parameters of equilibration run.
Stretch_data=zeros(Nstretch+1,4);
for index=0:Nstretch
	dirname=sprintf('stretch_evolve_%d',index);
	checkname1=sprintf('%s/connectivity_analysis_numeric.txt',dirname);
	checkname2=sprintf('%s/cluster_sizes.txt',dirname);
	temp=isfile(checkname1) && isfile(checkname2);
    if (temp>0)
		load(sprintf('%s/connectivity_analysis_numeric.txt',dirname));
		load(sprintf('%s/cluster_sizes.txt',dirname));
		Stretch_data(index+1,1)=connectivity_analysis_numeric(1);
		Stretch_data(index+1,2)=connectivity_analysis_numeric(5);
		Stretch_data(index+1,3)=std(cluster_sizes);
		Stretch_data(index+1,4)=connectivity_analysis_numeric(7);
	else
		Stretch_data(index+1, 1:4)=[0;0;0;0];
	end

end

stretchtime=(0:Nstretch)*timestep*dumpfreq;

fp=fopen('stretch.txt','w');
for i=1:length(stretchtime)
    fprintf(fp,'%d %f\n',stretchtime(i),Stretch_data(i,1));
end
for i=1:length(stretchtime)
    fprintf(fp,'%d %f\n',stretchtime(i),Stretch_data(i,2));
end
for i=1:length(stretchtime)
    fprintf(fp,'%d %f\n',stretchtime(i),Stretch_data(i,4));
end
fclose(fp);

if (plot_indi>0)

	figure7=figure(7)
	plot(stretchtime,Stretch_data(:,1),'ro-');
	xlabel('\tau');
	ylabel('Number of Beta Crystals');
	xlim([0 max(stretchtime)*1.1]);
	%ylim([0 max(Stretch_data(:,1))*1.1]);
	ylim([0 inf]); %<--
	print(figure7,'-dpng','num_Beta_stretch.png');
	close(figure7);

	figure8=figure(8)
	plot(stretchtime,Stretch_data(:,2),'ro-');
	%errorbar(stretchtime,Stretch1_data(:,2),Stretch1_data(:,3),'bo-');
	xlabel('\tau');
	ylabel('Average Crystal Size');
	xlim([0 max(stretchtime)*1.1]);
	print(figure8,'-dpng','ave_size_Crystal_stretch.png');
	close(figure8);

	figure9=figure(9)
	plot(stretchtime,Stretch_data(:,4),'ro-');
	xlabel('\tau');
	ylabel('Number of Connections');
	xlim([0 max(stretchtime)*1.1]);
	%ylim([0 max(Stretch_data(:,4))*1.1]);
	ylim([0 inf]); %<--
	print(figure9,'-dpng','num_Conn_stretch.png');
	close(figure9);

	end

%disp(sprintf('# of beta crystal cross-links ==>\n E %d Shear %d Stretch %d',...
    %(Shear_data(1,1)),((Stretch_data(1,1))),((Stretch_data(5,1)))));
%disp(sprintf('average beta crystal sizes ==>\n E %d Shear %d Stretch %d',...
    %(Shear_data(1,2)),((Stretch_data(1,2))),((Stretch_data(5,2)))));
%disp(sprintf('number of likes between crystals ==>\n E %d Shear %d Stretch %d',...
    %(Shear_data(1,4)),((Stretch_data(1,4))),((Stretch_data(5,4)))));