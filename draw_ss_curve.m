totalsteps = 100000;
volume = 96000;
volume_after_evaporation = volume*0.8*0.8*0.8;
maxstrain = 1;

equil=0
shear=0
evaporation=0

%----------------all_stress_after_equil
if (equil==1)
    fid=fopen('all_stress_after_equil.txt');
    nLines=0;
    while (fgets(fid)~=-1)
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
    end

%----------------all_stress_after_shear
if (shear==1):
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
end

%----------------all_stress_after_evaporation
if (evaporation==1)
    fid=fopen('all_stress_after_evaporation.txt');
    nLines=0;
    while (fgets(fid)~=-1)
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
end

