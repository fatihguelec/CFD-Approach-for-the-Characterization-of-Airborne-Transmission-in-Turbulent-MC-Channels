function [xi, no_of_rec_particles, particle_diameter, reception_time] = turbulent_stat_fun(no_of_files, v_air, dist)


    j = 1;
    for i_trial = 1:no_of_files
    %     file = sprintf('C:\Fatih\Research\Molecular_Communication\MC_for_Airborne_Pathogen_Transmission\Turbulent_Flow_Modeling\CFD\Ansys_Fluent_Files\Data\outlet_v_air_%.1f_dist_%.1f_tr_%d.dpm', v_air, dist, i_trial);             
        file = ['Data/outlet_v_air_', num2str(v_air), '_dist_',  num2str(dist), '_tr_', num2str(i_trial), '.dpm'];
%         file = ['G:\Other computers\Dell_Eckford_Lab\Research\Molecular_Communication\MC_for_Airborne_Pathogen_Transmission\Turbulent_Flow_Modeling\CFD\Ansys_Fluent_Files\Data\outlet_v_air_', num2str(v_air), '_dist_',  num2str(dist), '_tr_', num2str(i_trial), '.dpm'];
        frm = '%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]';
        fid = fopen(file); 
        A = textscan( fid, frm, 'HeaderLines', 2, 'Delimiter', {'(('}, 'MultipleDelimsAsOne', 1);
        fclose(fid);
        no_of_rec_particles(i_trial,:) = size(A{7},1);
        temp = j + size(A{7},1);
        particle_diameter(j:temp-1) = A{7};
        reception_time(j:temp-1)= A{12};
        j = temp;   
    end
    
    %% Kernel smoothing function estimate of received particles - only xi is used
    [f,xi] = ksdensity(no_of_rec_particles); %pdf estimate is f and xi is the equally-spaced points that cover the range of the data
    xi = xi(xi >= 0); %eliminate the values below 0
    
    no_of_rec_particles = no_of_rec_particles';
    particle_diameter = particle_diameter';
    reception_time = reception_time';
end


