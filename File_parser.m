clc
clear
% close all

no_of_files = 10; %no of data files in the big appended file
index_start = 501; %starting index of the trial
index_end = index_start + no_of_files - 1; %ending index of the trial
v_air = 0.1; %ambient air velocity in the room
dist = 1.5; % distance between the TX and RX

propagation_time = zeros(no_of_files*1000, 1); %Initialize

file = ['C:\Fatih\Research\Molecular_Communication\MC_for_Airborne_Pathogen_Transmission\Turbulent_Flow_Modeling\CFD\Ansys_Fluent_Files\Data\outlet_v_air_', num2str(v_air), '_dist_',  num2str(dist), '_tr_new10.dpm'];
frm = '%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]';
fid = fopen(file); 
A = textscan( fid, frm, 'HeaderLines', 2, 'Delimiter', {'(('}, 'MultipleDelimsAsOne', 1);
fclose(fid);
propagation_time = A{12};
diff_t = diff(propagation_time);
figure; plot(diff_t);
thr = -0.2; %determine the threshold according to the derivative figure
ind_parse = [0; find(diff_t < thr); length(propagation_time)]; %find the last index of each trial file
A_d = cell2mat(A); 

i = 1;
for i_trial = index_start:index_end
    file_wr = ['C:\Fatih\Research\Molecular_Communication\MC_for_Airborne_Pathogen_Transmission\Turbulent_Flow_Modeling\CFD\Ansys_Fluent_Files\Data\outlet_v_air_', num2str(v_air), '_dist_',  num2str(dist), '_tr_', num2str(i_trial), '.dpm'];
    fileID = fopen(file_wr,'w');
    fprintf(fileID,'(outlet 13)\n(          x           y            z            u            v            w     diameter            t  parcel-mass         mass  n-in-parcel         time    flow-time)\n');
%     fprintf(fileID,'%6s %12s\n','x','exp(x)');    
    fprintf(fileID,'((%11.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e)\n', A_d(ind_parse(i)+1:ind_parse(i+1),:)');
    fclose(fileID);
    i = i + 1;
end

