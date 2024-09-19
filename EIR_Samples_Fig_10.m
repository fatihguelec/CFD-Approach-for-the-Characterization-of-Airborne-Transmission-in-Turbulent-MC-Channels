close all
clear
clc

%% Parameters
dist = 1.5; % distance between the TX and RX
no_of_files = [68 109 203]; v_air = 0.1; % data file # for v_air = 0.1 - Fig 10-a
% no_of_files = [38 103 245]; v_air = 0.3;  % data file # for v_air = 0.3
% - Fig 10-b
% no_of_files = [1 10 166]; v_air = 0.5;  % data file # for v_air = 0.5 - Fig 10-c

%% End-to-end impulse response (EIR) Plots
% for i_trial = 1:no_of_files
for i_trial = no_of_files
    file = ['Data/outlet_v_air_', num2str(v_air), '_dist_',  num2str(dist), '_tr_', num2str(i_trial), '.dpm'];
    frm = '%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]';
    fid = fopen(file); 
    A = textscan( fid, frm, 'HeaderLines', 2, 'Delimiter', {'(('}, 'MultipleDelimsAsOne', 1);
    fclose(fid);

    edges_x = 2:0.1:18; 
    [hist, edges] = histcounts(A{12}, edges_x);
    h = figure; 
    plot(edges, [0 hist], 'b-', 'LineWidth', 1.5); %v_air = 0.1
%     plot(edges, [0 hist], 'r-.', 'LineWidth', 1.5); %v_air = 0.3
    % plot(edges, [0 hist], 'g-', 'LineWidth', 1.5); %v_air = 0.5
    grid on; 
    xlabel('Time (s)'); ylabel('h(t)');
    title(['v_{air} = ', num2str(v_air), ' m/s']);
    fontsize(gca, 14,'points');
%     xlim([3.4 18]); 

    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(h,sprintf('EIR_sample_plot_v_air_%.1f_trial_%d.pdf', v_air, i_trial),'-dpdf','-r0')%save as pdf 
end


