close all
clear
clc

%% Lindsley et.al., 2012
load('Aerosol_Avg.mat')
A = table2array(Aerosol_Avg);
d_cough_L = 10^-6.*A(:,1); % Aerosol diameters
N_cough_L = ceil(A(:,2)); % Aerosol numbers 
N_cough_total = sum(N_cough_L); % Avg Total number of aerosols

data_cough_L = [];
% d_lim = find(d_cough_L == 300*10^-6);
for i = 1:size(A,1)
    data_cough_L = [data_cough_L ones(1,N_cough_L(i))*d_cough_L(i)];
end

h = figure;
bar(d_cough_L, N_cough_L);

h2 = figure;
histogram(data_cough_L, length(d_cough_L), 'Normalization','pdf');

pd = fitdist(data_cough_L', 'Weibull')
k_s = pd.B; %shape parameter - spread of the particle size
lambda_s = pd.A; %scale parameter - mean particle size

h = figure;
hf = histfit(data_cough_L, length(d_cough_L)-8, 'weibull');
title('Aerosol Size Distribution');
xlabel('Aerosol Size (m)'); ylabel('Number of Aerosols');
label1 = 'Histogram of measurements'; label2 = ['Fitted Weibull distribution ($\lambda_a = $', num2str(lambda_s), ', $k_a = $', num2str(k_s), ')'];
legend([hf(1) hf(2)],{label1, label2},'Interpreter','latex'); grid on;
% ylim([0 200]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'Aerosol_size_distr','-dpdf','-r0')%save as pdf 


% Cough Mass Flow Rate Calculation
rho_a = 1.704; % Aerosol density (kg/l)
Q_v = [3.906 6.949 7.558 4.453 5.329 5.771 5.386 4.480 7.337 3.871 9.086 7.638 5.824 4.792 4.465 7.025 6.599 6.411 4.947 4.474 6.517 6.432 4.727]; % Cough peak volumetric flow rate for 23 patients, Session 1 , Cough 1 (l/s) - Suppl-Inf-2
Q = rho_a.*Q_v; % Mass flow rate (kg/s)
Q_avg = mean(Q); % Average mass flow rate to be used in Fluent





