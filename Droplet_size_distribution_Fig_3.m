close all
clear
clc

%% Xie, 2009
d_cough_X = 10^-6.*[10 15 20 30 35 40 45 50 75  100 150 200 250 300 350 400 450 500 1000 1500] ; % Droplet diameters (Xie,2009, Avg numbers in Table 2)
N_cough_X =        [11 4  6  27 12 27 48 23 179 176 150 56  27  4   23  2   4   3   19   2] ; % Droplet numbers (Xie,2009, Avg numbers of four subjects in Table 2)

data_cough_X = [];
d_lim = find(d_cough_X == 300*10^-6);
for i = 1:d_lim
    data_cough_X = [data_cough_X ones(1,N_cough_X(i))*d_cough_X(i)];
end

% h = figure;
% bar(d_cough_X, N_cough_X);
% 
% h2 = figure;
% histogram(data_cough_X, length(d_cough_X), 'Normalization','pdf');
pd = fitdist(data_cough_X', 'Weibull')
k_s = pd.B; %shape parameter - spread of the particle size
lambda_s = pd.A; %scale parameter - mean particle size

h = figure;
hf = histfit(data_cough_X, length(d_cough_X)-9, 'weibull');
title('Droplet Size Distribution');
xlabel('Droplet Size (m)'); ylabel('Number of Droplets');
label1 = 'Histogram of measurements'; label2 = ['Fitted Weibull distribution ($\lambda_d = $', num2str(lambda_s), ', $k_d = $', num2str(k_s), ')'];
legend([hf(1) hf(2)],{label1, label2},'Interpreter','latex'); grid on;
ylim([0 200]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'Droplet_size_distr','-dpdf','-r0')%save as pdf 



% % hold on
% h4 = figure;
% % lambda_s2 = 80*10^-6; k_s2 = 8; 
% y = wblpdf(hf(2).XData, lambda_s, k_s);
% plot(hf(2).XData, y/sum(y));






% %% Duguid, 1946
% d_cough_D = 10^-6.*[2 4 8 16 24 32 40 50 75 100 125 150 200 250 500 1000 2000]; %Droplet diameters (Duguid,1946)
% N_cough_D = [50 290 970 1600 870 420 240 110 140 85 48 38 35 29 34 12 2]; %Droplet numbers (Duguid,1946)
% 
% data_cough_D = [];
% d_lim = find(d_cough_D == 2000*10^-6);
% for i = 1:d_lim
%     data_cough_D = [data_cough_D ones(1,N_cough_D(i))*d_cough_D(i)];
% end
% 
% h = figure;
% bar(d_cough_D, N_cough_D);
% 
% h2 = figure;
% histogram(data_cough_D, length(d_cough_D), 'Normalization','probability');
% 
% h3 = figure;
% histfit(data_cough_D, 10*length(d_cough_D), 'weibull');
% 
% pd = fitdist(data_cough_D', 'Weibull')
% k_s = pd.B; %shape parameter - spread of the particle size
% lambda_s = pd.A; %scale parameter - mean particle size

