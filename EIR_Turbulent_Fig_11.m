% End-to-end impulse response estimation and analysis
close all
clear
clc

tic
%% Parameters
no_of_files = 500; %no of data files for each ambient air velocity
v_air = [0.1 0.3 0.5]; %vambient air velocity in the room in the same direction with cough
dist = 1.5; % distance between the TX and RX
gamma = 0:1500:45000; % detection threshold
num_of_bins = 50; %vnumber of bins in histogram

options = statset('MaxIter',10000,'MaxFunEvals',20000); %options for mle estimation

%% Number of received particles
for i = 1:length(v_air)
    if i == 1
        [xi, no_of_rec_particles, particle_diameter_1, reception_time_1] = turbulent_stat_fun(no_of_files, v_air(i), dist);
        clearvars xi no_of_rec_particles y1 y2 y3 y3_m
    elseif i == 2
        [xi, no_of_rec_particles, particle_diameter_2, reception_time_2] = turbulent_stat_fun(no_of_files, v_air(i), dist);
        clearvars xi no_of_rec_particles y1 y2 y3 y3_m
    else
        [xi, no_of_rec_particles, particle_diameter_3, reception_time_3] = turbulent_stat_fun(no_of_files, v_air(i), dist);     
    end  

end

%% V_AIR = 0.1 m/s
% Histogram
h = figure;
hist_t1 = histogram(reception_time_1,100,'Normalization','pdf');
title(['v_{air} = ', num2str(v_air(1)), ' m/s' ]); grid on;
xlabel('Time (s)'); ylabel('Probability density'); hold on;
xlim([2 18]);
% MLE for pdf
pdf_mixture = @(reception_time_1, sigma_1, mu_1, sigma_2, mu_2, sigma_3, mu_3, w_1, w_2) w_1*pdf('Normal',reception_time_1,mu_1,sigma_1) +...
    w_2*pdf('Normal',reception_time_1,sigma_2,mu_2) + (1-w_1-w_2)*pdf('Normal',reception_time_1,sigma_3,mu_3);
t_01 = [0.18 4.5 0.2 12.65 0.5 13.4 0.05 0.55]; %initial guess for v_air=0.1 m/s
%      (1-w_1)*pdf('Stable',reception_time_1, alpha_2, beta_2, c_2, mu_2s); % + (1-w_1-w_2)*pdf('Normal',reception_time_1,mu_2,sigma_2);
% t_01 = [0.2 4.5 1 0.81 0.17 12.7 0.1]; %initial guess for v_air=0.1 m/s


% t_1 = mle(reception_time_1, 'pdf',pdf_mixture,'Start',t_01,'Options',options,'LowerBound',[0 4 0 0 0 0 0], 'UpperBound',[1 5 2 2 2 14 1]); % normal+stable
t_1 = mle(reception_time_1,'pdf',pdf_mixture,'Start',t_01,'Options',options, 'LowerBound',[0.15 4 0.15 12.6 0.4 13.4 0 0.55], 'UpperBound',[1 5 0.21 12.7 0.55 13.6 0.1 0.6]); % 3 normal
t_1(9) = 1 - t_1(7) - t_1(8);
% t_1(8) = 1 - t_1(7);
% t_1(6) = 1 - t_1(5);

pd1_1 = makedist('Normal','mu',t_1(2),'sigma',t_1(1));
pd1_2x = makedist('Normal','mu',t_1(4),'sigma',t_1(3));
pd1_2 = makedist('Normal','mu',t_1(6),'sigma',t_1(5));
% pd1_2 = makedist('Stable','alpha',t_1(3),'beta',t_1(4), 'gam', t_1(5), 'delta', t_1(6));

[f,ti] = ksdensity(reception_time_1);
ti_pos = ti(ti >= 0); %eliminate the negative values 
yt1_1 = t_1(7)*pdf(pd1_1,ti);
yt1_2x = t_1(8)*pdf(pd1_2x,ti);
yt1_2 = t_1(9)*pdf(pd1_2,ti);
yt1 = yt1_1 + yt1_2x + yt1_2;
hist_nv_t = hist_t1.BinCounts./(hist_t1.BinWidth*length(reception_time_1)); %calculate the exact normalized values in pdf histogram
for i_m = 1:length(hist_nv_t)
    i_ht(i_m) = mean(hist_t1.BinEdges(i_m:i_m+1));
    i_yt(i_m) = find(min(abs(ti_pos - i_ht(i_m))));
end
yt1_m = yt1(i_yt);
ms_error = mean((yt1_m - hist_nv_t).^2);
% plot(ti_pos,yt1_1, 'g.-', LineWidth=1.5);
% plot(ti_pos,yt1_2x, 'm.-', LineWidth=1.5);
% plot(ti_pos,yt1_2, 'k--', LineWidth=1.5);
plot(ti_pos,yt1, 'r-', LineWidth=1.5);
% legend('Histogram', '$w_1 \mathcal{N}$($\mu_1,\sigma_1^2$)',  '$w_2 \mathcal{N}$($\mu_2,\sigma_2^2$)',...
%            '$w_2 \mathcal{N}$($\mu_3,\sigma_3^2$)', '$w_1 \mathcal{N}$($\mu_1,\sigma_1^2$) + $w_2 \mathcal{N}$($\mu_2,\sigma_2^2$) + $w_2 \mathcal{N}$($\mu_3,\sigma_3^2$)', 'interpreter', 'latex');
legend('Histogram', '$w_1 \mathcal{N}$($\mu_1,\sigma_1^2$) + $w_2 \mathcal{N}$($\mu_2,\sigma_2^2$) + $w_2 \mathcal{N}$($\mu_3,\sigma_3^2$)', 'interpreter', 'latex');
title(['v_{air} = ', num2str(v_air(1)), ' m/s,', ' MSE = ', num2str(ms_error)]); grid on;
fontsize(gca, 14,'points');

%% V_AIR = 0.3 m/s
h = figure;
hist_t2 = histogram(reception_time_2,100,'Normalization','pdf');
title(['v_{air} = ', num2str(v_air(2)), ' m/s' ]); grid on;
xlabel('Time (s)'); ylabel('Probability density'); hold on;
% MLE for pdf
% pdf_mixture = @(reception_time_2, k_1, lambda_1) pdf('Weibull',reception_time_2,lambda_1, k_1);
% t_02 = [4.5 16]; %initial guess for v_air=0.3 m/s
% t_2 = mle(reception_time_2, 'pdf',pdf_mixture,'Start',t_02,'Options',options,'LowerBound',[0 0], 'UpperBound',[100 100]); % for v_air=0.3 m/s
% pd2 = makedist('Weibull','A',t_2(2),'B',t_2(1));

% Addition of two normal distributions
pdf_mixture = @(reception_time_2, sigma_1, mu_1, sigma_2, mu_2, w_1) w_1*pdf('Normal',reception_time_2, mu_1, sigma_1) +...
    (1-w_1)*pdf('Normal',reception_time_2, mu_2, sigma_2);
t_02 = [0.5 4.5 0.5 4.7 0.5]; %initial guess for v_air=0.3
t_2 = mle(reception_time_2, 'pdf',pdf_mixture,'Start',t_02,'Options',options, 'LowerBound', [0 0 0 0 0], 'UpperBound', [10 10 10 10 1]); % for v_air=0.3 m/s
t_2(6) = 1 - t_2(5); % w_2

pd2_1 = makedist('Normal','mu',t_2(2),'sigma',t_2(1));
pd2_2 = makedist('Normal','mu',t_2(4),'sigma',t_2(3));

% pdf_mixture = @(reception_time_2, sigma_1, mu_1) pdf('Normal',reception_time_2,mu_1, sigma_1);
% t_02 = [0.2 4.5]; %initial guess for v_air=0.3 m/s
% t_2 = mle(reception_time_2, 'pdf',pdf_mixture,'Start',t_02,'Options',options,'LowerBound',[0 0], 'UpperBound',[100 100]); % for v_air=0.3 m/s
% pd2 = makedist('Normal','mu',t_2(2),'sigma',t_2(1));

[f,ti] = ksdensity(reception_time_2);
ti_pos = ti(ti >= 0); %eliminate the negative values
yt2_1 = t_2(5)*pdf(pd2_1,ti);
yt2_2 = t_2(6)*pdf(pd2_2,ti);
yt2 = yt2_1 + yt2_2;
% yt2 = pdf(pd2,ti); % one normal distribution
hist_nv_t = hist_t2.BinCounts./(hist_t2.BinWidth*length(reception_time_2)); %calculate the exact normalized values in pdf histogram
for i_m = 1:length(hist_nv_t)
    i_ht(i_m) = mean(hist_t2.BinEdges(i_m:i_m+1));
    i_yt(i_m) = find(min(abs(ti_pos - i_ht(i_m))));
end
yt2_m = yt2(i_yt);
ms_error = mean((yt2_m - hist_nv_t).^2);
plot(ti,yt2, 'r-', LineWidth=1.5);
legend('Histogram', '$w_1 \mathcal{N}$($\mu_1,\sigma_1^2$) + $w_2 \mathcal{N}$($\mu_2,\sigma_2^2$)', 'interpreter', 'latex');
title(['v_{air} = ', num2str(v_air(2)), ' m/s,', ' MSE = ', num2str(ms_error)]); grid on;
fontsize(gca, 14,'points');

%% V_AIR = 0.5 m/s
h = figure;
hist_t3 = histogram(reception_time_3,100,'Normalization','pdf');
title(['v_{air} = ', num2str(v_air(3)), ' m/s' ]); grid on;
xlabel('Time (s)'); ylabel('Probability density');hold on;
% MLE for pdf
% pdf_mixture = @(reception_time_3, alpha_1, beta_1, c_1, mu_1) pdf('Stable',reception_time_3, alpha_1, beta_1, c_1, mu_1);
% t_03 = [1.47 0.99 0.03 2.79]; %initial guess for v_air=0.3 m/s
% t_3 = mle(reception_time_3, 'pdf',pdf_mixture,'Start',t_03,'Options',options); % for v_air=0.3 m/s
% pd3 = makedist('Stable','alpha',t_3(1),'beta',t_3(2),'gam',t_3(3),'delta',t_3(4));

% pdf_3 = fitdist(reception_time_3,'InverseGaussian');
pdf_3 = fitdist(reception_time_3,'Stable');
[f,ti] = ksdensity(reception_time_3);
ti_pos = ti(ti >= 0); %eliminate the negative values 
yt3 = pdf(pdf_3,ti);
hist_nv_t = hist_t3.BinCounts./(hist_t3.BinWidth*length(reception_time_3)); %calculate the exact normalized values in pdf histogram
for i_m = 1:length(hist_nv_t)
    i_ht(i_m) = mean(hist_t3.BinEdges(i_m:i_m+1));
    i_yt(i_m) = find(min(abs(ti_pos - i_ht(i_m))));
end
yt3_m = yt3(i_yt);
ms_error = mean((yt3_m - hist_nv_t).^2);
plot(ti,yt3, 'r-', LineWidth=1.5);
legend('Histogram', '$\mathcal{S}$($\alpha_1,\beta_1, c_1, m_1$)', 'interpreter', 'latex');
title(['v_{air} = ', num2str(v_air(3)), ' m/s,', ' MSE = ', num2str(ms_error)]); grid on;
fontsize(gca, 14,'points');

%% V_AIR - AVG Velocity Comparison
v_avg(1,:) = dist./[t_1(2) t_1(4) t_1(6)];
v_avg(2,:) = dist/t_2(2);
v_avg(3,:) = dist/pdf_3.delta;

toc
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,sprintf('EIR_pdf_plot_v_air_%.1f.pdf', v_air(1)),'-dpdf','-r0')%save as pdf 
% savefig(h,sprintf('SIR_TLW_SIR_plot_dinf_%.1f_alpha_r_%.1f_MC_%d.fig',d_inf,alpha_r,MC)); %save the figure file
