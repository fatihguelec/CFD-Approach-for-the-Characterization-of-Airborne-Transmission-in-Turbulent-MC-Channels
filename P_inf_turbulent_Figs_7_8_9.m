close all
clear
clc

%% Parameters
no_of_files = 500; %vno of data files for each ambient air velocity
v_air = [0.1 0.3 0.5]; %vambient air velocity in the room in the same direction with cough
dist = 1.5; % distance between the TX and RX
gamma = 0:1500:45000; % detection threshold
num_of_bins = 50; %vnumber of bins in histogram

options = statset('MaxIter',10000,'MaxFunEvals',20000); %options for mle estimation

%% Number of received particles
for i = 1:length(v_air)
    if i == 1
        [xi, no_of_rec_particles, particle_diameter_1, reception_time_1] = turbulent_stat_fun(no_of_files, v_air(i), dist);
        % MAXIMUM LIKELIHOOD ESTIMATION
        %  Two Gaussian = w_1*N(mu_1,sigma_1) + (1-w_1)*N(mu_2,sigma_2)
        pdf_mixture = @(no_of_rec_particles, sigma_1, mu_1, sigma_2, mu_2, w_1) w_1*normpdf(no_of_rec_particles,mu_1, sigma_1) +...
            (1-w_1)*normpdf(no_of_rec_particles,mu_2, sigma_2);
        x_01 = [500 6500 500 9000 0.5]; %initial guess for v_air=0.1 m/s

        %             (1-w_1)*pdf('Uniform', no_of_rec_particles, a_2,  b_2); 
%         x_01 = [500 2500 0 10000 0.5]; %initial guess for v_air=0.1 m/s, x_0 = [sigma_1 mu_1 w_1]
        
        x_1 = mle(no_of_rec_particles, 'pdf',pdf_mixture,'Start',x_01,'Options',options, 'LowerBound', [0 0 0 0 0], 'UpperBound', [4000 10000 1000 100000 1]); % for v_air=0.1 m/s
        x_1(6) = 1-x_1(5); % for v_air=0.1 m/s [sigma_1 mu_1 sigma_2 mu_2  w_1 w_2]
%         x_1(5) = x_1(3); x_1(3) = min(no_of_rec_particles); x_1(4) = max(no_of_rec_particles); x_1(6) = 1-x_1(5); % for v_air=0.1 m/s [sigma_1, mu_1, a_2 b_2 w_1 w_2]


        % Generate the Parametric functions
        pd1 = makedist('Normal','mu',x_1(2),'sigma',x_1(1));
        pd2 = makedist('Normal','mu',x_1(4),'sigma',x_1(3));

        % Get the probability values
        y1 = pdf(pd1,xi)*x_1(5); % x(5) is the weight of pdf1
        y2 = pdf(pd2,xi)*x_1(6); % x(6) is the weight of pdf2
        y3 = y1 + y2; %final estimated pdf

        %Histogram
        h1 = figure;
        hist = histogram(no_of_rec_particles,num_of_bins,'Normalization','pdf');
       hold on;
        % Mean Squared Error between estimated pdf and histogram data
        hist_nv = hist.BinCounts./(hist.BinWidth*length(no_of_rec_particles)); %calculate the exact normalized values in pdf histogram
        for i_m = 1:length(hist_nv)
            i_h(i_m) = mean(hist.BinEdges(i_m:i_m+1));
            i_y(i_m) = find(min(abs(xi - i_h(i_m))));
        end
        y3_m = y3(i_y);
        ms_error = mean((y3_m - hist_nv).^2);

%         plot(xi,y1, 'g.-', LineWidth=1.5); %first pdf
%         plot(xi,y2, 'm--', LineWidth=1.5); %second pdf
        plot(xi,y3, 'r-', 'LineWidth', 1.5); %final pdf
        legend('Histogram', '$w_1 \mathcal{N}$($\mu_1,\sigma_1^2$) + $w_2 \mathcal{N}$($\mu_2,\sigma_2^2$)', 'interpreter', 'latex');
        title(['v_{air} = ', num2str(v_air(i)), ' m/s,', ' MSE = ', num2str(ms_error)]); grid on;
        xlabel('Number of received particles'); ylabel('Probability density'); ylim([0 6e-04]);
        %fontsize(gca, 14,'points');

        clearvars xi no_of_rec_particles y1 y2 y3 y3_m
    elseif i == 2
        [xi, no_of_rec_particles, particle_diameter_2, reception_time_2] = turbulent_stat_fun(no_of_files, v_air(i), dist);
        % MAXIMUM LIKELIHOOD ESTIMATION
        pdf_mixture = @(no_of_rec_particles,sigma_1, mu_1, sigma_2, mu_2, k_2, lambda_2, w_1, w_2) w_1*normpdf(no_of_rec_particles,mu_1, sigma_1) +... 
            w_2*normpdf(no_of_rec_particles,mu_2, sigma_2) + (1-w_1-w_2)*pdf('Weibull',no_of_rec_particles,lambda_2, k_2);
        x_02 = [500 6500 500 8000 50 15000 0.25 0.25]; %initial guess for v_air=0.3 m/s [sigma_1, mu_1, sigma_2, mu_2, k_2, lambda_2, w_1, w_2]
        
        x_2 = mle(no_of_rec_particles,'pdf',pdf_mixture,'Start',x_02,'Options',options);
        x_2(9) = 1 - x_2(7) - x_2(8); % w_3 for 2N + 1W, x_2 = [sigma_1, mu_1, sigma_2, mu_2, k_2, lambda_2, w_1, w_2, w_3]

        % Generate the Parametric functions
        pd1 = makedist('Normal','mu',x_2(2),'sigma',x_2(1));
        pd2 = makedist('Normal','mu',x_2(4),'sigma',x_2(3));
        pd3 = makedist('Weibull','A',x_2(6),'B',x_2(5));

        % Get the probability values
        y1 = pdf(pd1,xi)*x_2(7) + pdf(pd2,xi)*x_2(8); % x(5) is the weight of pdf1
        y2 = pdf(pd3,xi)*x_2(9); % x(6) is the weight of pdf2
        y3 = y1 + y2; %final estimated pdf

        %Histogram
        h2 = figure;
        hist = histogram(no_of_rec_particles,num_of_bins,'Normalization','pdf');
        xlabel('Number of received particles'); ylabel('Probability density');
        hold on;

        % Mean Squared Error between estimated pdf and histogram data
        hist_nv = hist.BinCounts./(hist.BinWidth*length(no_of_rec_particles)); %calculate the exact normalized values in pdf histogram
        for i_m = 1:length(hist_nv)
            i_h(i_m) = mean(hist.BinEdges(i_m:i_m+1));
            i_y(i_m) = find(min(abs(xi - i_h(i_m))));
        end
        y3_m = y3(i_y);
        ms_error = mean((y3_m - hist_nv).^2);

%         plot(xi,y1, 'g.-', LineWidth=1.5); %first pdf
%         plot(xi,y2, 'm--', LineWidth=1.5); %second pdf
        plot(xi,y3, 'r-', 'LineWidth', 1.5); %final pdf
        legend('Histogram',  '$w_1 \mathcal{N}$($\mu_1,\sigma_1^2$) + $w_2 \mathcal{N}$($\mu_2,\sigma_2^2$) + $w_3 \mathcal{W}$($k_3,\lambda_3$)', 'interpreter', 'latex');
        title(['v_{air} = ', num2str(v_air(i)), ' m/s,', ' MSE = ', num2str(ms_error)]);
        xlim([0.05 20000]); grid on;
        %fontsize(gca, 14,'points');

        clearvars xi no_of_rec_particles y1 y2 y3 y3_m
    else
        [xi, no_of_rec_particles, particle_diameter_3, reception_time_3] = turbulent_stat_fun(no_of_files, v_air(i), dist);
        % Weibull + Gaussian = w_1*W(k_1,lambda_1) + (1-w_1)*N(mu_2,sigma_2)
        pdf_mixture = @(no_of_rec_particles,k_1, lambda_1, sigma_2, mu_2, w_1) w_1*pdf('Weibull',no_of_rec_particles,lambda_1, k_1)+... 
            (1-w_1)*normpdf(no_of_rec_particles,mu_2, sigma_2);
        x_03 = [15 15000 500 33000 0.75]; %initial guess for v_air=0.5 m/s, x_0 = [k_1 lambda_1 sigma_2 mu_2 w_1]

        x_3 = mle(no_of_rec_particles,'pdf',pdf_mixture,'Start',x_03,'Options',options); % for v_air=0.5 m/s 
        x_3(6) = 1 - x_3(5); % x = [k_1 lambda_1 sigma_2 mu_2 w_1 w_2] where w_2 = 1 - w_1

        pd1 = makedist('Weibull','A',x_3(2),'B',x_3(1));
        pd2 = makedist('Normal','mu',x_3(4),'sigma',x_3(3));

        % Get the probability values
        y1 = pdf(pd1,xi)*x_3(5); % x(5) is the weight of pdf1
        y2 = pdf(pd2,xi)*x_3(6); % x(6) is the weight of pdf2
        y3 = y1 + y2; %final estimated pdf     

        %Histogram
        h3 = figure;
        hist = histogram(no_of_rec_particles,num_of_bins,'Normalization','pdf');
        xlabel('Number of received particles'); ylabel('Probability density');
        hold on;

        % Mean Squared Error between estimated pdf and histogram data
        hist_nv = hist.BinCounts./(hist.BinWidth*length(no_of_rec_particles)); %calculate the exact normalized values in pdf histogram
        for i_m = 1:length(hist_nv)
            i_h(i_m) = mean(hist.BinEdges(i_m:i_m+1));
            i_y(i_m) = find(min(abs(xi - i_h(i_m))));
        end
        y3_m = y3(i_y);
        ms_error = mean((y3_m - hist_nv).^2);

%         plot(xi,y1, 'g.-', LineWidth=1.5); %first pdf
%         plot(xi,y2, 'm--', LineWidth=1.5); %second pdf
        plot(xi,y3, 'r-', 'LineWidth', 1.5); %final pdf
        legend('Histogram', '$w_1 \mathcal{W}$($k_1,\lambda_1$) + $w_2 \mathcal{N}$($\mu_2,\sigma_2^2$)', 'interpreter', 'latex');
        title(['v_{air} = ', num2str(v_air(i)), ' m/s,', ' MSE = ', num2str(ms_error)]); grid on;
        %fontsize(gca, 14,'points');
%         xlim([0.05 20000]);         
    end    

end


%% Probability of infection (P_inf) 
% P_inf(1,:) = x_1(5).*qfunc( (gamma - x_1(2)) ./ x_1(1) ) + x_1(6).* ((x_1(4) - gamma)./(x_1(4) - x_1(3)) ); %for v_air=0.1 m/s
P_inf(1,:) = x_1(5).*qfunc( (gamma - x_1(2)) ./ x_1(1) ) + x_1(6).* qfunc( (gamma - x_1(4)) ./ x_1(3) ); %for v_air=0.1 m/s
P_inf(2,:) = x_2(7).*qfunc( (gamma - x_2(2)) ./ x_2(1) ) + x_2(8).*qfunc( (gamma - x_2(4)) ./ x_2(3)) + x_2(9).*exp(-(gamma./x_2(6)).^x_2(5) ); %for v_air=0.3 m/s
P_inf(3,:) = x_3(5).*exp(-(gamma./x_3(2)).^x_3(1) ) + x_3(6) .* qfunc( (gamma - x_3(4)) ./ x_3(3) ); %for v_air=0.5 m/s

h = figure;
plot(gamma, P_inf(1,:), 'b-o', LineWidth=1.25); hold on;
plot(gamma, P_inf(2,:), 'r-.*', LineWidth=1.25);
plot(gamma, P_inf(3,:), 'g--d', LineWidth=1.25);
xlabel('Threshold'); ylabel('Probability of infection');
legend('v_{air} = 0.1 m/s', 'v_{air} = 0.3 m/s', 'v_{air} = 0.5 m/s'); grid on;
axis([-500 45000 -0.05 1.05]); 


%% Particle Diameters
h = figure;
histogram(particle_diameter_1,50,'Normalization','probability');
title(['v_{air} = ', num2str(v_air(1)), ' m/s' ]); grid on;
xlabel('Diameter of received particles (m)'); ylabel('Probability'); 
% fontsize(gca, 14,'points');
axis([-1e-6 4.5e-5 0 0.18]); 

h = figure;
histogram(particle_diameter_2,50,'Normalization','probability');
title(['v_{air} = ', num2str(v_air(2)), ' m/s' ]); grid on;
xlabel('Diameter of received particles (m)'); ylabel('Probability'); 
% fontsize(gca, 14,'points');
% axis([-1e-6 4.5e-5 0 0.18]); 
% 
% ax = gca; ax.FontSize = 13;
% magnifyOnFigure;

h = figure;
histogram(particle_diameter_3,50,'Normalization','probability');
title(['v_{air} = ', num2str(v_air(3)), ' m/s' ]); grid on;
xlabel('Diameter of received particles (m)'); ylabel('Probability'); 
% fontsize(gca, 14,'points');
axis([-1e-6 4.5e-5 0 0.18]); 

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,sprintf('P_diameter_v_air_%.1f_mag.pdf', v_air(3)),'-dpdf','-r0')%save as pdf 
% % print(h,sprintf('P_inf_all.pdf'),'-dpdf','-r0')%save as pdf 
% savefig(h,sprintf('SIR_TLW_SIR_plot_dinf_%.1f_alpha_r_%.1f_MC_%d.fig',d_inf,alpha_r,MC)); %save the figure file


set(h3,'Units','Inches');
pos = get(h3,'Position');
set(h3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h3,sprintf('pdf_v_air_0.5.pdf'),'-dpdf','-r0')%save as pdf 
% savefig(h,sprintf('SIR_TLW_SIR_plot_dinf_%.1f_alpha_r_%.1f_MC_%d.fig',d_inf,alpha_r,MC)); %save the figure file



