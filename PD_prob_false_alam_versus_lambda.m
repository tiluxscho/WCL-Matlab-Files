function PD_prob_false_alam_versus_lambda()
%% Description: 
% This function plots the probability of false alarm obtained
% from the Monte-Carlo simulation and the closed-form expression. 
% Note that the function is written for the investigated power detector (PD)

%% Author: Tilahun M. Getu  

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018 

%% Matlab code:   

m_1=2;
% The SOI fading severity parameter set as per [1, Table 1].
av_snr=[0,5]; 
% The average signal-to-noise ratios (SNRs) in dB.
lambda=1:10;
% Decision thresholds.  
num_iter=10000; 
% The number of Monte-Carlo iterations set as per [1, Table 1]. 
N=10000; 
% The number of intercepted samples (per a realization) set as per [1, Table 1].  
FA_simulation=zeros(length(av_snr),length(lambda)); 
% Initialization of the Monte-Carlo simulation results to be evaluated for
% different average SNRs and decision thresholds
FA_closed_form=zeros(length(av_snr),length(lambda)); 
% Initialization of the closed-form results to be evaluated for
% different average SNRs and decision thresholds
for jj=1:length(av_snr)
    for ii=1:length(lambda)
        FA_simulation(jj,ii)=PD_MC_prob_false_alarm(lambda(ii), m_1, av_snr(jj), N, num_iter); 
        % Results obtained via the Monte-Carlo simulation 
        FA_closed_form(jj,ii)=PD_analytical_prob_false_alarm(lambda(ii), m_1, av_snr(jj)); 
        % Results obtained via the derived closed-form expressions [1, eqs. (12)] 
    end 
end

%% Plots
plot(lambda,FA_simulation(1,:), '-ro', 'LineWidth',1); 
xlabel('$\lambda$', 'interpreter', 'latex','FontSize',12);
ylabel('\textbf{Probability of false alarm}', 'interpreter', 'latex','FontSize',12); 
xlim([1,10]);
ylim([0, 1]); 
grid on; 
hold on;  
plot(lambda,FA_closed_form(1,:), '--b*', 'LineWidth',1); 
plot(lambda,FA_simulation(2,:), '-gs', 'LineWidth',1); 
plot(lambda,FA_closed_form(2,:), '--kd', 'LineWidth',1); 
h=legend('$\textbf{Monte-Carlo: }\bar{\gamma}_{snr}=\textbf{0 dB}$', '$\textbf{Closed-form (12): }\bar{\gamma}_{snr}=\textbf{0 dB}$',..........
    '$\textbf{Monte-Carlo: }\bar{\gamma}_{snr}=\textbf{5 dB}$', '$\textbf{Closed-form (12): }\bar{\gamma}_{snr}=\textbf{5 dB} $');  
set(h,'interpreter','latex','FontSize',12); 


end

