function CROC_simulation_for_PD()
%% Description: 
% This function simulates the complementary receiver operating characteristics (CROC)
% of the investigated power detector (PD) for Q=1   
% Note that the CROC is going to be plotted via the Monte-Carlo (MC) simulation
% and the closed-form expressions whose results are going to overlap with the
% the MC simulation results 

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Matlab code: 

m_1=2;
% The SOI fading severity parameter set as per [1, Table 1].
m_2=2;
% Fading severity parameter of an RFI set as per [1, Table 1].
m=[m_1, m_2]; 
% A vector of the fading severity parameters of the SOI and RFI 
av_snr=0;
% The average signal-to-noise ratio (SNR) in dB
av_inr=[5,10]; 
% The average interference-to-noise ratios (INRs) for an RFI in dB
lambda=5:-0.25:1.5;
% Decision thresholds 
N=10000; 
% The number of intercepted samples (per a realization) set as per [1, Table 1]. 
num_iter=50000;
% The number of Monte-Carlo iterations 
PM_closed_form=zeros(length(av_inr),length(lambda)); 
% Initialization of the closed-form results for the probability of miss
% (P_m) evaluated for different average INRs and decision thresholds. 
% for P_d being the probability of detection, the probability of miss is
% computed as P_m=1-P_d 
FA_closed_form=zeros(length(av_inr),length(lambda)); 
% Initialization of the closed-form results for the probability of false
% alarm evaluated for different average INRs and decision thresholds. 
PM_Monte_Carlo=zeros(length(av_inr),length(lambda)); 
% Initialization of the Monte-Carlo simulation results for the probability of miss
% (P_m) evaluated for different average INRs and decision thresholds. 
% Note again: P_m=1-P_d  
FA_Monte_Carlo=zeros(length(av_inr),length(lambda)); 
% Initialization of the Monte-Carlo simulation results for the probability of false
% alarm evaluated for different average INRs and decision thresholds. 

for jj=1:length(av_inr)
    for ii=1:length(lambda)
        PM_closed_form(jj,ii)=1-PD_analytical_prob_RFI_detection(lambda(ii), m, av_snr, av_inr(jj)); 
        % Results obtained via the derived closed-form expression (via [1, eq. (9)] )
        FA_closed_form(jj,ii)=PD_analytical_prob_false_alarm(lambda(ii), m_1, av_snr); 
        % Results obtained via the derived closed-form expressions (via [1, eq. (12)] ) 
        PM_Monte_Carlo(jj,ii)=1-PD_MC_prob_RFI_detection(lambda(ii), m, av_snr, av_inr(jj), N, num_iter); 
        % Results obtained via the Monte-Carlo simulation   
        FA_Monte_Carlo(jj,ii)=PD_MC_prob_false_alarm(lambda(ii), m_1, av_snr, N, num_iter); 
        % Results obtained via the Monte-Carlo simulation   
    end 
end

%% CROC plots for different average INRs

loglog(FA_closed_form(1,:),PM_closed_form(1,:), '-ro', 'LineWidth',1); 
xlabel('\textbf{Probability of false alarm}', 'interpreter', 'latex','FontSize',12);
ylabel('\textbf{Probability of miss}', 'interpreter', 'latex','FontSize',12); 
xlim([0.005,0.733]); 
% Ranges of x-axis (FAR axis) for the aforementioned ranges of lambda) 
ylim([1e-3,1]); 
grid on; 
hold on;  
loglog(FA_Monte_Carlo(1,:),PM_Monte_Carlo(1,:), '--kx', 'LineWidth',1); 
 % Monte-Carlo
loglog(FA_closed_form(2,:),PM_closed_form(2,:), '-gd', 'LineWidth',1); 
 % Closed-form
loglog(FA_Monte_Carlo(2,:),PM_Monte_Carlo(2,:), '--bv', 'LineWidth',1); 
 % Closed-form 
h=legend('$\textbf{Closed-form (9) \& (12): }\bar{\gamma}_{inr}^1=\textbf{5 dB}$','$\textbf{Monte-Carlo: }\bar{\gamma}_{inr}^1=\textbf{5 dB}$', .......
    '$\textbf{Closed-form (9) \& (12): }\bar{\gamma}_{inr}^1=\textbf{10 dB}$', '$\textbf{Monte-Carlo:: }\bar{\gamma}_{inr}^1=\textbf{10 dB}$');  
set(h,'interpreter','latex','FontSize',12); 

end

         