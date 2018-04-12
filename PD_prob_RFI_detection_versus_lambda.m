function PD_prob_RFI_detection_versus_lambda()
%% Description: 
% This function plots the probability of RFI detection obtained
% from the Monte-Carlo simulation and the closed-form expressions. 
% Note that this function plots the probability of detection for Q=1 and Q=2. 
% Note also that this function is written for the investigated power detector (PD)

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Matlab code:  

m_1=2;
% The SOI fading severity parameter set as per [1, Table 1].
m_2=2;
% Fading severity parameter of the first RFI set as per [1, Table 1].
m_3=2; 
% Fading severity parameter of the second RFI set as per [1, Table 1].
av_snr=0;
% The average signal-to-noise ratio (SNR) in dB
lambda=1:10;
% Decision thresholds 
av_inr_1=[8,10]; 
% The average interference-to-noise ratios (INRs) for the first RFI in dB 
av_inr_2=7; 
% The average INR for the second RFI in dB
num_iter=10000; 
% The number of Monte-Carlo iterations set as per [1, Table 1]. 
N=10000; 
% The number of intercepted samples (per a realization) set as per [1, Table 1]. 

%% The single interferer (Q=1) case  
m=[m_1, m_2]; 
% A vector of the fading severity parameters of the SOI and the first RFI
% Note that m is set as per [1, Table 1].  
av_inr=av_inr_1; 
% The average INR of the first RFI in dB  
PD_simulation=zeros(length(av_inr),length(lambda)); 
% Initialization of the Monte-Carlo simulation results to be evaluated for
% different average INRs and decision thresholds  
PD_closed_form=zeros(length(av_inr),length(lambda)); 
% Initialization of the closed-form results to be evaluated for
% different average INRs and decision thresholds  

for jj=1:length(av_inr)
    for ii=1:length(lambda)
        PD_simulation(jj,ii)=PD_MC_prob_RFI_detection(lambda(ii), m, av_snr, av_inr(jj), N, num_iter); 
        % Results obtained via the Monte-Carlo simulation  
        PD_closed_form(jj,ii)=PD_analytical_prob_RFI_detection(lambda(ii), m, av_snr, av_inr(jj)); 
        % Results obtained via the derived closed-form expressions [1, eqs. (9) and (10)]   
    end 
end

%% The double interferers (Q=2) case  
m=[m_1, m_2, m_3]; 
% A vector of the fading severity parameters of the SOI and the two RFIs
% Note that m is set as per [1, Table 1]. 
av_inr=[av_inr_1(1),av_inr_2(1)];
% The average INR (in dB) of the first and the second RFIs
PD_simulation_mult_Q=zeros(1,length(lambda)); 
% Initialization of the Monte-Carlo simulation results to be evaluated for
% different average INRs and decision thresholds
PD_closed_form_mult_Q=zeros(1,length(lambda)); 
% Initialization of the closed-form results to be evaluated for
% different average INRs and decision thresholds
for ii=1:length(lambda)
    PD_simulation_mult_Q(1,ii)=PD_MC_prob_RFI_detection(lambda(ii), m, av_snr, av_inr, N, num_iter); 
    % Results obtained via the Monte-Carlo simulation 
    PD_closed_form_mult_Q(1,ii)=PD_analytical_prob_RFI_detection(lambda(ii), m, av_snr, av_inr); 
    % Results obtained via the derived closed-form expressions [1, eqs. (9) and (10)] 
end

%% Plots for Q=1
plot(lambda,PD_simulation(1,:), '-ro', 'LineWidth',1); 
xlabel('$\lambda$', 'interpreter', 'latex','FontSize',12);
ylabel('\textbf{Probability of RFI detection}', 'interpreter', 'latex','FontSize',12); 
xlim([1,10]);
ylim([0, 1]); 
grid on; 
hold on; 
plot(lambda,PD_closed_form(1,:), '--b*', 'LineWidth',1); 
plot(lambda,PD_simulation(2,:), '-gs', 'LineWidth',1); 
plot(lambda,PD_closed_form(2,:), '--kd', 'LineWidth',1); 

%% Plots for Q=1
plot(lambda,PD_simulation_mult_Q(1,:), '-.cx', 'LineWidth',1); 
plot(lambda,PD_closed_form_mult_Q(1,:), '-.mv', 'LineWidth',1); 

h=legend('$\textbf{Monte-Carlo: }\bar{\gamma}_{inr}^1=\textbf{8 dB}$', '$\textbf{Closed-form (9): }\bar{\gamma}_{inr}^1=\textbf{8 dB}$',..........
    '$\textbf{Monte-Carlo: }\bar{\gamma}_{inr}^1=\textbf{10 dB}$', '$\textbf{Closed-form (9): }\bar{\gamma}_{inr}^1=\textbf{10 dB} $',.......
    '$\textbf{Monte-Carlo: }\bar{\gamma}_{inr}^1=\textbf{8 dB}, \hspace{1mm} \bar{\gamma}_{inr}^2=\textbf{7 dB} $',........
    '$\textbf{Closed-form (10): }\bar{\gamma}_{inr}^1=\textbf{8 dB}, \hspace{1mm} \bar{\gamma}_{inr}^2=\textbf{7 dB} $');  
set(h,'interpreter','latex','FontSize',12); 




end

 