function PD_versus_KD_in_single_RFI_detection_same_FAR_diff_INR()
%% Description: 
% This function plots the probability of RFI detection exhibited by the
% investigated power detector (PD) and kurtosis detector (KD). 
% Note that the detection performance comparison between KD and PD is made
% with respect to the desired false alarm rate (FAR) of 0.1 and different
% values of the interference-to-noise ratios (INRs). 

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Matlab code:

m_1=2;
% The SOI fading severity parameter set as per [1, Table 1].
m_2=2;
% The fading severity parameter of the RFI set as per [1, Table 1].
m=[m_1, m_2]; 
%  A vector of the fading severity parameters of the SOI and RFI 
av_snr=-5; 
% The average signal-to-noise ratio (SNR) in dB (for the assessment of PD)
av_snr_KD=[-5, -3]; 
% The average SNRs in dB (for the assessments of KD)
av_inr=0:1:10; 
% The average interference-to-noise ratios (INRs) for the RFI in dB
N=100000; 
% The number of intercepted samples (per a realization) 
num_iter=10000; 
% The number of Monte-Carlo iterations set as per [1, Table 1].
lambda_PD=1.615; 
% The decision threshold of PD rendering the desired FAR of 0.1 evaluated 
% at an average SNR of -5 dB (Q=1: single interferer scenario)  
PD_simulation=zeros(length(av_inr),length(av_snr)); 
% Initialization of PD's Monte-Carlo simulation results to be evaluated for
% different average INRs and an average SNR of -5 dB  
KD_simulation=zeros(length(av_inr),length(av_snr_KD));
% Initialization of KD's Monte-Carlo simulation results to be evaluated for
% different average INRs and average SNRs (i.e., av_snr_KD=[-5, -3] dB)  
PD_closed_form_det=zeros(length(av_inr),length(av_snr)); 
% Initialization of the PD's closed-form detection expression [1, eq. (9)] 
% to be evaluated for different average INRs and average SNR of -5 dB

for jj=1:length(av_inr)
    for ii=1:length(av_snr_KD)
        if ii==1
            PD_simulation(jj,ii)=MC_PD_wrt_FAR(m, av_snr(ii), av_inr(jj), N, num_iter); 
            % Monte-Carlo simulation for PD 
            PD_closed_form_det(jj,ii)=PD_analytical_prob_RFI_detection(lambda_PD, m, av_snr(ii), av_inr(jj));  
            % Closed-form expression [1, eq. (9)] for PD  
        end
    
        KD_simulation(jj,ii)=MC_KD_wrt_FAR(m, av_snr_KD(ii), av_inr(jj), N, num_iter); 
        % Monte-Carlo simulation for the KD. 
        % Note that the scenarios SNR=[-5, -3] dB are simulated for KD.   
    end 
end

%% Plots for PD and KD 

plot(av_inr,PD_simulation(:,1), '-ro', 'LineWidth',1); 
xlabel('$\bar{\gamma}_{inr}^1\hspace{1mm} [\textnormal{dB}]$', 'interpreter', 'latex','FontSize',12);
ylabel('\textbf{Probability of RFI detection}', 'interpreter', 'latex','FontSize',12); 
xlim([0,10]);
ylim([0,1]); 
grid on; 
hold on; 
plot(av_inr,PD_closed_form_det(:,1), '--gp', 'LineWidth',1); 
plot(av_inr,KD_simulation(:,1), '-k*', 'LineWidth',1);
plot(av_inr,KD_simulation(:,2), '-.bs', 'LineWidth',1);

h=legend('$\textbf{PD: }P_f=\textbf{0.1},\hspace{1mm} \bar{\gamma}_{snr}=-5\hspace{0.5mm}\textbf{dB}$',.........
    '$\textbf{Closed-form (9): }P_f=\textbf{0.1},\hspace{1mm} \bar{\gamma}_{snr}=-5\hspace{0.5mm}\textbf{dB}$',.........
    '$\textbf{KD [1], [15]: }P_f=\textbf{0.1}, \hspace{1mm}\bar{\gamma}_{snr}=-5\hspace{0.5mm}\textbf{dB}$',........
    '$\textbf{KD [1], [15]: }P_f=\textbf{0.1}, \hspace{1mm}\bar{\gamma}_{snr}=-3\hspace{0.5mm}\textbf{dB}$');  
set(h,'interpreter','latex','FontSize',12); 

end

       