function P_f = PD_MC_prob_false_alarm(lambda, m_1, av_snr, N, num_iter)
%% Description: 
% This function is the Monte-Carlo simulation of the probability of false alarm
% exhibited by the investigated power detector (PD). 

%% Input parameters: 
%   lambda: Decision threshold 
%   m1: The SOI fading severity parameter 
%   av_snr: The average signal-to-noise ratio (SNR) in dB
%   N: The number of intercepted samples per a realization 
%   num_iter: The number of Monte-Carlo iterations 

%% Output parameter: 
% P_f: The probability of false alarm 

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Matlab code:

av_snr=10^(0.1*av_snr);
% Conversion from the logarithmic to the linear scale. 
P=10;
% The power of the SOI set as per [1, Table 1]. 
sigma=1; 
% The square root of the noise power set as per [1, Table 1].
FAR=0; 
% Initialization for the number of false alarm instances 
for k=1:num_iter
    hs_bar=(av_snr*sigma^2)/P; 
    % The local mean received power for the SOI channel determined
    % via [1, eq. (5)]. 
    % N.B.: A BPSK modulated SOI is considered as per the setting in [1, Sec. V]. 
    h=sqrt(gamrnd(m_1,hs_bar/m_1)); 
    % Generation of the Nakagami-m distributed SOI channel gain via the gamma
    % distribution.    
    Y=0; 
    % Initialization of the mean received power.
    for n=1:N
        Y=Y+(h*sqrt(P)*(2*randi([0,1])-1)+sigma*randn)^2; 
        % Computation of the average received power as in [1, Fig. 1]. 
        % N.B.: A zero mean AWGN is simulated through ''randn'' 
        % N.B.: A BPSK modulated SOI is considered as per the simulation setting in [1, Sec. V].  
    end
    Y=Y/N; 
    % Averaging of N squared samples to approximate the expectation operation.  
    if Y>lambda
        FAR=FAR+1; 
    end 
  
end 
P_f=FAR/num_iter; 
% The probability of false alarm computed through a Monte-Carlo simulation
% which averages over ''num_iter'' channel initializations.


end

  