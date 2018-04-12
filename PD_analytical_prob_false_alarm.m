function P_f = PD_analytical_prob_false_alarm(lambda, m_1, av_snr )
%% Description: 
% This function computes the probability of false alarm (analytical expression)
% exhibited by the investigated power detector (PD). 
%% Input parameters: 
%   lambda: Decision threshold 
%   m1: The SOI fading severity parameter 
%   av_snr: The average signal-to-noise ratio (SNR) in dB

%% Output parameter: 
% P_f: The probability of false alarm 

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Matlab code: 

sigma=1; 
% The square root of the noise power set as per [1, Table 1].   
av_snr=10^(0.1*av_snr);
% Conversion from the logarithmic to the linear scale. 
P_f=igamma(m_1, m_1*(lambda/sigma^2-1)/av_snr)/factorial(m_1-1);
% The derived analytical false alarm expression [1, eq. (12)].  

end
