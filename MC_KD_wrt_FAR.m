function P_d = MC_KD_wrt_FAR(m, av_snr, av_inr, N, num_iter)
%% Description: 
% This function is the Monte-Carlo simulation of the probability of RFI detection
% exhibited by the kurtosis detection (KD) algorithm. 
% Note that the function assumes a single RFI (Q=1). However, it can
% straightforwardly be extended to the detection of multiple RFIs by
% changing the received signal model (to the one which incorporates
% multiple interferers) 

%% Input parameters: 
%   m: A vector of the fading severity parameters of the SOI and RFI
%   m_1=m(1): The SOI fading severity parameter 
%   m_2=m(2): Fading severity parameter of an RFI 
%   av_snr: The average signal-to-noise ratio (SNR) in dB
%   av_inr: The average interference-to-noise ratio (INR) in dB 
%   N: The number of intercepted samples per a realization 
%   num_iter: The number of Monte-Carlo iterations 

%% Output parameter: 
%   P_d: The probability of RFI detection. 

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Corresponding references:  
% [2] S. Misra, P. Mohammed, B. Guner, C. Ruf, J. Piepmeier, and J. Johnson, ''Microwave radiometer radio-frequency interference detection 
% algorithms: A comparative study,'' IEEE Trans. Geosci. Remote Sens., vol. 47, no. 11, pp. 3742–3754, Nov. 2009. 

% [3] R. D. De Roo, S. Misra and C. S. Ruf, ''Sensitivity of the Kurtosis Statistic as a Detector of Pulsed Sinusoidal RFI,'' in IEEE Transactions 
% on Geoscience and Remote Sensing, vol. 45, no. 7, pp. 1938-1946, Jul. 2007.

%% Matlab code: 

av_snr=10^(0.1*av_snr);
% Conversion from the logarithmic to the linear scale. 
av_inr=10^(0.1*av_inr);  
% Conversion from the logarithmic to the linear scale. 
P=10; 
% The power of the SOI set as per [1, Table 1].
P_1=10; 
% The power of the RFI set as per [1, Table 1].
sigma=1; 
% The square root of the noise power set as per [1, Table 1].
z=1.6449; 
% The normalized magnitude of the standard deviation of the kurtosis rendering 
% a false alarm rate (FAR) of 0.1, as computed using [2, eq. (5)].  
sigma_R0=sqrt(24/N);  
% The standard deviation of the kurtosis statistic (in the absence of RFI)
% valid for a very large sample setting (N>50000) rendering a normal distributed
% kurtosis [3, Sec. III]
m_1=m(1); 
% The SOI fading severity parameter set as per [1, Table 1].
m_2=m(2); 
% Fading severity parameter of the RFI set as per [1, Table 1]. 
detection=0; 
% Initialization for the number of RFI detection instances 
for k=1:num_iter 
    hs_bar=(av_snr*sigma^2)/P; 
    % The local mean received power for the SOI channel determined
    % via [1, eq. (5)]. 
    % N.B.: A BPSK modulated SOI is considered as per the setting in [1, Sec. V].  
    gs_bar=(av_inr*sigma^2)/P_1; 
    % The local mean received power for the RFI channel determined
    % via [1, eq. (5)]. 
    h=sqrt(gamrnd(m_1,hs_bar/m_1)); 
    % Generation of the Nakagami-m distributed SOI channel gain from the gamma
    % distribution 
    g=sqrt(gamrnd(m_2,gs_bar/m_2)); 
    % Generation of the Nakagami-m distributed RFI channel gain from the gamma
    % distribution 
    received_sampled=zeros(1,N); 
    % Initialization of the sampled received signal for the duration of N*T.
    % (N is the number of intercepted samples per a realization and T is the
    % symbol duration.) 
    for ii=1:N  
        alpha_i=sqrt(P)*(2*randi([0,1])-1);
        % The baseband SOI sample 
        % N.B.: A BPSK modulated SOI is considered as per the setting in [1, Sec. V]. 
        beta_i=sqrt(P_1)*randn;
        % The baseband RFI sample 
         % N.B.: A zero mean and unit variance Gaussian RFI is simulated through ''randn''
        epsilon_i=sigma*randn; 
        % AWGN of unit power (as sigma=1, as per [1, Table 1]).   
        received_sampled(ii)=h*alpha_i+g*beta_i+epsilon_i; 
        % Computation of the value of the received signal samples 
    end

    received_mean=sum(received_sampled)/N; 
    % Approximation of the expectation operation to compute the mean of 
    % received signal (averaging over N intervals)  
    rec_samples_minus_rec_mean=received_sampled-received_mean; 
    % The subtraction of the mean from the received signal samples  
    Num=sum(rec_samples_minus_rec_mean.^4)/N; 
    % Numerator of the kurtosis expression (via [2, eqs. (1) and (2)])
    % as averaged through an operation on N samples
    Den=(sum(rec_samples_minus_rec_mean.^2)/N)^2; 
    % Denominator of the kurtosis expression (via [2, eqs. (1) and (2)]) 
    % as averaged through an operation on N samples 
    kurtosis=Num/Den; 
    % The computed kurtosis (cf., [2, eq. (2)]) 

    if kurtosis<3-z*sigma_R0 || kurtosis>3+z*sigma_R0  
        % An RFI detection threshold inspired by the underneath threshold: 
        % An RFI-free detection threshold [2]: 3-z*sigma_R0<=kurtosis<=3+z*sigma_R0 
        detection=detection+1;   
    end 


end 
P_d=detection/num_iter; 
% The probability of RFI detection computed through a Monte-Carlo simulation
% which averages over ''num_iter'' channel initializations.  

end

     