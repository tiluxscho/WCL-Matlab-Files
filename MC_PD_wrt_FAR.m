function P_d=MC_PD_wrt_FAR(m, av_snr, av_inr, N, num_iter)
%% Description: 
% This function calculates the Monte-Carlo simulation of the probability of 
% RFI detection exhibited by the investigated power detector (PD).  
% Note that the simulated detection probability is with respect to the
% desired false alarm rate (FAR)

%% Input parameters: 
%   m: A vector of the fading severity parameters of the SOI and RFI(s)
%   m_1=m(1): The SOI fading severity parameter 
%   m_2=m(2): Fading severity parameter of the first RFI 
%   m_3=m(3): Fading severity parameter of the second RFI
%   av_snr: The average signal-to-noise ratio (SNR) in dB
%   av_inr: A vector of average interference-to-noise ratios (INRs) 
%   av_inr_1=av_inr(1): The average INR of the first RFI in dB 
%   av_inr_2=av_inr(2): The average INR of the second RFI in dB
%   N: The number of intercepted samples per a realization 
%   num_iter: The number of Monte-Carlo iterations 

%% Output parameter: 
%   P_d: The probability of RFI detection (Q=1 and Q=2). 

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
P_1=10; 
% The power of the first RFI set as per [1, Table 1].
P_2=10; 
% The power of the second RFI set as per [1, Table 1].
sigma=1; 
% The square root of the noise power set as per [1, Table 1].
Q=length(av_inr); 
% The number of impinging RFIs.
lambda=1.615; 
% The decision threshold rendering a false alarm rate (FAR) of 0.1 
% Note that this value is evaluated at an average SNR of -5 dB using [1, eq. (12)].  
detection=0; 
% Initialization for the number of RFI detection instances 
if Q>2
    error('The function is written for one or two impinging RFIs');
    return;
else
    if Q==1
        av_inr=10^(0.1*av_inr);
        % The average INR of the impinging RFI.  
        % Conversion from the logarithmic to the linear scale. 
        % N.B.: av_inr=av_inr(1), as we got only one RFI 
        m_1=m(1); 
        % The SOI fading severity parameter set as per [1, Table 1].
        m_2=m(2); 
        % Fading severity parameter of the first RFI set as per [1, Table 1].
        for k=1:num_iter
            hs_bar=(av_snr*sigma^2)/P; 
            % The local mean received power for the SOI channel determined
            % via [1, eq. (5)]. 
            % N.B.: A BPSK modulated SOI is considered as per the setting in [1, Sec. V]. 
            gs_bar=(av_inr*sigma^2)/P_1; 
            % The local mean received power for the first RFI channel determined
            % via [1, eq. (5)].  
            h=sqrt(gamrnd(m_1,hs_bar/m_1)); 
            % Generation of the Nakagami-m distributed SOI channel gain via the gamma
            % distribution.
            g=sqrt(gamrnd(m_2,gs_bar/m_2)); 
            % Generation of the Nakagami-m distributed first RFI channel gain via the gamma
            % distribution.
            Y=0; 
            % Initialization of the mean received power.
            for n=1:N
                Y=Y+(h*sqrt(P)*(2*randi([0,1])-1)+g*sqrt(P_1)*randn+sigma*randn)^2; 
                % Computation of the average received power as in [1, Fig. 1]. 
                % N.B.: A zero mean and unit variance Gaussian RFI is simulated through ''randn'' 
                % N.B.: A BPSK modulated SOI is considered as per the simulation setting in [1, Sec. V].
            end
            Y=Y/N; 
            % Averaging of N squared samples to approximate the expectation operation. 
            if Y>lambda
                detection=detection+1; 
            end 
  
        end 
        P_d=detection/num_iter; 
        % The probability of RFI detection computed through a Monte-Carlo simulation
        % which averages over ''num_iter'' channel initializations. 
    end
    if Q==2
        av_inr_1=10^(0.1*av_inr(1));
        % The average INR of the first RFI.
        % Conversion from the logarithmic to the linear scale.
        av_inr_2=10^(0.1*av_inr(2));
        % The average INR of the second RFI.
        % Conversion from the logarithmic to the linear scale.
        m_1=m(1); 
         % The SOI fading severity parameter set as per [1, Table 1].
        m_2=m(2); 
        % Fading severity parameter of the first RFI set as per [1, Table 1].
        m_3=m(3); 
        % Fading severity parameter of the second RFI set as per [1, Table 1].
        for k=1:num_iter
            hs_bar=(av_snr*sigma^2)/P; 
            % The local mean received power for the SOI channel determined
            % via [1, eq. (5)].  
            % N.B.: A BPSK modulated SOI is considered as per the setting in [1, Sec. V].   
            gs_bar_1=(av_inr_1*sigma^2)/P_1; 
            % The local mean received power for the first RFI channel determined
            % via [1, eq. (5)]. 
            gs_bar_2=(av_inr_2*sigma^2)/P_2; 
            % The local mean received power for the second RFI channel determined
            % via [1, eq. (5)]. 
            h=sqrt(gamrnd(m_1,hs_bar/m_1)); 
            % Generation of the Nakagami-m distributed SOI channel gain via the gamma
            % distribution.  
            g_1=sqrt(gamrnd(m_2,gs_bar_1/m_2)); 
            % Generation of the Nakagami-m distributed first RFI channel gain via the gamma
            % distribution.  
            g_2=sqrt(gamrnd(m_3,gs_bar_2/m_3)); 
            % Generation of the Nakagami-m distributed second RFI channel gain via the gamma
            % distribution.       
            Y=0; 
            % Initialization of the mean received power.
            for n=1:N 
                Y=Y+(h*sqrt(P)*(2*randi([0,1])-1)+g_1*sqrt(P_1)*randn+g_2*sqrt(P_2)*randn+sigma*randn)^2; 
                % Computation of the average received power as in [1, Fig. 1]. 
                % N.B.: A zero mean and unit variance Gaussian RFIs are simulated through ''randn''
                % N.B.: A BPSK modulated SOI is considered as per the simulation setting in [1, Sec. V].
            end
            Y=Y/N; 
        %    Averaging of N squared samples to approximate the expectation operation. 
            if Y>lambda
                detection=detection+1; 
            end 
  
        end 
        P_d=detection/num_iter; 
        % The probability of RFI detection computed through a Monte-Carlo simulation
        % which averages over ''num_iter'' channel initializations. 
        
        
    end
    

end


end

    