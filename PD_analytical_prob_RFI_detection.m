function P_d = PD_analytical_prob_RFI_detection(lambda, m, av_snr, av_inr )
%% Description: 
% This function computes the probability of RFI detection (via the analytical expressions)
% exhibited by the investigated power detector (PD). Note that the function can 
% be used to compute the probability of an RFI detection (Q=1) and the probability of
% the detection of two RFIs (Q=2). 

%% Input parameters:  
%   lambda: Decision threshold 
%   m: A vector of the fading severity parameters of the SOI and RFI(s)
%   m_1=m(1): The SOI fading severity parameter 
%   m_2=m(2): Fading severity parameter of the first RFI 
%   m_3=m(3): Fading severity parameter of the second RFI
%   av_snr: The average signal-to-noise ratio (SNR) in dB
%   av_inr: A vector of average interferene-to-noise ratios (INRs) 
%   av_inr_1=av_inr(1): The average INR of the first RFI in dB 
%   av_inr_2=av_inr(2): The average INR of the second RFI in dB

%% Output parameter:  
%   P_d: The probability of RFI detection (Q=1 and Q=2). 

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
Q=length(av_inr); 
% The number of impinging RFIs.
if Q>2
    error('The function is written for one RFI or two impinging RFIs');
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
        eta_1=(sigma^2*av_snr)/m_1; 
        % Definition of eta_1 as defined beneath [1, eq. (7)].   
        eta_2=(sigma^2*av_inr)/m_2; 
         % Definition of eta_2 as defined beneath [1, eq. (7)]. 
        PD_sum=0; 
        % Initialization.
        for ii=1:2
            for k=1:m(ii)
                second_argum=m(ii)*(lambda/sigma^2-1)/(U(1-ii)*av_snr+U(ii-2)*av_inr);  
                % The second argument in the (upper) incomplete gamma function of [1, eq. (9)].  
                PD_sum=PD_sum+Chi2(ii, k, m_1, m_2, eta_1, eta_2)*igamma(k,second_argum)/factorial(k-1); 
                % Accumulation of the double summation summands to evaluate [1, eq. (9)].
            end
        end

            P_d=PD_sum;
    end
     
    if Q==2
        av_inr_1=10^(0.1*av_inr(1));
        % The average INR of the first RFI.
        % Conversion from the logarithmic to the linear scale.
        av_inr_2=10^(0.1*av_inr(2));
        % The average INR of the second RFI.
        % Conversion from the logarithmic to the linear scale. 
        av_inr=[av_inr_1, av_inr_2]; 
        % The average INR of the first and second RFIs 
        m_1=m(1); 
        % The SOI fading severity parameter set as per [1, Table 1]. 
        m_2=m(2); 
        % Fading severity parameter of the first RFI set as per [1, Table 1]. 
        m_3=m(3); 
        % Fading severity parameter of the second RFI set as per [1, Table 1]. 
        eta_1=(sigma^2*av_snr)/m_1; 
        % Definition of eta_1 as defined beneath [1, eq. (7)]. 
        eta_2=(sigma^2*av_inr_1)/m_2; 
        % Definition of eta_2 as defined beneath [1, eq. (7)]. 
        eta_3=(sigma^2*av_inr_2)/m_3; 
        % Definition of eta_3 as defined beneath [1, eq. (7)]. 
        PD_sum=0; 
        % Initialization
        for ii=1:3
            for k=1:m(ii)
                if ii==1
                    second_argum=m(ii)*(lambda/sigma^2-1)/(U(1-ii)*av_snr);
                     % The second argument in the (upper) incomplete gamma function of [1, eq. (10)]. 
                else
                    second_argum=m(ii)*(lambda/sigma^2-1)/(U(1-ii)*av_snr+U(ii-2)*av_inr(ii-1));
                    % The second argument in the (upper) incomplete gamma function of [1, eq. (10)].  
                end 
                PD_sum=PD_sum+Chi3(ii, k, m_1, m_2, m_3, eta_1, eta_2, eta_3)*igamma(k,second_argum)/factorial(k-1); 
                % Accumulation of the double summation summands to evaluate [1, eq. (10)].
            end
        end

            P_d=PD_sum;
     end
end



end



 