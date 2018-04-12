function y=Chi2(ii, k, m_1, m_2, eta_1, eta_2)
%% Description: 
% This function computes the weight ''Chi2'' defined in [2, eq. (A-5)].

%% Input parameters:  
%   ii: Index counting toward the total number of random variables (RVs)
%   (in this case, there are only two RVs, i.e., the SOI and the RFI) 
%   k: Index counting toward the fading severity parameter of the ith RV 
%   m_1: The SOI fading severity parameter 
%   m_2: Fading severity parameter of the RFI 
%   eta_1: Shape parameter (with respect to the SOI) as defined beneath [1, eq. (7)].   
%   eta_2: Shape parameter (with respect to the RFI) as defined beneath [1, eq. (7)].


%% Output parameter:  
%   y: The computed weight as per [2, eq. (A-5)]  

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Corresponding reference:
% [2] G. K. Karagiannidis, N. C. Sagias, and T. A. Tsiftsis, ''Closed-form statistics for the sum of squared Nakagami-m variates and its applications,'' 
% IEEE Trans. Commun., vol. 54, no. 8, pp. 1353–1359, Aug 2006

%% Matlab code:

m_vec=[m_1, m_2];
% A vector of the fading severity parameters of the SOI and RFI.
eta_vec=[eta_1, eta_2];
% A vector of the shape parameters with respect to the SOI and RFI. 
R2=sum(m_vec); 
% The summation of the fading severity parameters of the SOI and RFI.
Num=(eta_vec(ii)^k)*factorial(sum(m_vec)-k-1)*(1/eta_vec(ii)-1/eta_vec(1+U(1-ii)))^(k-sum(m_vec));
% The numerator expression of [2, eq. (A-5)] (without the multiplying constant)   
Den=(eta_1^m_1)*(eta_2^m_2)*factorial(m_vec(1+U(1-ii))-1)*factorial(m_vec(ii)-k);
% The Denominator expression of [2, eq. (A-5)] 
y=(-1)^(R2-m_vec(ii))*Num/Den; 
% The overall expression of [2, eq. (A-5)]   


end

