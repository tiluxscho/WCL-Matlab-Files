function weight=Chi3(ii, k, m_1, m_2, m_3, eta_1, eta_2, eta_3)
%% Description: 
% This function computes the weight ''Chi3'' defined in [2, eq. (A-8)].

%% Input parameters:  
%   ii: Index counting toward the total number of random variables (RVs)
%   (in this case, there are only three RVs, i.e., the SOI, the first RFI, 
%   and the second RFI) 
%   m_1: The SOI fading severity parameter 
%   m_2: Fading severity parameter of the first RFI 
%   m_3: Fading severity parameter of the second RFI 
%   eta_1: Shape parameter (with respect to the SOI) as defined beneath [1, eq. (7)].   
%   eta_2: Shape parameter (with respect to the first RFI) as defined beneath [1, eq. (7)].
%   eta_3: Shape parameter (with respect to the second RFI) as defined beneath [1, eq. (7)]. 

%% Output parameter:  
%   weight: The computed weight as per [2, eq. (A-8)]  

%% Author: Tilahun M. Getu 

%% Corresponding paper: 
% [1] Tilahun M. Getu, W. Ajib, and Rene Jr. Landry, ''Power-based broadband RF interference detector for wireless communication systems,'' 
% IEEE Wireless Commun. Lett., submitted, Apr. 2018.
% Date: Apr. 2018

%% Corresponding reference:
% [2] G. K. Karagiannidis, N. C. Sagias, and T. A. Tsiftsis, ''Closed-form statistics for the sum of squared Nakagami-m variates and its applications,'' 
% IEEE Trans. Commun., vol. 54, no. 8, pp. 1353–1359, Aug 2006

%% Matlab code:

m_vec=[m_1, m_2, m_3];
% A vector of the fading severity parameters of the SOI, the first RFI, and the second RFI.
eta_vec=[eta_1, eta_2, eta_3];
% A vector of the shape parameters with respect to the SOI, the first RFI, and the second RFI. 
R3=sum(m_vec); 
% The summation of the fading severity parameters of the SOI, the first RFI, and the second RFI. 
y=0; 
% Initialization of the weight summand 
for l_1=k:m_vec(ii)    
    Num_1=(eta_vec(ii)^k)*factorial(m_vec(ii)+m_vec(1+U(1-ii))-l_1-1)*(1/eta_vec(ii)-1/eta_vec(1+U(1-ii)))^(l_1-m_vec(ii)-m_vec(1+U(1-ii))); 
    Den_1=(eta_1^m_1)*(eta_2^m_2)*(eta_3^m_3)*factorial(m_vec(1+U(1-ii))-1)*factorial(m_vec(ii)-l_1); 
    % Numerator and denominator for the first expression of [2, eq. (A-8)]
    Num_2=factorial(l_1+m_vec(2+U(2-ii))-k-1)*(1/eta_vec(ii)-1/eta_vec(2+U(2-ii)))^(k-l_1-m_vec(2+U(2-ii)));  
    Den_2=factorial(m_vec(2+U(2-ii))-1)*factorial(l_1-k);  
    % Numerator and denominator for the second expression of [2, eq. (A-8)]
    y=y+(-1)^(R3-m_vec(ii))*(Num_1/Den_1)*(Num_2/Den_2); 
    % The overall expression of [2, eq. (A-8)]   
end

weight=y; 
% The computed weight as per [2, eq. (A-8)]  

end

