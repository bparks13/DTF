function [Theta_hat,SigmaZ_hat,sigma_hat]=estimate_ar_coefficients(X,m)
%% [Theta_hat,SigmaZ_hat,sigma_hat]=estimate_ar_coefficients(X,m)
%
%  Given a signal and the model order desired, returns the AR coefficients (Theta_hat) and
%  the variance of the error term from the estimated parameters (SigmaZ_hat)
%
%   Inputs:
%    - X: Vector of signals to estimate. Size is [N x 1], where N is the number of samples
%    - m: AR model order to estimate
%
%   Outputs:
%    - Theta_hat: Estimated AR coefficients, size is [m x 1]
%    - SigmaZ_hat: Estimated variance of the error term 
%    - sigma_hat: Estimated variance of the signal (sigma_squared)
%

[rxx,lags] = xcov(X,m,'unbiased');
r_n = rxx(lags>0);
idx0 = find(lags==0);
R_n = zeros(m,m);

for i=1:m
    tmpidx = idx0 -(i-1);
    R_n(:,i) = rxx(tmpidx:tmpidx+m-1);
end

Theta_hat = R_n\r_n; % Predictor for AR coefficients
SigmaZ_hat = rxx(lags>=0)'*[1; -Theta_hat];
sigma_hat = rxx(idx0);

end