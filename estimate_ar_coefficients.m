function [AR,C]=estimate_ar_coefficients(X,m)
%% [AR,SigmaZ_hat]=estimate_ar_coefficients(X,m)
%
%  Given a signal and the model order desired, returns the AR coefficients (AR) and
%  the variance of the error term from the estimated parameters (C)
%
%   Inputs:
%    - X: Vector of signals to estimate. Size is [N x 1], where N is the number of samples
%    - m: AR model order to estimate
%
%   Outputs:
%    - AR: Estimated AR coefficients, size is [m x 1]
%    - C: Estimated variance of the error term 
%
%  See also: mvar, dtf, estimate_residuals
%

[rxx,lags] = xcov(X,m,'unbiased');
r_n = rxx(lags>0);
idx0 = find(lags==0);
R_n = zeros(m,m);

for i=1:m
    tmpidx = idx0 -(i-1);
    R_n(:,i) = rxx(tmpidx:tmpidx+m-1);
end

AR = R_n\r_n; % Predictor for AR coefficients
C = rxx(lags>=0)'*[1; -AR];
% Sigma_hat = rxx(idx0);

end