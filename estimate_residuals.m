function [E,x_hat]=estimate_residuals(X,AR)
%% [E,x_hat]=estimate_residuals(X,AR)
%
%  Calculate the residuals between the AR model and the data. Used for calculating the
%  log-likelihood
%
%   Inputs:
%    - X: Vector of signals to estimate. Size is [N x 1], where N is the number of samples
%    - AR: AR coefficients for the model, size is [m x 1], where m is the model order
%   
%   Outputs:
%    - E: Residuals between the data given (X) and the model (AR), size is [N - m x 1],
%       where the first m points of the data are used to begin the model
%    - x_hat: Output of the model that estimates the signal
%
%  See also: estimate_ar_coefficients, mvar, test_model
%

m=length(AR);

x_hat=zeros(length(X)-m,1);

for i=m+1:length(X)
    x_hat(i-m)=sum(X(i-1:-1:i-m).*AR);
end

E=X(m+1:end)-x_hat;

end