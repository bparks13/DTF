function [E,C,x_hat]=estimate_residuals(X,AR)
%% [E,C,x_hat]=estimate_residuals(X,AR)
%
%  Calculate the residuals between the AR model and the data. Used for calculating the
%  log-likelihood
%
%   Inputs:
%    - X: Vector of signals to estimate. Size is [N x d], where N is the number of
%       samples, and d is the number of series
%    - AR: AR coefficients. Size is [m x 1] for univariate, and [d x d x m] for
%       multivariate, where m is the model order 
%   
%   Outputs:
%    - E: Residuals between the data given (X) and the model (AR), size is [N - m x d],
%       where the first m points of the data are used to begin the model
%    - C: Covariance matrix defined by the residuals, size is [d x d]
%    - x_hat: Output of the model that estimates the signal
%
%  See also: estimate_ar_coefficients, mvar, test_model
%

if size(X,2) > 1
    isUnivariate=false; 
else
    isUnivariate=true;
end

if isUnivariate
    m=length(AR);

    x_hat=zeros(length(X)-m,1);

    for i=m+1:length(X)
        x_hat(i-m)=sum(X(i-1:-1:i-m).*AR);
    end

    E=X(m+1:end)-x_hat;
else
    m=size(AR,3);
    
    x_hat=zeros(size(X,1)-m,size(X,2));
    
    for i=m+1:length(X)
        for j=1:m
            x_hat(i-m,:)=x_hat(i-m,:)+X(i-j,:)*(AR(:,:,j));
        end
    end
    
    E=X(m+1:end,:)-x_hat;
end

C=(E' * E) / length(E);

end