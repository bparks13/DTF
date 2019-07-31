function [AR,C]=estimate_ar_coefficients(X,m)
%% [AR,C]=estimate_ar_coefficients(X,m)
%
%  Given a signal and the model order desired, returns the AR coefficients (AR) and
%  the variance of the error term from the estimated parameters (C)
%
%   Inputs:
%    - X: Vector of signals to estimate. Size is [N x d], where N is the number of samples
%       and d is number of series. Can either be univariate or multivariate.
%    - m: AR model order to estimate
%
%   Outputs:
%    - AR: Estimated AR coefficients, size is [m x 1] for univariate signals, and 
%       [d x d x m] for multivariate signals
%    - C: Estimated variance of the error term. Size is scalar for univariate signals, and
%       [d x d] for multivariate signals
%
%  See also: mvar, dtf, estimate_residuals
%

if size(X,2) > 1
    isUnivariate=false;
    numSeries=size(X,2);   
else
    isUnivariate=true;
end

if isUnivariate
    [rxx,lags] = xcov(X,m,'unbiased');
    r_n = rxx(lags>0);
    idx0 = find(lags==0);
    R_n = zeros(m,m);

    for i=1:m
        tmpidx = idx0 -(i-1);
        R_n(:,i) = rxx(tmpidx:tmpidx+m-1);
    end

    AR = R_n\r_n; % Predictor for AR coefficients
    C = flip(rxx(lags<=0))'*[1; -AR];
else
    R_all=zeros(m*2+1,numSeries^2);

    for i=1:numSeries^2
        a=ceil(i/numSeries);
        b=mod(i-1,numSeries)+1;
%         fprintf('i = %d, a = %d, b = %d\n',i,a,b);
        [R_all(:,i),lags] = xcov(X(:,a),X(:,b),m,'unbiased');
    end

    r_n = zeros(m*numSeries,numSeries);
    
    for j=1:numSeries
        for i=1:m*numSeries
%             a=(j-1)*numSeries + mod(i-1,numSeries)+1;
%             b=-ceil(i/numSeries);
            a=j+mod(i-1,numSeries)*numSeries;
            b=ceil(i/numSeries);
%             fprintf('(%d,%d) a = %d, b = %d\n',i,j,a,b);
            r_n(i,j)=R_all(lags==b,a);
        end
    end

    R_n = zeros(m*numSeries);

    for i=1:m*numSeries
        for j=1:m*numSeries
            a=mod((i-1)*numSeries+mod(j-1,numSeries),numSeries^2)+1;
            b=-floor((j-1)/numSeries)+floor((i-1)/numSeries);
%             fprintf('(%d,%d) a = %d, b = %d\n',i,j,a,b);
            R_n(i,j)=R_all(lags==b,a);
        end
    end
    
    theta = (R_n \ r_n);
    AR=zeros(numSeries,numSeries,m);
    
    for i=1:m
        ind=(i-1)*numSeries+1:i*numSeries;
        AR(:,:,i) = theta(ind,:);
    end
     
    if nargout > 1 
        C=reshape(R_all(lags==0,:),numSeries,numSeries);

        for i=1:m
            C=C+reshape(R_all(lags==-i,:),numSeries,numSeries)' * -AR(:,:,i); 
        end
    end
end

end