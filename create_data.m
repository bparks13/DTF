function [X,a]=create_data(N,m,stdZ,a)
%% [X,a]=create_data(N,m,stdZ,a)
%
%  Create a vector of values that is based on a known AR model. Used for testing accuracy
%  of AR coefficient estimation. 
%  If no inputs are given, uses univariate N = 100000, m = 2, stdZ = 2, a = [0.5 0.2]'.
%
%   Inputs:
%    - N: Number of samples. If more than one value is given, second value is assumed to
%       be the number of channels to simulate [d] (multivariate case)
%    - m: AR model order
%    - stdZ: Standard deviation of the random Gaussian error injected into the signal
%    - a: Known AR coefficients for testing purposes. Size is [m x 1] for univariate, and
%       [d x d x m] for multivariate. If a is not given, or is empty, randomly simulates the
%       AR coefficients and returns them
%
%   Outputs:
%    - X: Vector of signals created by the AR model. Size is N, either univariate or
%       multivariate
%    - a: AR coefficients that have been simulated, if no coefficients were initially
%       given
%
%  See also: mvar, estimate_ar_coefficients, dtf, estimate_residuals
%

isUnivariate=true;

if nargin == 0
    N = 100000;   % Number of points
    m = 2;      % Desired order of the AR model
    stdZ = 2;   % Desired std of the random gaussian variable
    a = [0.5 0.2]';
elseif nargin < 4 || isempty(a)
    if length(N) > 1
        numSeries=N(2);
        N=N(1);
        isUnivariate=false;
    end

    if isUnivariate
        a=zeros(m,1);
        
        for i=1:m
            a(i)=(rand(1) - 0.5) / i;
        end
    else
        a = zeros(numSeries,numSeries,m);
        I = eye(numSeries);

        for i = 1:m
            randAR = rand(numSeries);
            a(:,:,i) = I.*((randAR) - 0.5) / i + (1 - I).*((randAR - 0.5) * 0.05);
        end
    end
end

if isUnivariate
    X = zeros(N+m,1);

    for i = m+1:N+m
        X(i) = sum(X(i-1:-1:i-m).*a) + stdZ.*randn(1);
    end
    
    X=X(m+1:end);
else
    X = zeros(N + m,numSeries);           

    for i = m+1:N+m
        for j = 1:m
            X(i,:) = X(i,:) + X(i-j,:) * a(:,:,j) + stdZ*randn(1,numSeries);
        end
    end
    
    X=X(m+1:end,:);
end

end