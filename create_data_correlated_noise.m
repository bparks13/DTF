function [X,noise,a,rho]=create_data_correlated_noise(N,m,a,stdZ,rho,stdRho)
%% [X,a,rho]=create_data_correlated_noise(N,m,a,stdZ,rho,stdRho)
%
%  Create a vector of values that is based on a known AR model. Used for testing accuracy
%  of AR coefficient estimation. Uses noise that is autocorrelated
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

numSamples=N(1);
numChannels=1;
bool_isUnivariate=true;

if length(N) > 1
    numChannels=N(2);
    bool_isUnivariate=false;
end

if isempty(a)
    if bool_isUnivariate
        a=zeros(m,1);
        
        for i=1:m
            a(i)=(rand(1) - 0.5) / i;
        end
    else
        a = zeros(numChannels,numChannels,m);
        I = eye(numChannels);

        for i = 1:m
            randAR = rand(numChannels);
            a(:,:,i) = I.*((randAR) - 0.5) / i + (1 - I).*((randAR - 0.5) * 0.05);
        end
    end
end

if isempty(rho)
    if bool_isUnivariate
        rho=zeros(m,1);
        
        for i=1:m
            rho(i)=(rand(1) - 0.5) / i;
        end
    else
        rho = zeros(numChannels,numChannels,m);
        I = eye(numChannels);

        for i = 1:m
            randAR = rand(numChannels);
            rho(:,:,i) = I.*((randAR) - 0.5) / i + (1 - I).*((randAR - 0.5) * 0.05);
        end
    end
end

noise=generate_correlated_noise(N,m,rho,stdRho);

if bool_isUnivariate
    X = zeros(numSamples + m,1);

    for i = m+1:numSamples+m
        X(i) = sum(X(i-1:-1:i-m).*a) + noise(i-m);
    end
    
    X=X(m+1:end);
else
    X = zeros(numSamples + m,numChannels);           

    for i = m+1:numSamples+m
        for j = 1:m
            X(i,:) = X(i,:) + X(i-j,:) * a(:,:,j);
        end
        X(i,:) = X(i,:) + noise(i-m,:);
    end
    
    X=X(m+1:end,:);
end

end