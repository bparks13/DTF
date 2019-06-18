function [X,noise,a,rho]=create_data_correlated_noise(N,m,a,rho,stdZ)
%% [X,noise,a,rho]=create_data_correlated_noise(N,m,a,rho,stdZ)
%
%  Create a vector of values that is based on a known AR model. Used for testing accuracy
%  of AR coefficient estimation. Uses noise that is autocorrelated
%
%   Inputs:
%    - N: Number of samples. If more than one value is given, second value is assumed to
%       be the number of channels to simulate [d] (multivariate case)
%    - m: Model orders for the AR coefficients and the correlated noise coefficients. If
%       there is a single value, it is assumed to be the same for both; otherwise, m(1) is
%       defined to be AR model order, and m(2) is defined to be the noise model order
%    - a: Known AR coefficients for testing purposes. Size is [m x 1] for univariate, and
%       [d x d x m] for multivariate. If a is not given, or is empty, randomly simulates the
%       AR coefficients and returns them
%    - rho: Known AR coefficients for correlated noise. Size is the same as a. If empty,
%       simulates the coefficients and returns them.
%    - stdZ: Standard deviation of the random Gaussian error injected into the signal
%
%   Outputs:
%    - X: Vector of signals created by the AR model. Size is N, either univariate or
%       multivariate
%    - a: AR coefficients that have been simulated, if no coefficients were initially
%       given
%    - Rho: AR coefficients for the correlated noise that have been simulated, if no
%       coefficients were initially given
%
%  See also: mvar, estimate_ar_coefficients, dtf, estimate_residuals,
%   comparing_ar_algorithms_correlated_noise
%

rng('shuffle');

numSamples=N(1);
numChannels=1;
bool_isUnivariate=true;

if length(N) > 1
    numChannels=N(2);
    bool_isUnivariate=false;
end

if length(m) == 1
    m=[m,m];
end

model_ar=m(1);
model_rho=m(2);

if isempty(a)
    if bool_isUnivariate
        a=zeros(model_ar,1);
        
        for i=1:model_ar
            a(i)=(rand(1) - 0.5) / i;
        end
    else
        a = zeros(numChannels,numChannels,model_ar);
        I = eye(numChannels);

        for i = 1:model_ar
            randAR = rand(numChannels);
            a(:,:,i) = I.*((randAR) - 0.5) / i + (1 - I).*((randAR - 0.5) * 0.05);
        end
    end
end

if isempty(rho) && model_rho ~= 0
    if bool_isUnivariate
        rho=zeros(model_rho,1);
        
        for i=1:model_rho
            rho(i)=(rand(1) - 0.5) / i;
        end
    else
        rho = zeros(numChannels,numChannels,model_rho);
        I = eye(numChannels);

        for i = 1:model_rho
            randAR = rand(numChannels);
            rho(:,:,i) = I.*((randAR) - 0.5) / i + (1 - I).*((randAR - 0.5) * 0.05);
        end
    end
elseif rho==0 || model_rho == 0
    rho=zeros(size(a));
end

noise=generate_correlated_noise(N,model_rho,rho,stdZ);

if bool_isUnivariate
    X = zeros(numSamples + model_ar,1);

    for i = model_ar+1:numSamples+model_ar
        X(i) = sum(X(i-1:-1:i-model_ar).*a) + noise(i-model_ar);
    end
    
    X=X(model_ar+1:end);
else
    X = zeros(numSamples + model_ar,numChannels);           

    for i = model_ar+1:numSamples+model_ar
        for j = 1:model_ar
            X(i,:) = X(i,:) + X(i-j,:) * a(:,:,j);
        end
        X(i,:) = X(i,:) + noise(i-model_ar,:);
    end
    
    X=X(model_ar+1:end,:);
end

if nargout == 2
    clear a rho
end

end