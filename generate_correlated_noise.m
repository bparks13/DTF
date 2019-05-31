function noise=generate_correlated_noise(N,m,rho,stdZ)
%% noise=generate_correlated_noise(N,m,rho,stdZ)
%
%  Generates noise that is correlated within the model order defined, for a given number
%  of samples. Based on the correlated noise model defined within
%  https://www.osapublishing.org/boe/abstract.cfm?uri=boe-4-8-1366. Correlated noise model
%  is e_t = rho_1 * e_(t-1) + v_t, for an AR(1) model, where v_t is uncorrelated white
%  noise.
%
%   Inputs:
%    - N: Number of samples to generate; if two numbers are defined this will generate
%       multiple channels of correlated noise, but with no cross-correlation
%    - m: Model order, defines the AR(m) model
%    - rho: Correlated noise AR coefficients; size should match the model order and the
%       number of channels defined. size is [m x c], where m is the model order and c is
%       the number of channels
%    - stdZ: Standard deviation of the white noise
%
%   Outputs:
%    - noise: Correlated noise for the given channels
%

if length(rho) ~= m
    disp('ERROR: Number of coefficients does not match model order');
    noise=nan;
    return;
end

if size(rho,1) < size(rho,2)
    rho=rho';
end

numSamples=N(1);
numChannels=1;
bool_isUnivariate=true;

if length(N) > 1
    numChannels=N(2);
    bool_isUnivariate=false;
end

noise=zeros(numSamples+m,numChannels);

if bool_isUnivariate
    for i=m+1:numSamples+m
        noise(i,:)=sum(noise(i-1:-1:i-m).*rho) + stdZ*randn(1);
    end
else
    for i=m+1:numSamples+m
        for j=1:m
            noise(i,:)=noise(i,:) + noise(i-j,:) * rho(:,:,j);
        end
        noise(i,:)=noise(i,:) + stdZ*randn(1);
    end
end

noise=noise(m+1:end,:);

end