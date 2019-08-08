function phi_hat=calculate_serial_coefficients(E,order)
%% phi_hat=calculate_serial_coefficients(E,order)
%
%   Inputs:
%    - E: Residuals from the AR model. Matrix with size [n x c], where n is the number of
%       samples and c is the number of channels
%    - order: Model order (number of lags) to calculate the coefficients for
%
%   Outputs:
%    - phi_hat: Serial correlated coefficients of the residuals. Matrix of size [c x c x o]
%       where o is the model order, and c is the number of channels
%
%  See also: filter_serial_correlation, filter_signal
%

numChannels=size(E,2);

phi_hat=zeros(numChannels,numChannels,order);

for i=1:numChannels
    phi_hat(i,i,:)=internal_calculate_acf(E(:,i),order);
end

%%
    function acf=internal_calculate_acf(x,maxLag)
        N=length(x);
        nFFT = 2^(nextpow2(N)+1);
        F = fft(x-mean(x),nFFT);
        F = F.*conj(F);
        acf = ifft(F);
        acf = acf(1:(maxLag+1)); % Retain non-negative lags
        acf = acf./acf(1); % Normalize
        acf = real(acf);
        acf=acf(2:end);
    end
end