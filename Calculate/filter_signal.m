function sig_filt=filter_signal(sig,phi)
%% sig_filt=filter_signal(sig,phi)
%
%  Filter a signal based on the serially correlated residuals, given by phi, where the new
%  signal is sig_filt = filter([1 phi(1) phi(2) ... phi(3)],1,sig) for each individual
%  channel (no cross correlation correction performed)
%
%   Inputs:
%    - sig: Signal contained in a matrix of size [n x c], where n is the number of samples
%       and c is the number of channels
%    - phi: Coefficients of the serially correlated errors, in a [o x c] matrix, where o
%       is the number of lags (model order) and c is the number of channels
%
%   Outputs:
%    - sig_filt: Filtered signal obtained by convolving the coefficients of phi with the
%       given signal. Same size as sig.
%   

numChannels=size(sig,2);
numCoeff=size(phi,1);

sig_filt=zeros(size(sig));

for i=1:numChannels
    f=zeros(1,numCoeff+1);
    f(1)=1;
    
    for j=1:numCoeff
        f(j+1)=phi(j,i);
    end
    
    sig_filt(:,i)=filter(f,1,sig(:,i));
end

% sig_filt=conv2(-phi,sig);

end