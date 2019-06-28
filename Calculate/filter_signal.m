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
%    - phi: Coefficients of the serially correlated errors, in a [c x c x m] matrix, where
%       c is the number of channels, and m is the model order. If there is only one
%       channel, this becomes a [m x 1] vector
%
%   Outputs:
%    - sig_filt: Filtered signal obtained by convolving the coefficients of phi with the
%       given signal. size is [(n-p) x c], where p is the number of lags
%   

numSamples=size(sig,1);
numChannels=size(sig,2);
bool_isUnivariate=false;

if size(phi,3) == 1 && numChannels == 1
    order=size(phi,1);
    bool_isUnivariate=true;   
    phi=phi';
else
    order=size(phi,3);
end

if bool_isUnivariate
    f=[1,-phi];
%     f=[1,phi];

    tmp_sig=filter(f,1,sig);
    sig_filt=tmp_sig(1+order:end);
else
    sig_filt=zeros(numSamples-order,numChannels);

    for i=1+order:numSamples
        sig_filt(i-order,:)=sig(i,:);
        
        for j=1:order
            sig_filt(i-order,:)=sig_filt(i-order,:)-sig(i-j,:)*phi(:,:,j);
%             sig_filt(i-order,:)=sig_filt(i-order,:)+sig(i-j,:)*phi(:,:,j);
        end
    end
end

end