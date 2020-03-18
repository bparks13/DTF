function sig_filt=filter_signal(sig,phi)
%% sig_filt=filter_signal(sig,phi)
%
%  Filter a signal based on the serially correlated residuals, given by phi, where the new
%  signal is sig_filt = filter([1 -phi(1) -phi(2) ... -phi(3)],1,sig) for each individual
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
%       given signal. size is [n x c]. The first m points are calculated using the 
%   
%  See also: filter_serial_correlation
%

numSamples=size(sig,1);
numChannels=size(sig,2);
bool_isUnivariate=false;

if numChannels == 1
    if size(phi,3)==1
        order=size(phi,1);
        phi=phi';
    else
        order=size(phi,3);
        phi=squeeze(phi)';
    end
    
    bool_isUnivariate=true;   
else
    order=size(phi,3);
end

if bool_isUnivariate
    f=[1,-phi];

    tmp_sig=filter(f,1,sig);
    sig_filt=zeros(size(tmp_sig));
    
    for i=1:order
        prevSig=-phi(1:i-1)*sig(1:i-1);
        
        if isempty(prevSig)
            sig_filt(i)=sig(i)*sqrt(1-phi(i)^2);
        else
            sig_filt(i)=sig(i)*sqrt(1-phi(i)^2)+prevSig;
        end
    end
    
    sig_filt(1+order:end)=tmp_sig(1+order:end);
else
    sig_filt=zeros(numSamples,numChannels);
    I=eye(numChannels);
    
    for i=1:order
        sig_filt(i,:)=sig(i,:)*sqrt(I-phi(:,:,i).^2);
        
        for j=i-1:-1:1
            sig_filt(i,:)=sig_filt(i,:)-sig(j,:)*phi(:,:,j);
        end
    end

    for i=1+order:numSamples
        sig_filt(i,:)=sig(i,:);
        
        for j=1:order
            sig_filt(i,:)=sig_filt(i,:)+sig(i-j,:)*(-phi(:,:,j));
        end
    end
end

end