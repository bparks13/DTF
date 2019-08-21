function [connectivity,H]=dtf(mdl,freqRange,fs,normalize)
%% [connectivity,H]=dtf(mdl,freqRange,fs,normalize)
%
%  Using the estimated AR model, calculate the connectivity of the series in question. Can
%  return either the raw connections or the normalized connections as specified in config
%
%   Inputs:
%    - mdl: Estimated AR model struct as returned by mvar
%       AR: Autoregressive coefficients as found by estimate_ar_coefficients
%       order: Model order that is found to have the lowest information criterion
%       numSeries: Number of channels/series used 
%    - freqRange: Vector containing the specific frequencies to calculate the connectivity
%       over
%    - fs: Sampling frequency in Hz
%    - normalize: Boolean defining whether to normalize [true, default] or not [false]
%
%   Outputs:
%    - connectivity: Directed transfer function values for all combinations of series and
%       frequencies; either normalized [default] or not, depending on if config is set
%    - H: System transfer function
%
%  See also: mvar, plot_connectivity
%

if nargin < 4
    normalize=true;
end

nFreqs=length(freqRange);
numSeries=mdl.numSeries;
modelOrder=mdl.order;
I=eye(numSeries);

AR=mdl.AR;

tmp_pdc=zeros(numSeries,numSeries,nFreqs);

% Calculate PDC
for i=1:nFreqs
    tmp_pdc(:,:,i)=I-reshape(sum(bsxfun(@times,reshape(AR,numSeries^2,modelOrder),...
        exp(-(2*pi*1i*(1:modelOrder))*(freqRange(i)/fs))),2),numSeries,numSeries);
end

% Calculate H (system transfer function) from the inverse of the PDC above
H=zeros(size(tmp_pdc));

for i=1:nFreqs
    H(:,:,i)=inv(tmp_pdc(:,:,i));
end

% Calculate theta^2

theta=zeros(size(H));

for i=1:size(H,1)
    for j=1:size(H,2)
        for k=1:size(H,3)
            theta(i,j,k)=abs(H(i,j,k))^2;
        end
    end
end

if ~normalize
    connectivity=theta;
    return
end

% Calculate gamma^2

gamma=zeros(size(H));

for k=1:size(H,3)
    for i=1:size(H,1)
        den=sum(abs(H(i,:,k)).^2);
        for j=1:size(H,2)
            gamma(i,j,k)=theta(i,j,k)/den;
        end
    end
end

connectivity=gamma;

end