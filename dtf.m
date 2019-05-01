function [connectivity]=dtf(estMdl,freqRange,fs,config)
%% [connectivity]=dtf(estMdl,freqRange,fs,config)
%
%  Using the estimated AR model, calculate the connectivity of the series in question. Can
%  return either the raw connections or the normalized connections as specified in config
%
%   Inputs:
%    - estMdl: Estimated AR model struct as returned by varm and estimate
%    - freqRange: Vector containing the specific frequencies to calculate the connectivity
%       over
%    - fs: Sampling frequency in Hz
%    - config: Struct containing additional optional parameters
%       normalize: Boolean defining whether to normalize [default] or not
%
%   Outputs:
%    - connectivity: Directed transfer function values for all combinations of series and
%       frequencies; either normalized [default] or not, depending on if config is set
%
%  See also: varm, estimate, mvar
%

normalize=true;

if nargin > 3 && isstruct(config)
    if isfield(config,'normalize')
        if islogical(config.normalize)
            normalize=config.normalize;
        end
    end
end

nFreqs=length(freqRange);
numSeries=estMdl.NumSeries;
modelOrder=estMdl.P;
I=eye(numSeries);

AR=zeros([numSeries,numSeries,modelOrder]);  % Extract the coefficients from the model, [numseries x numseries x modelorder]
for i=1:length(estMdl.AR)
    AR(:,:,i)=estMdl.AR{i};
end

tmp_pdc=zeros(numSeries,numSeries,nFreqs);

% Calculate PDC
for i=1:nFreqs
    tmp_pdc(:,:,i)=I-reshape(sum(bsxfun(@times,reshape(AR,numSeries^2,modelOrder),...
        exp(-(2*pi*1i/fs)*(1:modelOrder)*freqRange(i))),2),numSeries,numSeries);
    % This is the matrix calculation from SIFT
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
            gamma(i,j,k)=abs(H(i,j,k))^2/den;
        end
    end
end

connectivity=gamma;

end