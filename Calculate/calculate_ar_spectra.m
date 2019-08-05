function S=calculate_ar_spectra(AR,freqRange,fs,C,norm)
%% [S]=calculate_ar_psd(AR,freqRange,fs,C,norm)
%
%  Given the AR coefficients, calculate the normalized spectra. Normalized based on the
%  number of frequencies that are calculated. 
%
%   Inputs:
%    - AR: 3D matrix of autoregressive coefficient values. Size is [c x c x o], where
%       c is the number of channels, and o is the model order
%    - freqRange: Vector of frequencies to calculate the PSD over
%    - fs: Sampling frequency in Hz
%    - C: Covariance matrix of the AR model
%    - norm: Optional bool input, defining whether or not to normalize the spectra by the
%       number of frequencies calculated
%
%   Outputs:
%    - S: Estimated power 
%
%   See also: plot_psd, estimate_ar_coefficients
%

if nargin==4
    norm=false;
end

nFreqs=length(freqRange);
numSeries=size(AR,1);
I=eye(numSeries);
modelOrder=size(AR,3);

tmp_pdc=nan(numSeries,numSeries,nFreqs);

% Calculate PDC
for i=1:nFreqs
    tmp_pdc(:,:,i)=I-reshape(sum(bsxfun(@times,reshape(AR,numSeries^2,modelOrder),...
        exp(-(2*pi*1i/fs)*(1:modelOrder)*freqRange(i))),2),numSeries,numSeries);
end

% Calculate H (system transfer function) from the inverse of the PDC above
H=nan(size(tmp_pdc));

for i=1:nFreqs
    H(:,:,i)=inv(tmp_pdc(:,:,i));
end

% Calculate S (spectra)

S=nan(size(H));

for i=1:nFreqs
    S(:,:,i)=H(:,:,i) * C * H(:,:,i)';
end

if norm
    S=S./(sqrt(nFreqs) * numSeries);
end

end