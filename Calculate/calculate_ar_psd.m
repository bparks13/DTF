function [pxx]=calculate_ar_psd(ar,freqRange,fs)
%% [pxx]=calculate_ar_psd(ar,freqRange,fs)
%
%  Given the AR coefficients, calculate the PSD
%
%   Inputs:
%    - ar: Vector of AR coefficients. Size is [m x 1], where m is the model order
%    - freqRange: Vector of frequencies to calculate the PSD over
%    - fs: Sampling frequency in Hz
%
%   Outputs:
%    - pxx: Estimated power 
%
%   See also: plot_psd, estimate_ar_coefficients, calculate_ar_spectra
%

pxx=nan;
return

if size(ar,1)==1
    ar=squeeze(ar);
end

nFreqs=length(freqRange);
modelOrder=length(ar);

pxx=nan(nFreqs,1);
    
for i=1:nFreqs
    pxx(i)=abs(1-[1;ar]'*exp(-2*pi*1i*(i/fs)*(0:modelOrder)'));
end

end