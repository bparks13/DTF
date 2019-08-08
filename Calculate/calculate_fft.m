function [P1,spectral_range]=calculate_fft(x,fs,norm)
%% [P1,spectral_range]=calculate_fft(x,fs,norm)
%
%  Calculate a simple FFT of the given signal, at the maximum frequency resolution. Return
%  the FFT values as well as the frequency values as determined by the sampling frequency
%
%   Inputs:
%    - x: Signal to calculate FFT for. Size is [n x s], where n is the number of samples
%       and s is the number of series
%    - fs: Sampling frequency in Hz
%    - norm: Optional input, if true will normalize the FFT by the number of samples
%
%   Outputs:
%    - P1: One-sided power spectrum as calculated by FFT, normalized by the number of
%       samples. Size is [f x s], where f is equal to N/2, and s is the number of series
%    - spectral_range: Vector of frequencies over which the FFT was calculated
%    
%  See also: calculate_ar_spectra, mvar
%

if nargin==2
    norm=true;
end

N=length(x);
                
spectral_range=fs*(0:(N/2))/N;

Y=fft(x);

if norm
    P2=abs(Y/N);
else
    P2=abs(Y);
end

P1=P2(1:floor(N/2)+1,:);
P1(2:end,:)=2*P1(2:end,:);

end