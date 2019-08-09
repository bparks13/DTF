function [num,den]=CreateLPF_butter(fs,order,cutoff,visualize)
%% [num,den]=CreateLPF_butter(fs,order,cutoff,visualize)
%
% Creates the numerator and denominator coefficients for a low pass Butterworth filter,
% based on the given sampling frequency, filter order, and cutoff frequency
%
%   Inputs: 
%    - fs: Sampling frequency [Hz]
%    - order: Filter order
%    - cutoff: Cutoff frequency [Hz]
%    - visualize: If nonempty, plots the magnitude and phase of the filter
%
%   Outputs:
%    - num: Numerator coefficients
%    - den: Denominator coefficients
%

if nargin<4
    visualize=false;
elseif nargin==4
    visualize=true;
end

Ny=fs/2;        % Nyquist
Wn=cutoff/Ny;   % Cutoff, between 0 and 1.0

[num,den]=butter(order,Wn,'low');

if visualize
    freqz(num,den)
end


end