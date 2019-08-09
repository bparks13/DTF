function [num,den]=CreateBSF_butter(fs,order,cutoff,visualize)
%% [num,den]=CreateBSF_butter(fs,order,cutoff,visualize)
%
% Creates the numerator and denominator coefficients for a band-stop Butterworth filter,
% based on the given sampling frequency, filter order, and cutoff frequencies
%
%   Inputs: 
%    - fs: Sampling frequency [Hz]
%    - order: Filter order
%    - cutoff: Cutoff frequencies [Hz], 1x2 vector of frequencies 
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
Wn=cutoff./Ny;   % Cutoff, between 0 and 1.0

[num,den]=butter(order,Wn,'stop');

if visualize || nargout==0
    freqz(num,den)
end


end