function plot_acf(x,lag,numSD)
%% plot_acf(x,lag,numSD)
%
%  Function to plot the autocorrelation up to the lag given. Based on the 'autocorr'
%  function in Matlab, reiterated here to bypass the licensing issue.
%
%   Inputs:
%    - x: Signal to autocorrelate. Should be a vector of values.
%    - lag: Max lag to plot the ACF up to. Default is 20
%    - numSD: Number of standard deviations to plot for the confidence interval. Default
%       is 2
%
%   Outputs:
%    Figure containing the autocorrelation values at lags up to lag, with confidence
%       intervals defined by numSD
%
%   See also: autocorr
%

if size(x,1) == 1
    x=squeeze(x);
end

maxLag=20;
numSTD=2;
numMA=0;

if nargin==2
    maxLag=lag;
elseif nargin==3
    maxLag=lag;
    numSTD=numSD;
end

N=length(x);
lags=(0:maxLag)';
numLags=length(lags);

nFFT = 2^(nextpow2(N)+1);
F = fft(x-mean(x),nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf(1:(maxLag+1)); % Retain non-negative lags
acf = acf./acf(1); % Normalize
acf = real(acf);

sigmaNMA = sqrt((1+2*(acf(2:numMA+1)'*acf(2:numMA+1)))/N);  
bounds = sigmaNMA*[numSTD;-numSTD];

figure;

lineHandles = stem(lags,acf,'filled','r-o');
set(lineHandles(1),'MarkerSize',4)
grid('on')
xlabel('Lag')
ylabel('Sample Autocorrelation')
title('Sample Autocorrelation Function')
hold('on')
plot([numMA+0.5 numMA+0.5; numLags numLags],[bounds([1 1]) bounds([2 2])],'-b');
plot([0 numLags],[0 0],'-k');
hold('off')
% a = axis;
% axis([a(1:3) 1]);
xlim([0 numLags]);

end