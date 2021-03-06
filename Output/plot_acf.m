function plot_acf(x,lag,numSD,hAx)
%% plot_acf(x,lag,numSD,hAx)
%
%  Function to plot the autocorrelation up to the lag given. Based on the 'autocorr'
%  function in Matlab, reiterated here to bypass the licensing issue. If a matrix is
%  given, each channel will be plotted in a different figure
%
%   Inputs:
%    - x: Signal to autocorrelate. Can be a vector or a matrix. Size is [n x c], where n
%       is the number of samples and c is the number of channels
%    - lag: Max lag to plot the ACF up to. Default is 20. Can be empty to use default
%    - numSD: Number of standard deviations to plot for the confidence interval. Default
%       is 2. Can be empty to use default
%    - hAx: Optional input containing the handle to an axis to plot in, instead of
%       creating a new figure. This input is ignored if there is a matrix of values given
%
%   Outputs:
%    Figure(s) containing the autocorrelation values at lags up to lag, with confidence
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
bool_newAxis=true;

if nargin==1
    
elseif nargin==2
    maxLag=lag;
elseif nargin==3
    maxLag=lag;
    numSTD=numSD;
elseif nargin==4
    if ~isempty(lag)
        maxLag=lag;
    end
    
    if ~isempty(numSD)
        numSTD=numSD;
    end
    
    if ~isempty(hAx)
        bool_newAxis=false;
    end
end

if size(x,2) > 1
    for i=1:size(x,2)
        plot_acf(x(:,i),maxLag,numSTD);
    end
    
    return
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

if bool_newAxis
    figure;
    lineHandles = stem(lags,acf,'filled','r-o'); hold on;
    plot([numMA+0.5 numMA+0.5; numLags numLags],[bounds([1 1]) bounds([2 2])],'-b');
    plot([0 numLags],[0 0],'-k');
else
    lineHandles = stem(hAx,lags,acf,'filled','r-o'); hold on;
    plot(hAx,[numMA+0.5 numMA+0.5; numLags numLags],[bounds([1 1]) bounds([2 2])],'-b');
    plot(hAx,[0 numLags],[0 0],'-k');
end

set(lineHandles(1),'MarkerSize',4)
grid('on')
xlabel('Lag')
ylabel('Sample Autocorrelation')
title('Sample Autocorrelation Function')
xlim([0 numLags]);

end