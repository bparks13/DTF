function plot_compare(x_orig,x_filt)
%% plot_compare(x_orig,x_filt)
%
%  Given two signals, plot each channel in its own subplot to compare before and after
%  serially decorrelating the signal
%
%   Inputs:
%    - x_orig: Original signal, for now is assumed to be a single trial, with size
%       [s x c], where s is the number of samples and c is the number of channels
%    - x_filt: Filtered signal, same size as x_orig
%
%   Outputs:
%    Figure containing multiple subplots, with each channel showing both the original and
%       the decorrelated signal
%
%  See also: filter_serial_correlation
%

figure;

numSamples=size(x_orig,1);
numChannels=size(x_orig,2);

t=0:numSamples-1;

for i=1:numChannels
    subplot(numChannels,1,i);
    plot(t,x_orig(:,i),'b',t,x_filt(:,i),'r');
    legend('Original','Filtered');
end

end