function plot_time_series(x,x_hat,e,fs)
%% plot_time_series(x,x_hat,e,fs)
%
%  In a new figure, plots the original signal, the estimated signal, and the residuals for
%  all channels given
%
%   Inputs:
%    - x: Original signal, size is [s x c], where s is the number of samples and c is the
%       number of channels
%    - x_hat: Estimated signal, where size is [(s - o) x c], where o is the model order
%    - e: Residuals, size is [(s - o) x c]
%    - fs: Sampling frequency
%
%   Outputs:
%    Figure contianing c subplots, with each subplot containing the original/estimated
%       signal, and the residuals 
%
%  See also: estimate_residuals

figure;
numChannels=size(x,2);
numSamples=length(x);
t=(0:(numSamples-1))/fs;
colors=linspecer(3);
modelOrder=numSamples-length(e);
ax=zeros(numChannels,1);

if ~isempty(x_hat) && ~isempty(e)
    for i=1:numChannels
        ax(i)=subplot(numChannels,1,i);
        plot(t,x(:,i),'Color',colors(1,:)); hold on;
        plot(t(modelOrder+1:end),x_hat(:,i),'Color',colors(2,:));
        plot(t(modelOrder+1:end),e(:,i),'Color',colors(3,:));
        xlim([t(1) t(end)])
        legend('x','x\_hat','error');
    end
else
    for i=1:numChannels
        ax(i)=subplot(numChannels,1,i);
        plot(t,x(:,i),'Color',colors(1,:)); 
        xlim([t(1) t(end)])
    end    
end

linkaxes(ax,'x');

end