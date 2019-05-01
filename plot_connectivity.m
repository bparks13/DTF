function plot_connectivity(connectivity,freqRange,labels)
%% plot_connectivity(connectivity,freqRange,labels)
%
%  Given the DTF values of connectivity, plot the connectivities across all frequencies in
%  the off-diagonal subplots, and the individual power spectra of the series on the diagonal.
%
%   Inputs:
%    - connectivity: DTF values of connectivity (either gamma (normalized) or theta) for
%       all series
%    - freqRange: Vector of the range of frequencies over which the conenctivity is measured
%    - labels: Labels of the series, used for the titles to indicate the directionality of
%       the connections. Should bea cell array of strings corresponding to the signals used.
%
%   Outputs:
%    - Figure containing subplots with PSD on the diagonal, and DTF connectivity
%       measurements in the off-diagonal plots
%
%  See also: varm, estimate, mvar, dtf
%

figure;

numSeries=size(connectivity,1);

ax_diag=zeros(numSeries);
ax_offdiag=zeros(numSeries * numSeries - numSeries);

for i=1:numSeries
    for j=1:numSeries
        currSubPlot=(i-1)*numSeries+j;
        
        if i == j
            ax_diag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
            title(labels{i});
            % Plot the PSD's here
        elseif i ~= j 
            ax_offdiag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
            plot(freqRange,squeeze(connectivity(i,j,:)));
            title(sprintf('%s %c %s',labels{i},8594,labels{j}));
        end
    end
end

linkaxes(ax_diag); xlim([freqRange(1) freqRange(end)]); ylim([-40 10]);
subplot(numSeries,numSeries,numSeries); linkaxes(ax_offdiag); 
xlim([freqRange(1) freqRange(end)]); ylim([0 1]);

end