function plot_connectivity(connectivity,series,freqRange,labels,config)
%% plot_connectivity(connectivity,series,freqRange,labels,config)
%
%  Given the DTF values of connectivity, plot the connectivities across all frequencies in
%  the off-diagonal subplots, and the individual power spectra of the series on the diagonal.
%
%   Inputs:
%    - connectivity: DTF values of connectivity (either gamma (normalized) or theta) for
%       all series
%    - series: Either the original signal used, AR coefficients that are estimated, or a PSD;
%       these will then be used to calculate the PSD and plotted on the diagonals. Which
%       input is used is defined in config, or the original signal is assumed to be given
%    - freqRange: Vector of the range of frequencies over which the conenctivity is measured
%    - labels: Labels of the series, used for the titles to indicate the directionality of
%       the connections. Should bea cell array of strings corresponding to the signals used.
%    - config: Optional struct containing additional parameters
%       seriesType: Int defining which type of series is given. 1 (original signal
%           [default]), 2 (AR coefficients), 3 (pre-calculated PSD)
%       fs: Sampling frequency in Hz. If not defined, assumed to be 2400 Hz
%       hFig: Handle to an existing figure, does not create a new figure anymore
%
%   Outputs:
%    - Figure containing subplots with PSD on the diagonal, and DTF connectivity
%       measurements in the off-diagonal plots
%
%  See also: varm, estimate, mvar, dtf
%

fs=2400;
bool_calcOrigPSD=true;
bool_calcARPSD=false;
bool_newFig=true;

if nargin > 4 && isstruct(config)
    if isfield(config,'fs')
        fs=config.fs;
    end
    
    if isfield(config,'seriesType')
        if config.seriesType == 1
            bool_calcOrigPSD=true;
            bool_calcARPSD=false;
        elseif config.seriesType == 2
            bool_calcOrigPSD=false;
            bool_calcARPSD=true;
        elseif config.seriesType == 3
            bool_calcOrigPSD=false;
            bool_calcARPSD=false;
            psd=series;
        end
    end
    
    if isfield(config,'hFig')
        bool_newFig=false;
    end
end

if bool_calcOrigPSD
    window=round(fs);
    overlap=round(window/2);
    
    psd=pwelch(series,window,overlap,freqRange,fs);
elseif bool_calcARPSD
    disp('Unable to calculate PSD from AR coefficients currently');
    psd=[];
end

if bool_newFig
    figure;
else
    figure(config.hFig);
end

numSeries=size(connectivity,1);

ax_diag=zeros(numSeries);
ax_offdiag=zeros(numSeries * numSeries - numSeries);

for i=1:numSeries
    for j=1:numSeries
        currSubPlot=(i-1)*numSeries+j;
        
        if i == j
            ax_diag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
            plot(freqRange,10*log10(psd(:,i)),'k'); title(labels{i}); hold on;
        elseif i ~= j 
            ax_offdiag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
            plot(freqRange,squeeze(connectivity(i,j,:)),'b'); hold on;
            title(sprintf('%s %c %s',labels{i},8594,labels{j}));
        end
    end
end

linkaxes(ax_diag); xlim([freqRange(1) freqRange(end)]); 
subplot(numSeries,numSeries,numSeries); linkaxes(ax_offdiag); 
xlim([freqRange(1) freqRange(end)]); ylim([0 1]);

end