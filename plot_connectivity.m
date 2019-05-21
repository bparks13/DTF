function [avgPSD,avgConnectivity,stdPSD,stdConnectivity]=plot_connectivity(connectivity,series,freqRange,labels,config)
%% [avgPSD,,avgConnectivity,stdPSD,stdConnectivity]=plot_connectivity(connectivity,series,freqRange,labels,config)
%
%  Given the DTF values of connectivity, plot the connectivities across all frequencies in
%  the off-diagonal subplots, and the individual power spectra of the series on the diagonal.
%
%   Inputs:
%    - connectivity: DTF values of connectivity (either gamma (normalized) or theta) for
%       all series. Size is [s x s x f x t], where s is the number of series, f is the
%       number of frequencies being analyzed, and t is the number of trials
%    - series: Either the original signal used, AR coefficients that are estimated, or a PSD;
%       these will then be used to calculate the PSD and plotted on the diagonals. Which
%       input is used is defined in config, or the original signal is assumed to be given.
%       For the signal, size is [n x s x t], where n is the number of samples, s is the
%       number of series, and t is the number of trials. 
%       For the AR coefficients, should be a struct that contains the mdl field, from
%       which the coefficients can be extracted. Base struct should have length equal to
%       the number of trials, and each AR subfield should be [s x s x m], where s is the
%       number of series, and m is the model order
%       For the PSD, size is [f x s x t], where f is the number of frequencies being
%       analyzed, s is the number of series, and t is the number of trials
%    - freqRange: Vector of the range of frequencies over which the connectivity is measured
%    - labels: Labels of the series, used for the titles to indicate the directionality of
%       the connections. Should bea cell array of strings corresponding to the signals used.
%    - config: Optional struct containing additional parameters
%       seriesType: Int defining which type of series is given. 1 (original signal
%           [default]), 2 (AR coefficients), 3 (pre-calculated PSD)
%       fs: Sampling frequency in Hz. If not defined, assumed to be 2400 Hz
%       hFig: Handle to an existing figure, does not create a new figure anymore
%       figTitle: String containing a figure title, containing for example the condition
%           being tested, or the patient/date/run combo, or all of the above
%       plotType: String denoting what to plot. 'ind' plot all individual traces, 'avg'
%           plot averages, 'avgerr' plot averages with shaded error bars. Default behavior
%           is to plot individual traces. If plotType is 'avg' or 'avgerr', returns the
%           average PSDs and the average gamma values if there are output variables to
%           give them to 
%       h: Only applicable for plotType = 'ind'; denotes whether the null hypothesis of no
%           autocorrelation among the residuals is kept (h = 0) or rejected (h = 1). Plots
%           individual traces where the null hypothesis is rejected in red. Size is 
%           [t x s], where t is the number of trials, and s is the number of series
%
%   Outputs:
%    - Figure containing subplots with PSD on the diagonal, and DTF connectivity
%       measurements in the off-diagonal plots
%    - avgPSD: If plotType is 'avg' or 'avgerr', returns the average PSD calculated
%    - avgConnectivity: If plotType is 'avg' or 'avgerr', returns the average gamma values
%       calculated 
%    - stdPSD: If plotType is 'avgerr', returns the standard deviation of the PSDs
%    - stdConnectivity: If plotType is 'avgerr', returns the standard deviation of the
%       connectivities
%
%  See also: varm, estimate, mvar, dtf, estimate_ar_coefficients, estimate_residuals
%

fs=2400;
bool_calcOrigPSD=true;
bool_calcARPSD=false;
bool_newFig=true;
plotType='ind';
figTitle='';
bool_showRejectedNull=false; % whether or not to plot the rejected null hypothesis trials in red

numSeries=size(connectivity,1);
numTrials=size(connectivity,4);

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
            pxx=series;
        end
    end
    
    if isfield(config,'hFig')
        bool_newFig=false;
    end
    
    if isfield(config,'figTitle')
        figTitle=config.figTitle;
    end
    
    if isfield(config,'plotType')
        plotType=config.plotType;
    end
    
    if isfield(config,'h') && strcmp(plotType,'ind')
        bool_showRejectedNull=true;
        h=config.h;
    end
end

% Calculate the PSD from the original signal 

if bool_calcOrigPSD
    window=round(fs);
    overlap=round(window/2);
    
    pxx=nan(length(freqRange),numSeries,numTrials);
    
    for i=1:numSeries
        for j=1:numTrials
            pxx(:,i,j)=pwelch(series(:,i,j),window,overlap,freqRange,fs);
        end
    end
elseif bool_calcARPSD
    pxx=nan(length(freqRange),numSeries,numTrials);
    
    for i=1:numSeries
        for j=1:numTrials
            pxx(:,i,j)=calculate_ar_psd(series(j).mdl.AR(i,i,:),freqRange,fs);
        end
    end
end

% Create a new figure, or focus on an existing one

if bool_newFig
    figure;
else
    figure(config.hFig);
end

% If averages need to be calculated

if strcmp(plotType,'avg')
    avgPSD=mean(pxx,3);
    avgConnectivity=mean(connectivity,4);
elseif strcmp(plotType,'avgerr')
    avgPSD=mean(pxx,3);
    avgConnectivity=mean(connectivity,4);
    stdPSD=std(10*log10(pxx),0,3);  
    % Standard deviation is the standard deviation of the logarithmic power, not the raw
    % power
    stdConnectivity=std(connectivity,0,4);    
end

% Plotting

if strcmp(plotType,'ind')
    for i=1:numTrials
        [ax_diag,ax_offdiag]=plotting(pxx(:,:,i),connectivity(:,:,:,i),[],[],h(i,:));
    end
    yLimits=[min(min(10*log10(pxx)))-5,max(max(10*log10(pxx)))+5];
    set_figure(ax_diag,ax_offdiag,yLimits);
elseif strcmp(plotType,'avg')
    [ax_diag,ax_offdiag]=plotting(avgPSD,avgConnectivity);
    set_figure(ax_diag,ax_offdiag,[min(min(10*log10(avgPSD)))-5,max(max(10*log10(avgPSD)))+5]);
elseif strcmp(plotType,'avgerr')
    [ax_diag,ax_offdiag]=plotting(avgPSD,avgConnectivity,stdPSD,stdConnectivity);
    set_figure(ax_diag,ax_offdiag,[min(min(10*log10(avgPSD)))-5,max(max(10*log10(avgPSD)))+5]);
end

if ~isempty(figTitle)
    tmp_hFig=gcf;
    tmp_hFig.Name=figTitle;
end

drawnow;

%% Internal plotting function to deal with a single trials worth of signals

    function [ax_diag,ax_offdiag]=plotting(pxx,conn,pxx_std,conn_std,h)
        
        if nargin == 4
            bool_plotErrorBars=true;
        else
            bool_plotErrorBars=false;
        end
        
        ax_diag=zeros(numSeries);
        ax_offdiag=zeros(numSeries * numSeries - numSeries);

        for k=1:numSeries
            for l=1:numSeries
                currSubPlot=(k-1)*numSeries+l;

                if k == l
                    ax_diag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
                    
                    if bool_plotErrorBars
                        shadedErrorBar(freqRange,10*log10(pxx(:,k)),pxx_std(:,k));
                    else
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,10*log10(pxx(:,k)),'r'); 
                        else
                            plot(freqRange,10*log10(pxx(:,k)),'k'); 
                        end
                    end
                    
                    title(labels{k}); hold on;
                elseif k ~= l 
                    ax_offdiag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
                    
                    if bool_plotErrorBars
                        shadedErrorBar(freqRange,conn(k,l,:),conn_std(k,l,:));
                    else
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,squeeze(conn(k,l,:)),'r'); 
                        else
                            plot(freqRange,squeeze(conn(k,l,:)),'b'); 
                        end
                    end
                    
                    hold on;
                    title(sprintf('%s %c %s',labels{k},8594,labels{l}));
                end
            end
        end
    end

%% Internal function to linkaxes and set limits for the subplots
    function set_figure(ax_diag,ax_offdiag,yLimits)
       
        linkaxes(ax_diag); 
        xlim([freqRange(1) freqRange(end)]); 
        ylim([yLimits(1),yLimits(2)])
        subplot(numSeries,numSeries,numSeries); 
        
        linkaxes(ax_offdiag); 
        xlim([freqRange(1) freqRange(end)]); 
        ylim([0 1]);

    end
end