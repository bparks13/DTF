function [avgPSD,avgConn,stdPSD,stdConn]=plot_connectivity(conn,series,freqRange,labels,config)
%% [avgPSD,avgConn,stdPSD,stdConn]=plot_connectivity(conn,series,freqRange,labels,config)
%
%  Given the DTF values of connectivity, plot the connectivities across all frequencies in
%  the off-diagonal subplots, and the individual power spectra of the series on the diagonal.
%
%   Inputs:
%    - conn: DTF values of connectivity (either gamma (normalized) or theta) for
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
%       To plot both the PSD of the original signal and the estimated signal, series can
%       be given as a struct, with the following fields;
%           signal: A matrix that follows the format of the signal as above
%           estimated: A struct that contains the subfield 'mdl', with t trials worth of
%               this field. Within 'mdl', there should be both a matrix 'x_hat' and a
%               matrix 'pxx', corresponding to the estimated signal and estimated power of
%               the estimated signal, respectively. Based on options chosen in config, the
%               pre-calculated values can be used, or will be recalculated based on
%               additional parameters given
%    - freqRange: Vector of the range of frequencies over which the connectivity is measured
%    - labels: Labels of the series, used for the titles to indicate the directionality of
%       the connections. Should bea cell array of strings corresponding to the signals used.
%    - config: Optional struct containing additional parameters
%       seriesType: Int defining which type of series is given. 1 (signal [default]), 2
%           (AR coefficients), 3 (pre-calculated PSD). If series is given as a struct
%           containing both the original and estimated signal, then both will be treated
%           the same way; either both will use the signal given to calculate pxx, or the
%           given pxx will be used for plotting
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
%       freqLims: If defined, changes the limits of the plots to only show the frequency
%           range specified. Can be a vector of the entire range, or just the min and max
%           values to be used
%
%   Outputs:
%    - Figure containing subplots with PSD on the diagonal, and DTF connectivity
%       measurements in the off-diagonal plots
%    - avgPSD: If plotType is 'avg' or 'avgerr', returns the average PSD calculated
%    - avgConn: If plotType is 'avg' or 'avgerr', returns the average gamma values
%       calculated 
%    - stdPSD: If plotType is 'avgerr', returns the standard deviation of the PSDs
%    - stdConn: If plotType is 'avgerr', returns the standard deviation of the
%       connectivities
%
%  See also: varm, estimate, mvar, dtf, estimate_ar_coefficients, estimate_residuals
%

% TODO: Add an option to change the pwelch parameters in config (i.e. give specific
%           window, overlap, etc. values) 
%       Add an option to plot the original signal traces instead of PSD
%       Add a legend to the diagonal plots when both the original and estimated signal PSD
%           are plotted

fs=2400;
bool_newFig=true;
bool_changeFreqRange=false;
plotType='ind';
figTitle='';
seriesType=1;

bool_showRejectedNull=false; % whether or not to plot the rejected null hypothesis trials in red

numSeries=size(conn,1);
numTrials=size(conn,4);

h=zeros(numTrials,numSeries);

freqLims=[];

if nargin > 4 && isstruct(config)
    if isfield(config,'fs')
        fs=config.fs;
    end
    
    if isfield(config,'seriesType')
        seriesType=config.seriesType;
    end
    
    if isfield(config,'hFig')
        if ~isempty(config.hFig)
            bool_newFig=false;
        end
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
    
    if isfield(config,'freqLims')
        bool_changeFreqRange=true;
        freqLims=[config.freqLims(1),config.freqLims(end)];
    end
end

% Create a new figure, or focus on an existing one

if bool_newFig
    figure;
else
    figure(config.hFig);
end

% Calculate the PSDs, or grab them from the given variables

window=round(fs);
overlap=round(window/2);
    
if isstruct(series)
    if seriesType == 1
        pxx_sig=nan(length(freqRange),numSeries,numTrials);
        pxx_est=nan(length(freqRange),numSeries,numTrials);
    
        for i=1:numSeries
            for j=1:numTrials
                pxx_sig(:,i,j)=pwelch(series.signal(:,i,j),window,overlap,freqRange,fs);
                pxx_est(:,i,j)=pwelch(series.estimated(j).mdl.x_hat(:,i),window,overlap,freqRange,fs);
            end
        end
    elseif seriesType == 2
        pxx_sig=nan(length(freqRange),numSeries,numTrials);
        pxx_est=nan(length(freqRange),numSeries,numTrials);
    
        for i=1:numSeries
            for j=1:numTrials
                pxx_sig(:,i,j)=pwelch(series.signal(:,i,j),window,overlap,freqRange,fs);
                pxx_est(:,i,j)=calculate_ar_psd(series.estimated(j).mdl.AR(i,i,:),freqRange,fs);
            end
        end
    elseif seriesType == 3
        pxx_sig=series.signal;
        pxx_est=nan(length(freqRange),numSeries,numTrials);
        
        for i=1:numSeries
            for j=1:numTrials
                if ~isempty(series.esimated(j).mdl.pxx)
                    pxx_est(:,i,j)=series.esimated(j).mdl.pxx(:,i);
                else
                    pxx_est(:,i,j)=pwelch(series.estimated(j).mdl.x_hat(:,i),window,overlap,freqRange,fs);
                end
            end
        end
    end
    
    % Calculate averages and standard deviations if needed, then plot everything
    
    if strcmp(plotType,'ind')
        for i=1:numTrials
            [~,~]=plotting(pxx_sig(:,:,i),conn(:,:,:,i),[],[],h(i,:),'-b');
            [ax_diag,ax_offdiag]=plotting(pxx_est(:,:,i),[],[],[],h(i,:),'-g');
        end

        yLimits=[min(min(10*log10(pxx_sig)))-5,max(max(10*log10(pxx_sig)))+5];

        if bool_changeFreqRange
            set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
        else
            set_figure(ax_diag,ax_offdiag,yLimits);
        end
    elseif strcmp(plotType,'avg')
        avgPSD_sig=mean(pxx_sig,3);
        avgPSD_est=mean(pxx_est,3);
        avgConn=mean(conn,4);
        
        [~,~]=plotting(avgPSD_sig,avgConn,[],[],[],'-b');
        [ax_diag,ax_offdiag]=plotting(avgPSD_est,[],[],[],[],'-g');
        
        yLimits=[min(min(10*log10(avgPSD_sig)))-5,max(max(10*log10(avgPSD_sig)))+5];
    
        if bool_changeFreqRange
            set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
        else
            set_figure(ax_diag,ax_offdiag,yLimits);
        end
    elseif strcmp(plotType,'avgerr')
        avgPSD_sig=mean(pxx_sig,3);
        stdPSD_sig=std(10*log10(pxx_sig),0,3);  
        avgPSD_est=mean(pxx_est,3);
        stdPSD_est=std(10*log10(pxx_est),0,3);  
        avgConn=mean(conn,4);
        stdConn=std(conn,0,4); 
        
        [~,~]=plotting(avgPSD_sig,avgConn,stdPSD_sig,stdConn,[],'-b');
        [ax_diag,ax_offdiag]=plotting(avgPSD_est,[],stdPSD_est,[],[],'-g');
        
        yLimits=[min(min(10*log10(avgPSD_sig)))-5,max(max(10*log10(avgPSD_sig)))+5];
    
        if bool_changeFreqRange
            set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
        else
            set_figure(ax_diag,ax_offdiag,yLimits);
        end   
    end
else
    if seriesType == 1
        pxx_sig=nan(length(freqRange),numSeries,numTrials);
    
        for i=1:numSeries
            for j=1:numTrials
                pxx_sig(:,i,j)=pwelch(series(:,i,j),window,overlap,freqRange,fs);
            end
        end
    elseif seriesType == 3
        pxx_sig=series;
    end
    
    % Calculate averages and standard deviations if needed, then plot everything
    
    if strcmp(plotType,'ind')
        for i=1:numTrials
            [ax_diag,ax_offdiag]=plotting(pxx_sig(:,:,i),conn(:,:,:,i),[],[],h(i,:));
        end

        yLimits=[min(min(10*log10(pxx_sig)))-5,max(max(10*log10(pxx_sig)))+5];

        if bool_changeFreqRange
            set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
        else
            set_figure(ax_diag,ax_offdiag,yLimits);
        end
    elseif strcmp(plotType,'avg')
        avgPSD=mean(pxx_sig,3);
        avgConn=mean(conn,4);
        
        [ax_diag,ax_offdiag]=plotting(avgPSD,avgConn);
        
        yLimits=[min(min(10*log10(avgPSD)))-5,max(max(10*log10(avgPSD)))+5];
    
        if bool_changeFreqRange
            set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
        else
            set_figure(ax_diag,ax_offdiag,yLimits);
        end
    elseif strcmp(plotType,'avgerr')
        avgPSD=mean(pxx_sig,3);
        avgConn=mean(conn,4);
        stdPSD=std(10*log10(pxx_sig),0,3);  
        stdConn=std(conn,0,4);  
        
        [ax_diag,ax_offdiag]=plotting(avgPSD,avgConn,stdPSD,stdConn);
        
        yLimits=[min(min(10*log10(avgPSD)))-5,max(max(10*log10(avgPSD)))+5];
    
        if bool_changeFreqRange
            set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
        else
            set_figure(ax_diag,ax_offdiag,yLimits);
        end
    end
end

if ~isempty(figTitle)
    tmp_hFig=gcf;
    tmp_hFig.Name=figTitle;
end

drawnow;

if nargout==0
    clear avgPSD avgConn stdPSD stdConn
elseif nargout == 2
    clear stdPSD stdConn
end

%% Internal plotting function to deal with a single trials worth of signals

    function [ax_diag,ax_offdiag]=plotting(pxx,conn,pxx_std,conn_std,h,lineProps)
        
        if nargin == 4 || (nargin == 6 && isempty(h) && ~isempty(pxx_std))
            bool_plotErrorBars=true;
        else
            bool_plotErrorBars=false;
        end
        
        if nargin < 6
            lineProps='-k';
        end
        
        if isempty(conn)
            conn=nan(numSeries,numSeries,length(freqRange));
            
            if isempty(conn_std)
                conn_std=conn;
            end
        end
        
        ax_diag=zeros(numSeries);
        ax_offdiag=zeros(numSeries * numSeries - numSeries);

        for k=1:numSeries
            for l=1:numSeries
                currSubPlot=(k-1)*numSeries+l;

                if k == l
                    ax_diag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
                    
                    if bool_plotErrorBars
                        shadedErrorBar(freqRange,10*log10(pxx(:,k)),pxx_std(:,k),'lineProps',lineProps);
                    else
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,10*log10(pxx(:,k)),'r'); 
                        elseif ~strcmp(lineProps,'-k')
                            plot(freqRange,10*log10(pxx(:,k)),lineProps); 
                        else
                            plot(freqRange,10*log10(pxx(:,k)),'k'); 
                        end
                    end
                    
                    title(labels{k}); hold on;
                elseif k ~= l 
                    ax_offdiag(currSubPlot)=subplot(numSeries,numSeries,currSubPlot);
                    
                    if bool_plotErrorBars
                        shadedErrorBar(freqRange,conn(k,l,:),conn_std(k,l,:),'lineProps',lineProps);
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
    function set_figure(ax_diag,ax_offdiag,yLimits,xLimits)
        
        if nargin==3
            xLimits=[freqRange(1) freqRange(end)];
        end
        
        linkaxes(ax_diag); 
        xlim([xLimits(1), xLimits(2)]); 
        ylim([yLimits(1),yLimits(2)])
        subplot(numSeries,numSeries,numSeries); 

        linkaxes(ax_offdiag); 
        xlim([xLimits(1), xLimits(2)]); 
        ylim([0 1]);

    end
end