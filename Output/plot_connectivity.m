function [avgPSD,avgConn,stdPSD,stdConn]=plot_connectivity(conn,series,freqRange,labels,config)
%% [avgPSD,avgConn,stdPSD,stdConn]=plot_connectivity(conn,series,freqRange,labels,config)
%
%  Given the DTF values of connectivity, plot the connectivities across all frequencies in
%  the off-diagonal subplots, and the individual power spectra of the channels on the
%  diagonal. 
%
%   Inputs:
%    - conn: DTF values of connectivity (either gamma (normalized) or theta) for
%       all series. Size is [c x c x f x t], where c is the number of channels, f is the
%       number of frequencies being analyzed, and t is the number of trials
%    - series: This input can be a number of things, and is therefore dependent on being
%       defined in config. The following are the possible inputs; 
%           1) Original signal [default]: Matrix assumed to be [n x c x t], where n is the
%              number of samples, c is the number of channels, and t is the number of
%              trials
%           2) Estimated signal: Matrix assumed to be [n x c x t], where n is the
%              number of samples, c is the number of channels, and t is the number of
%              trials 
%           3) PSD of the original signal: Matrix assumed to be [f x c x t], where f is
%              the number of frequencies being analyzed, c is the number of channels, and
%              t is the number of trials  
%           4) PSD of the estimated signal: Matrix assumed to be [f x c x t], where f is
%              the number of frequencies being analyzed, c is the number of channels, and
%              t is the number of trials  
%           5) Surrogate analysis values: Matrix containing all of the surrogate analysis
%              values for a particular condition. Size should be [c x c x f x t], where c
%              is the number of channels, f is the number of frequencies being analyzed,
%              and t is the number of iterations the surrogate analysis was run for
%           6) Struct containing any combination of the inputs above. The possible
%              subfields are given below;
%               'original' is a matrix similar to 1)
%               'estimated' is a matrix similar to 2)
%               'original_psd' is a matrix similar to 3)
%               'estimated_psd' is a matrix similar to 4)
%               'surrogate' is a matrix similar to 5)
%    - freqRange: Vector of the range of frequencies over which the connectivity is
%       measured 
%    - labels: Labels of the series, used for the titles to indicate the directionality of
%       the connections. Should be a cell array of strings corresponding to the signals
%       used. 
%    - config: Optional struct containing additional parameters. Not all of them need to
%       be defined
%           seriesType: Int defining which type of series is given. Number corresponds
%             with the series given above. 
%               NOTE: If series is a struct, '4' does not need to be defined. Rather,
%               seriesType should define how the values inside the struct are to be
%               treated 
%           fs: Sampling frequency in Hz. If not defined, assumed to be 2400 Hz
%           hFig: Handle to an existing figure, does not create a new figure anymore
%           figTitle: String containing a figure title, containing for example the
%             condition being tested, or the patient/date/run combo, or all of the above
%           plotType: String denoting what to plot. 
%             'ind' plot all individual traces [Default]
%             'avg' plot averages
%             'avgerr' plot averages with shaded error bars
%           NOTE: If plotType is 'avg' or 'avgerr', can return the average PSDs and the
%             average gamma values
%           h: Only applicable for plotType = 'ind'; denotes whether the null hypothesis
%             of no autocorrelation among the residuals is kept (h = 0) or rejected
%             (h = 1). Plots individual traces where the null hypothesis is rejected in
%             red. Size is [t x c], where t is the number of trials, and c is the number
%             of channels 
%           freqLims: If defined, changes the limits of the plots to only show the
%             frequency range specified. Can be a vector of the entire range, or just the
%             min and max values to be used
%           surr_params: Additional struct that can be defined to modify how the surrogate
%             data is shown. However, this is not a necessary field, as all defaults are
%             defined
%               threshold: Float defining what percentage of the distribution is
%                 considered to be a significant amount of connection. Default is 0.01, or
%                 1%
%
%   Outputs:
%    - Figure containing subplots with PSD on the diagonal, and DTF connectivity
%       measurements in the off-diagonal plots
%    - avgPSD: If plotType is 'avg' or 'avgerr', can return the average PSD calculated
%    - avgConn: If plotType is 'avg' or 'avgerr', can return the average gamma values
%       calculated 
%    - stdPSD: If plotType is 'avgerr', can return the standard deviation of the PSDs
%    - stdConn: If plotType is 'avgerr', can return the standard deviation of the
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

numChannels=size(conn,1);
numTrials=size(conn,4);

h=zeros(numTrials,numChannels);

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
        pxx_sig=nan(length(freqRange),numChannels,numTrials);
        pxx_est=nan(length(freqRange),numChannels,numTrials);
    
        for i=1:numChannels
            for j=1:numTrials
                pxx_sig(:,i,j)=pwelch(series.signal(:,i,j),window,overlap,freqRange,fs);
                pxx_est(:,i,j)=pwelch(series.estimated(j).mdl.x_hat(:,i),window,overlap,freqRange,fs);
            end
        end
    elseif seriesType == 2
        pxx_sig=nan(length(freqRange),numChannels,numTrials);
        pxx_est=nan(length(freqRange),numChannels,numTrials);
    
        for i=1:numChannels
            for j=1:numTrials
                pxx_sig(:,i,j)=pwelch(series.signal(:,i,j),window,overlap,freqRange,fs);
                pxx_est(:,i,j)=calculate_ar_psd(series.estimated(j).mdl.AR(i,i,:),freqRange,fs);
            end
        end
    elseif seriesType == 3
        pxx_sig=series.signal;
        pxx_est=nan(length(freqRange),numChannels,numTrials);
        
        for i=1:numChannels
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
    if seriesType == 1 && ~isempty(series)
        pxx_sig=nan(length(freqRange),numChannels,numTrials);
        
        for i=1:numChannels
            for j=1:numTrials
                pxx_sig(:,i,j)=pwelch(series(:,i,j),window,overlap,freqRange,fs);
            end
        end
    elseif seriesType == 2 && ~isempty(series)
        pxx_sig=nan(length(freqRange),numChannels,numTrials); % This is the estimated PSD from the AR coefficients. 
    
        for i=1:numChannels
            for j=1:numTrials
                pxx_sig(:,i,j)=calculate_ar_psd(series(j).mdl.AR(i,i,:),freqRange,fs);
            end
        end
    elseif seriesType == 3 && ~isempty(series)
        pxx_sig=series;
    else
        pxx_sig=nan(length(freqRange),numChannels,numTrials);
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
        
        if all(any(isnan(avgPSD)))
            yLimits=ylim;
        else
            yLimits=[min(min(10*log10(avgPSD)))-5,max(max(10*log10(avgPSD)))+5];
        end
        
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
            conn=nan(numChannels,numChannels,length(freqRange));
            
            if isempty(conn_std)
                conn_std=conn;
            end
        end
        
        ax_diag=zeros(numChannels);
        ax_offdiag=zeros(numChannels * numChannels - numChannels);

        for k=1:numChannels
            for l=1:numChannels
                currSubPlot=(k-1)*numChannels+l;

                if k == l
                    ax_diag(currSubPlot)=subplot(numChannels,numChannels,currSubPlot);
                    
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
                    ax_offdiag(currSubPlot)=subplot(numChannels,numChannels,currSubPlot);
                    
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
        xlim([xLimits(1),xLimits(2)]); 
        ylim([yLimits(1),yLimits(2)])
        subplot(numChannels,numChannels,numChannels); 

        linkaxes(ax_offdiag); 
        xlim([xLimits(1), xLimits(2)]); 
        ylim([0 1]);

    end
end