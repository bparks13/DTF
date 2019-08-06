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
%           6) Transfer function values: Matrix containing the output from the dtf.m
%              function. Size is [c x c x f], where c is the number of channels and f is
%              the number of frequencies (must match freqRange given)
%           7) Spectral matrix values: Matrix containing the spectral matrix defined as
%              S = H * C * H'. Size is [c x c x f], where c is the number of channels and
%              f is the number of frequencies (must match freqRange given)
%           8) Struct containing any combination of the inputs above. The possible
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
%       Add ability to specify lineProps in config

fs=2400;
bool_newFig=true;
bool_changeFreqRange=false;
plotType='ind';
figTitle='';
seriesType=1;
threshold=0.01;

bool_showRejectedNull=false; % whether or not to plot the rejected null hypothesis trials in red
bool_plotThreshold=false; % Plot the values given in the surrogate analysis for significance
bool_plotTransferFunction=false; % Plot the transfer function values instead of the Pxx values
bool_plotSpectralMatrix=false; % Plot the spectral matrix values instead of the Pxx values

numChannels=size(conn,1);
numTrials=size(conn,4);
numFrequencies=length(freqRange);

h=zeros(numTrials,numChannels);

freqLims=[];

if nargin > 4 && isstruct(config)
    if isfield(config,'fs')
        fs=config.fs;
    end
    
    if isfield(config,'seriesType')
        seriesType=config.seriesType;
        
        if isstruct(series)
            seriesType=6;
        end
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
    
    if isfield(config,'surr_params')
        if isfield(config.surr_params,'threshold')
            threshold=config.surr_params.threshold;
        end
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

if seriesType == 1 || seriesType == 2
    pxx=nan(numFrequencies,numChannels,numTrials);
    
    for i=1:numChannels
        for j=1:numTrials
            pxx(:,i,j)=pwelch(series(:,i,j),window,overlap,freqRange,fs);
        end
    end
    
    if strcmp(plotType,'ind')
        [ax_diag,ax_offdiag,yLimits]=plot_ind(pxx,conn,h);
    elseif strcmp(plotType,'avg')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(pxx,conn);
    elseif strcmp(plotType,'avgerr')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(pxx,conn);
    end
elseif seriesType == 3 || seriesType == 4
    if strcmp(plotType,'ind')
        [ax_diag,ax_offdiag,yLimits]=plot_ind(series,conn,h);
    elseif strcmp(plotType,'avg')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(series,conn);
    elseif strcmp(plotType,'avgerr')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(series,conn);
    end
elseif seriesType == 5
    bool_plotThreshold=true;
    [ax_diag,ax_offdiag,yLimits]=plot_surrogate(series,threshold);
    bool_plotThreshold=false;
    
    if strcmp(plotType,'ind')
        [~,~,~]=plot_ind([],conn,h);
    elseif strcmp(plotType,'avg')
        [~,~,~,avgPSD,avgConn]=plot_avg([],conn);
    elseif strcmp(plotType,'avgerr')
        [~,~,~,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr([],conn);
    end
elseif seriesType == 6
    bool_plotTransferFunction=true;
    
    if size(series,3) > 1
        series=resize_spectra(series);
    end
    
    if strcmp(plotType,'ind')
        [ax_diag,ax_offdiag,yLimits]=plot_ind(series,conn,h);
    elseif strcmp(plotType,'avg')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(series,conn);
    elseif strcmp(plotType,'avgerr')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(series,conn);
    end
elseif seriesType == 7
    bool_plotSpectralMatrix=true;
    
    if size(series,3) > 1
        series=resize_spectra(series);
    end
    
    if strcmp(plotType,'ind')
        [ax_diag,ax_offdiag,yLimits]=plot_ind(series,conn,h);
    elseif strcmp(plotType,'avg')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(series,conn);
    elseif strcmp(plotType,'avgerr')
        [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(series,conn);
    end
elseif seriesType == 8
    if ~isstruct(series)
        error('ERROR: Series is set as ''struct'' in config, but is not a struct');
    end
    
    yLimits=[];
    
    all_lineprops={'-b','-g','-r','-k','-c','-m','-y'};
    colorNum=2;
    
    if isfield(series,'original')
        pxx=nan(numFrequencies,numChannels,numTrials);
    
        for i=1:numChannels
            for j=1:numTrials
                pxx(:,i,j)=pwelch(series.original(:,i,j),window,overlap,freqRange,fs);
            end
        end

        if strcmp(plotType,'ind')
            [ax_diag,ax_offdiag,yLimits]=plot_ind(pxx,conn,h,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avg')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(pxx,conn,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avgerr')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(pxx,conn,all_lineprops{1},all_lineprops{colorNum});
        end
        
        colorNum=colorNum+1;
    end
    
    if isfield(series,'estimated')
        pxx=nan(numFrequencies,numChannels,numTrials);
    
        for i=1:numChannels
            for j=1:numTrials
                pxx(:,i,j)=pwelch(series.estimated(:,i,j),window,overlap,freqRange,fs);
            end
        end

        if strcmp(plotType,'ind')
            [ax_diag,ax_offdiag,yLimits]=plot_ind(pxx,conn,h,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avg')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(pxx,conn,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avgerr')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(pxx,conn,all_lineprops{1},all_lineprops{colorNum});
        end
        
        colorNum=colorNum+1;
    end
    
    if isfield(series,'original_psd')
        if strcmp(plotType,'ind')
            [ax_diag,ax_offdiag,yLimits]=plot_ind(series.original_psd,conn,h,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avg')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(series.original_psd,conn,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avgerr')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(series.original_psd,conn,all_lineprops{1},all_lineprops{colorNum});
        end
        
        colorNum=colorNum+1;
    end
    
    if isfield(series,'estimated_psd')
        if strcmp(plotType,'ind')
            [ax_diag,ax_offdiag,yLimits]=plot_ind(series.estimated_psd,conn,h,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avg')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn]=plot_avg(series.estimated_psd,conn,all_lineprops{1},all_lineprops{colorNum});
        elseif strcmp(plotType,'avgerr')
            [ax_diag,ax_offdiag,yLimits,avgPSD,avgConn,stdPSD,stdConn]=plot_avgerr(series.estimated_psd,conn,all_lineprops{1},all_lineprops{colorNum});
        end
        
        colorNum=colorNum+1;
    end
    
    if isfield(series,'surrogate')
        bool_plotThreshold=true;
        [ax_diag,ax_offdiag,~]=plot_surrogate(series.surrogate,threshold,all_lineprops{colorNum});
        bool_plotThreshold=false;
    end
    
    if isempty(yLimits)
        yLimits=ylim;
    end
end

if bool_changeFreqRange
    set_figure(ax_diag,ax_offdiag,yLimits,freqLims);
else
    set_figure(ax_diag,ax_offdiag,yLimits);
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

%% Internal function to plot individual trials

    function [ax_diag,ax_offdiag,yLim]=plot_ind(y_diag,y_offdiag,h,lineprops_offdiag,lineprops_diag)
        
        if nargin == 3 || (isempty(lineprops_diag) && sempty(lineprops_offdiag))
            lineprops_offdiag='-b';
            lineprops_diag='-k';
        elseif nargin == 3 || isempty(lineprops_diag)
            lineprops_diag='-k';
        end
        
        if isempty(y_diag)
            y_diag=nan(numFrequencies,numChannels,numTrials);
        end
            
        for k=1:numTrials
            [ax_diag,ax_offdiag]=plotting(y_diag(:,:,k),y_offdiag(:,:,:,k),[],[],h(k,:),lineprops_offdiag,lineprops_diag);
        end

        if bool_plotTransferFunction
            yLim=[0 max(max(y_diag.^2))+2];
        elseif bool_plotSpectralMatrix
            yLim=[0 max(max(y_diag))+5];
        else
            yLim=[min(min(10*log10(y_diag)))-5,max(max(10*log10(y_diag)))+5];
        end
    end

%% Internal function to plot averages

    function [ax_diag,ax_offdiag,yLimits,avg_diag,avg_offdiag]=plot_avg(y_diag,y_offdiag,lineprops_offdiag,lineprops_diag)
        
        if nargin == 2 || (isempty(lineprops_diag) && sempty(lineprops_offdiag))
            lineprops_offdiag='-b';
            lineprops_diag='-k';
        elseif nargin == 3 || isempty(lineprops_diag)
            lineprops_diag='-k';
        end
        
        if isempty(y_diag)
            y_diag=nan(numFrequencies,numChannels,numTrials);
        end
            
        avg_diag=mean(y_diag,3);
        avg_offdiag=mean(y_offdiag,4);
        
        [ax_diag,ax_offdiag]=plotting(avg_diag,avg_offdiag,[],[],[],lineprops_offdiag,lineprops_diag);
        
        if bool_plotTransferFunction
            yLimits=[0 max(max(avg_diag.^2))+2];
        elseif bool_plotSpectralMatrix
            yLimits=[0 max(max(avg_diag))+5];
        else
            yLimits=[min(min(10*log10(avg_diag)))-5,max(max(10*log10(avg_diag)))+5];
        end
    end

%% Internal function to plot averages with shaded error bars

    function [ax_diag,ax_offdiag,yLimits,avg_diag,avg_offdiag,std_diag,std_offdiag]=plot_avgerr(y_diag,y_offdiag,lineprops_offdiag,lineprops_diag)
        
        if nargin == 2 || (isempty(lineprops_diag) && sempty(lineprops_offdiag))
            lineprops_offdiag='-b';
            lineprops_diag='-k';
        elseif nargin == 3 || isempty(lineprops_diag)
            lineprops_diag='-k';
        end
        
        if isempty(y_diag)
            y_diag=nan(numFrequencies,numChannels,numTrials);
        end
            
        avg_diag=mean(y_diag,3);
        std_diag=std(10*log10(y_diag),0,3);  
        avg_offdiag=mean(y_offdiag,4);
        std_offdiag=std(y_offdiag,0,4); 
        
        [ax_diag,ax_offdiag]=plotting(avg_diag,avg_offdiag,std_diag,std_offdiag,[],lineprops_offdiag,lineprops_diag);
        
        if bool_plotTransferFunction
            yLimits=[0 max(max(avg_diag.^2+std_diag))+2];
        elseif bool_plotSpectralMatrix
            yLimits=[0 max(max(avg_diag+std_diag))+5];
        else
            yLimits=[min(min(10*log10(avg_diag)-std_diag))-5,max(max(10*log10(avg_diag)+std_diag))+5];
        end
    end

%% Internal function to plot the significance threshold from surrogate analysis

    function [ax_diag,ax_offdiag,yLim]=plot_surrogate(surrogate,threshold,lineprops_offdiag)
        
        if nargin == 2 || isempty(lineprops_offdiag)
            lineprops_offdiag='--k';
        end
        
        numFrequencies=size(surrogate,3);
        numIterations=size(surrogate,4);
        
        significant_values=ones(numChannels,numChannels,numFrequencies);
        
        for k=1:numChannels
            for l=1:numChannels
                if k ~= l
                    for m=1:numFrequencies
                        orderedConnections=sort(squeeze(surrogate(k,l,m,:)));
                        significant_values(k,l,m)=orderedConnections(round(numIterations*(1-threshold)));
                    end
                end
            end
        end
        
        [ax_diag,ax_offdiag]=plotting(nan(numFrequencies,numChannels),significant_values,[],[],[],lineprops_offdiag);
        yLim=ylim;
    end

%% Internal plotting function to deal with a single trials worth of signals

    function [ax_diag,ax_offdiag]=plotting(pxx,conn,pxx_std,conn_std,h,lineprops_offdiag,lineprops_diag)
        
        if nargin == 4 || (nargin == 7 && isempty(h) && (~isempty(pxx_std) || ~isempty(conn_std)))
            bool_plotErrorBars=true;
        else
            bool_plotErrorBars=false;
        end
        
        if nargin < 6
            lineprops_diag='-k';
            lineprops_offdiag='-b';
        elseif nargin < 7
            lineprops_diag='-k';            
        end
        
        if isempty(conn)
            conn=nan(numChannels,numChannels,numFrequencies);
            
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
                        shadedErrorBar(freqRange,10*log10(pxx(:,k)),pxx_std(:,k),'lineProps',lineprops_diag);
                    elseif bool_plotThreshold
                        ax=gca;
                        
                        if isempty(ax.Children)
                            axis off;
                        end
                    elseif bool_plotTransferFunction
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,pxx(:,k).^2,'r'); % Not actually pxx, is actually H 
                        else
                            plot(freqRange,pxx(:,k).^2,lineprops_diag); 
                        end
                    elseif bool_plotSpectralMatrix
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,pxx(:,k),'r'); % Not actually pxx, is actually H 
                        else
                            plot(freqRange,pxx(:,k),lineprops_diag); 
                        end
                    else
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,10*log10(pxx(:,k)),'r'); 
                        else
                            plot(freqRange,10*log10(pxx(:,k)),lineprops_diag); 
                        end
                    end
                    
                    title(labels{k}); hold on;
                elseif k ~= l 
                    ax_offdiag(currSubPlot)=subplot(numChannels,numChannels,currSubPlot);
                    
                    if bool_plotErrorBars
                        shadedErrorBar(freqRange,conn(k,l,:),conn_std(k,l,:),'lineProps',lineprops_offdiag);
                    else
                        if bool_showRejectedNull && (h(k) || h(l))
                            plot(freqRange,squeeze(conn(k,l,:)),'r'); 
                        elseif bool_plotThreshold
                            plot(freqRange,squeeze(conn(k,l,:)),':k','LineWidth',3);
                        else
                            plot(freqRange,squeeze(conn(k,l,:)),lineprops_offdiag); 
                        end
                    end
                    
                    hold on;
                    title(sprintf('%s %c %s',labels{l},8594,labels{k}));
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
        ylim([yLimits(1),yLimits(2)]);
        subplot(numChannels,numChannels,numChannels); 

        linkaxes(ax_offdiag); 
        xlim([xLimits(1), xLimits(2)]); 
        ylim([0 1]);

    end
end