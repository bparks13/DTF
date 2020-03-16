function plot_bar_with_error(freqBand,gamma,labels,significance,config)
%% plot_bar_with_error(freqBand,gamma,labels,significance,config)
%
%  Plots a bar graph, with error bars, of the average connectivity values over the
%  frequency band(s) given. This function can plot a single condition and a single
%  frequency band, a single condition and multiple frequency bands, or a single frequency
%  band and multiple conditions. It cannot plot multiple frequency bands and multiple
%  conditions
%
%   Inputs:
%    - freqBand: Indices pertaining to the frequencies in gamma to average over. Can
%       either be a vector containing one single frequency band, or a cell array
%       containing multiple frequency bands. If defined as a cell array, gamma cannot be a
%       struct
%    - gamma: Matrix containing all of the connectivity values calculated for this
%       condition. Size is [d x d x f x t], where d is the number of channels, f is the
%       frequencies, and t is the number of trials. A single combination of channels can
%       be given and will be plotted. This can be combined with a subplot axis handle for
%       plotting multiple channel combinations in a figure with multiple subplots.
%       Can also be a struct containing multiple conditions, which follows the same format
%       as the matrix above, just with multiple fields pertaining to each condition. If
%       defined as a struct, freqBand cannot be a cell array
%    - labels: Cell array containing the channel labels. If only one combination of
%       channels is given, this should only contain one entry, which is the already
%       formatted string containing the combination (i.e. Ch1 ? Ch2)
%    - significance: Optional input - Matrix defining what connectivity value is
%       considered to be significant. Depending on the type specified, the size of
%       significance can change. For 'invariant', the size is [c x c], where c is the
%       number of channels. For 'dependent', size is [c x c x f], where f is the number of
%       frequencies. Significance type is automatically detected from the matrix size
%    - config: Optional struct containing additional parameters
%    -- yLim: Limits of the y-axis. Default is [0 1]
%    -- title: String containing the title of the axis
%    -- figTitle: String containing the title of the figure window
%    -- legend: Given multiple conditions/frequency bands on the same plot, the legend is
%           used to name the components of the figure
%    -- figHandle: Figure handle to an existing figure to plot on
%    -- axHandle: Axis handle to an existing subplot to plot on
%    -- hideUselessConnections: Boolean value to denote whether or not to plot all
%       possible connections (false, default), or to hide the connections that never have
%       a significant relation in any frequency (true). The significance variable must be
%       defined for this to have any meaning
%    -- usefulConnections: Matrix of channel combinations that are useful.
%
%   Outputs:
%    Figure containing the bar plot, with error bars, and the channel labels on the
%     x-axis. The error bars are 2 times the standard error
%
%  See also: plot_connectivity.m, returnUsefulConnections.m
%

yLim=[0 1];
axTitle='Connectivity Values';
figTitle='';
bool_showLegend=false;
bool_multipleBands=false;
bool_multipleConds=false;
bool_singleChannel=false;
bool_newFigure=true;
bool_newAxis=true;
bool_useSignificance=false;
bool_hideUselessConnections=false;
sigType='';
plotLegend={};

if iscell(freqBand)
    bool_multipleBands=true;
end

if isstruct(gamma)
    bool_multipleConds=true;
    condLabels=fieldnames(gamma);
    bool_showLegend=true;
    config.legend=condLabels;
    
    numChannels=size(gamma.(condLabels{1}),1);
    numTrials=size(gamma.(condLabels{1}),4);
else
    numChannels=size(gamma,1);
    numTrials=size(gamma,4);
end

if bool_multipleBands && bool_multipleConds
    error('Cannot plot multiple frequency bands and multiple conditions simultaneously');
end

if nargin >= 4 && ~isempty(significance)
    bool_useSignificance=true;
    
    if size(significance.(condLabels{1}),3) > 1
        sigType='dependent';
    else
        sigType='invariant';
    end
    
end

if nargin == 5
    if isstruct(config)
        if isfield(config,'yLim') && ~isempty(config.yLim)
            yLim=config.yLim;
        end
        
        if isfield(config,'title') && ~isempty(config.title)
            axTitle=config.title;
        end
        
        if isfield(config,'figTitle') && ~isempty(config.figTitle)
            figTitle=config.figTitle;
        end
        
        if isfield(config,'legend') && ~isempty(config.legend)
            bool_showLegend=true;
            plotLegend=config.legend;
        end
        
        if isfield(config,'figHandle') && ~isempty(config.figHandle)
            bool_newFigure=false;
        end
        
        if isfield(config,'axHandle') && ~isempty(config.axHandle)
            bool_newAxis=false;
        end
        
        if isfield(config,'hideUselessConnections') && ~isempty(config.hideUselessConnections)
            bool_hideUselessConnections=config.hideUselessConnections;
        end
        
        if isfield(config,'usefulConnections') && ~isempty(config.usefulConnections)
            bool_hideUselessConnections=true;
        end
    end
end

% Extract individual band

if numChannels~=1
    if bool_hideUselessConnections
        combinations=config.usefulConnections;
    else
        combinations=[nchoosek(1:numChannels,2);nchoosek(numChannels:-1:1,2)];
    end
    
    numCombinations=size(combinations,1);
else
    combinations=[1,1];
    numCombinations=size(combinations,1);
    bool_singleChannel=true;
end

if bool_newFigure
    figure('Name',figTitle);
else
    figure(config.figHandle);
    tmp_fig=gcf;
    
    if isempty(tmp_fig.Name)
        tmp_fig.Name=figTitle;
    end
end

if bool_multipleBands
    numBands=length(freqBand);
    avgGamma_band=nan(numCombinations,numBands);
    stdGamma_band=nan(numCombinations,numBands);
    
    for l=1:numBands
        gamma_band=nan(numCombinations,numTrials);

        for i=1:numCombinations
            for j=1:numTrials
                gamma_band(i,j)=mean(gamma(combinations(i,1),combinations(i,2),freqBand{l},j));
            end
        end
        
        if bool_useSignificance
            if strcmp(sigType,'invariant')
                error('Multiple bands is not configured yet');
            else
                error('''invariant'' is the only type of significance currently programmed');
            end
        else
            avgGamma_band(:,l)=mean(gamma_band,2);
%             stdGamma_band(:,l)=std(gamma_band,[],2);
            stdGamma_band(:,l)=2*ste(gamma_band);
        end
    end
    
    internal_plot(avgGamma_band,stdGamma_band);
elseif bool_multipleConds
    numConds=length(condLabels);
    avgGamma_band=nan(numCombinations,numConds);
    stdGamma_band=nan(numCombinations,numConds);
    
    for l=1:numConds
        numTrials=size(gamma.(condLabels{l}),4);
        gamma_band=nan(numCombinations,numTrials);
            
        for i=1:numCombinations
            for j=1:numTrials
                gamma_band(i,j)=mean(gamma.(condLabels{l})(combinations(i,1),combinations(i,2),freqBand,j));
            end
        end
        
        if bool_useSignificance
            if strcmp(sigType,'invariant')
                for i=1:numCombinations
                    tmp_mean=mean(squeeze(gamma_band(i,:)));

                    if tmp_mean > significance.(condLabels{l})(combinations(i,1),combinations(i,2))
                        avgGamma_band(i,l)=tmp_mean;
%                         stdGamma_band(i,l)=std(gamma_band(i,:),[],2);
                        stdGamma_band(i,l)=2*ste(gamma_band(i,:));
                    end
                end
            else
                error('''invariant'' is the only type of significance currently programmed');
            end
        else
            avgGamma_band(:,l)=mean(gamma_band,2);
%             stdGamma_band(:,l)=std(gamma_band,[],2);
            stdGamma_band(:,l)=2*ste(gamma_band);
        end
    end
    
    internal_plot(avgGamma_band,stdGamma_band);
else
    gamma_band=nan(numCombinations,numTrials);

    for i=1:numCombinations
        for j=1:numTrials
            gamma_band(i,j)=mean(gamma(combinations(i,1),combinations(i,2),freqBand,j));
        end
    end

    avgGamma_band=mean(gamma_band,2);
%     stdGamma_band=std(gamma_band,[],2);
    stdGamma_band=2*ste(gamma_band);

    internal_plot(avgGamma_band,stdGamma_band);
end

%%
    function internal_plot(avgGamma_band,stdGamma_band)
        % Plot bar graph

        if bool_newAxis
            if bool_singleChannel
                pos=1:length(avgGamma_band);
                
                for k=1:length(pos)
                    bar(pos(k),avgGamma_band(k)); hold on;
                end
            else
                pos=1:numCombinations;
                bar(pos,avgGamma_band);
            end
        else
            if bool_singleChannel
                pos=1:length(avgGamma_band);
                
                for k=1:length(pos)   
                    bar(config.axHandle,pos(k),avgGamma_band(k)); hold on;  
                end
            else
                pos=1:numCombinations;
                bar(config.axHandle,pos,avgGamma_band);
            end
        end
        
        hold on;
        internal_errorbar(avgGamma_band,stdGamma_band,pos);
        
        ylim(yLim);

        title(axTitle); ylabel('Normalized Connectivity [DTF]');

        tmp_labels=cell(numCombinations,1);

        for k=1:numCombinations
            tmp_labels{k}=sprintf('%s %c %s',labels{combinations(k,2)},8594,labels{combinations(k,1)});
        end

        xticks(pos);
        ax=gca;
        
        if bool_singleChannel
            ax.XTick=[];
            title(labels);
        else
            ax.XTickLabel=tmp_labels;
        end
        
        ax.XTickLabelRotation=45;

        if bool_showLegend
            legend(plotLegend);
        end

        drawnow;
    end

%%
    function internal_errorbar(avgGamma_band,stdGamma_band,pos)
        
        if bool_multipleConds
            groupwidth = min(0.8, numConds/(numConds + 1.5));
            colors=linspecer(numConds);
            ax=gca;
            
            for k = 1:numConds
                if bool_singleChannel
                    x=pos;
                else
                    x = pos - groupwidth/2 + (2*k-1) * groupwidth / (2*numConds);
                end
                
                if bool_newAxis
                    errorbar(x,avgGamma_band(:,k),stdGamma_band(:,k),'.k');
                else
                    errorbar(config.axHandle,x,avgGamma_band(:,k),stdGamma_band(:,k),'.k');                    
                end
                
                ax.Children(end-k+1).FaceColor=colors(k,:);
            end
        elseif bool_multipleBands
            groupwidth = min(0.8, numBands/(numBands + 1.5));
            colors=linspecer(numBands);
            ax=gca;
            
            if bool_singleChannel
                x=pos;
                
                for k = 1:numBands
                    if bool_newAxis
                        errorbar(x(k),avgGamma_band(k),stdGamma_band(k),'.k');
                    else
                        errorbar(config.axHandle,x(k),avgGamma_band(k),stdGamma_band(k),'.k');                    
                    end

                    ax.Children(end-k+1).FaceColor=colors(k,:);
                end
            else
                for k = 1:numBands
                    x = pos - groupwidth/2 + (2*k-1) * groupwidth / (2*numBands);

                    if bool_newAxis
                        errorbar(x,avgGamma_band(:,k),stdGamma_band(:,k),'.k');
                    else
                        errorbar(config.axHandle,x,avgGamma_band(:,k),stdGamma_band(:,k),'.k');                    
                    end

                    ax.Children(end-k+1).FaceColor=colors(k,:);
                end
            end
        else
            if bool_newAxis
                errorbar(pos,avgGamma_band,avgGamma_band-stdGamma_band,avgGamma_band-stdGamma_band,...
                    'k','LineStyle','None');
            else
                errorbar(config.axHandle,pos,avgGamma_band,avgGamma_band-stdGamma_band,avgGamma_band-stdGamma_band,...
                    'k','LineStyle','None');
            end
        end
    end

end