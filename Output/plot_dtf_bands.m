%% plot_dtf_bands
%
%  Extract the DTF values from the gamma struct, and plot bar graphs comparing the values
%  in different frequency bands for all connections. Assumes the file has already been
%  opened
%

%% Load Variables 

numConditions=length(meta.vars.cond_labels);
freqBands=[4,8;8,12;12,20;20,30;30,45;45,70;70,100];
freqBandIndices=convert_from_hertz_to_indices(freqBands,meta.settings.freqForAnalysis);
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Low Gamma','Mid Gamma','High Gamma'};

significance=calculate_significance_from_surrogate(datastorage_dtf.surrogate.data,meta.settings.alpha,'invariant');

%% Plot one frequency band per figure, with all conditions in different colors, and reject non-significant values

config_plot=struct;
config_plot.hideUselessConnections=true;
config_plot.usefulConnections=returnUsefulConnections(freqBandIndices,datastorage_dtf.gamma,significance);

for i=1:size(freqBands,1)
    config_plot.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBandIndices{i},datastorage_dtf.gamma,meta.vars.contactNames,significance,config_plot);
end

%% Plot one figure per frequency band, with each conditions plotted in the same bar graph

config_plot=struct;

for i=1:length(freqBands)
    config_plot.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBands{i},gamma_filt,labels,[],config_plot);
end

%% Plot one figure per condition, with each frequency band plotted in the same bar graph

config_plot=struct;

for i=1:numConditions
    currCond=cond_labels{i};
    
%     config.figTitle=sprintf('Band Connectivity During %s',currCond);
    config_plot.title=sprintf('Band Connectivity During %s',currCond);
    config_plot.legend=freqBandLabel;
    plot_bar_with_error(freqBands,gamma_filt.(currCond),labels,[],config_plot);
end

%% Plot the bars in the same format as the [d x d] plot of all connectivity values

config_plot=struct;

numChannels=size(gamma_filt.Rest,1);
colors=linspecer(length(freqBandLabel));

for k=1:numConditions
    currCond=cond_labels{k};
    config_plot.figHandle=figure;
    config_plot.figTitle=sprintf('Band Connectivity During %s',currCond);

    for i=1:numChannels
        for j=1:numChannels
            currSubPlot=(i-1)*numChannels+j;
            
            if i~=j
                label=cellstr(sprintf('%s %c %s',labels{i},8594,labels{j}));
                config_plot.axHandle=subplot(numChannels,numChannels,currSubPlot);
                plot_bar_with_error(freqBands,gamma_filt.(currCond)(i,j,:,:),label,[],config_plot)
            else
                subplot(numChannels,numChannels,currSubPlot);
                axis off; hold on;
                
                for l=1:length(freqBandLabel)
                    plot(nan,nan,'LineWidth',10,'Color',colors(l,:)); 
                end
                
                legend(freqBandLabel,'Location','East');
            end
        end
    end
end

%% Plot individual figures with one frequency band per figure

% for i=1:numConditions
%     currCond=cond_labels{i};
%     
%     for j=1:length(freqBands)
%         config.figTitle=sprintf('%s Band Connectivity During %s',freqBandLabel{j},currCond);
%         plot_bar_with_error(freqBands{j},gamma.(currCond),labels,[],config);
%     end
% end

