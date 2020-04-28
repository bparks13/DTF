%% plot_dtf_bands
%
%  Extract the DTF values from the gamma struct, and plot bar graphs comparing the values
%  in different frequency bands for all connections
%

CCC;
dtf_startup;

%% Load data

% FILE='ET_CL_001__2017_05_17__run12__PSD__Z_SCORE.mat';
% 
% load(fullfile(get_root_path(),'Files',FILE));

numConditions=length(cond_labels);
% freqBands={4:8,8:12,12:20,20:30,30:100}; % Define the frequency values
freqBands=[4,8;8,12;12,20;20,30;30,45;45,70;70,100];
freqBandIndices=convert_from_hertz_to_indices(freqBands,freqForAnalysis);
% freqBandIndices={1:9,9:17,17:33,33:53,53:193};
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Low Gamma','Mid Gamma','High Gamma'};

% contactNames=get_structure_names(subjID);

%% Use the surrogate values to create a significance threshold

significance=calculate_significance_from_surrogate(surrogate_filt,alpha,'invariant');
% significance=calculate_significance_from_surrogate(surrogate_filt,0.01,'dependent');

%% Plot one frequency band per figure, with all conditions in different colors, and reject non-significant values

config=struct;
config.hideUselessConnections=true;
config.usefulConnections=returnUsefulConnections(freqBandIndices,gamma_filt,significance);

for i=1:size(freqBands,1)
    config.title=sprintf('%s Band Connectivity',freqBandLabel{i});
%     freqIndStart=find(freqForAnalysis==freqBands{i}(1),1);
%     freqIndEnd=find(freqForAnalysis==freqBands{i}(end),1);
    plot_bar_with_error(freqBandIndices{i},gamma_filt,contactNames,significance,config);
end

%% Plot one figure per frequency band, with each conditions plotted in the same bar graph

config=struct;

for i=1:length(freqBands)
    config.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBands{i},gamma_filt,labels,[],config);
end

%% Plot one figure per condition, with each frequency band plotted in the same bar graph

config=struct;

for i=1:numConditions
    currCond=cond_labels{i};
    
%     config.figTitle=sprintf('Band Connectivity During %s',currCond);
    config.title=sprintf('Band Connectivity During %s',currCond);
    config.legend=freqBandLabel;
    plot_bar_with_error(freqBands,gamma_filt.(currCond),labels,[],config);
end

%% Plot the bars in the same format as the [d x d] plot of all connectivity values

config=struct;

numChannels=size(gamma_filt.Rest,1);
colors=linspecer(length(freqBandLabel));

for k=1:numConditions
    currCond=cond_labels{k};
    config.figHandle=figure;
    config.figTitle=sprintf('Band Connectivity During %s',currCond);

    for i=1:numChannels
        for j=1:numChannels
            currSubPlot=(i-1)*numChannels+j;
            
            if i~=j
                label=cellstr(sprintf('%s %c %s',labels{i},8594,labels{j}));
                config.axHandle=subplot(numChannels,numChannels,currSubPlot);
                plot_bar_with_error(freqBands,gamma_filt.(currCond)(i,j,:,:),label,[],config)
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

