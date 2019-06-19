%% plot_dtf_bands
%
%  Extract the DTF values from the gamma struct, and plot bar graphs comparing the values
%  in different frequency bands for all connections
%

CCC;

%% Load data

FILE='ET_CL_001__2017_05_17__run12__PSD__Z_SCORE.mat';

load(fullfile(get_root_path(),'Files',FILE));

numConditions=length(cond_labels);
freqBands={4:8,8:12,12:20,20:30,30:100};
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Gamma'};

config=struct;

%% Plot individual figures with one frequency band per figure

% for i=1:numConditions
%     currCond=cond_labels{i};
%     
%     for j=1:length(freqBands)
%         config.figTitle=sprintf('%s Band Connectivity During %s',freqBandLabel{j},currCond);
%         plot_bar_with_error(freqBands{j},gamma.(currCond),labels,config);
%     end
% end

%% Plot one figure per condition, with each frequency band plotted in the same bar graph

for i=1:numConditions
    currCond=cond_labels{i};
    
    config.figTitle=sprintf('Band Connectivity During %s',currCond);
    config.legend=freqBandLabel;
    plot_bar_with_error(freqBands,gamma.(currCond),labels,config);
end

%% Plot one figure per frequency band, with each conditions plotted in the same bar graph

for i=1:length(freqBands)
    config.figTitle=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBands{i},gamma,labels,config);
end

