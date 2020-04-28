% Monday Morning Meeting - 04/13/2020?
%#ok<*UNRCH>

CCC;

print_figures=true;

%% Load file

cd(get_root_path);
load('Files\S02_D01_R01');

%% Plot all connectivity, resize, and print the figure

currCond=cond_labels{1};
config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity (Decorrelated)',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
series=struct('original',x_filt.(currCond),'surrogate',surrogate_filt.(currCond));
config_plot.surr_params.highlightSignificance=true;
config_plot.surr_params.threshold=alpha;
plot_connectivity(gamma_filt.(currCond),series,freqForAnalysis,contactNames,config_plot);

hFig=gcf;
hFig.Position=[hFig.Position(1:2)-300,900,650];

if print_figures
    PrintFigure('Pictures\MondayMeeting\S02_D01_R01_Rest_Connectivity');
end

%% Setup the barplots, plot them, resize, and print

numConditions=length(cond_labels);
freqBands=[4,8;8,12;12,20;20,30;30,45;45,70;70,100];
freqBandIndices=convert_from_hertz_to_indices(freqBands,freqForAnalysis);
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Low Gamma','Mid Gamma','High Gamma'};
significance=calculate_significance_from_surrogate(surrogate_filt,alpha,'invariant');

config=struct;
config.hideUselessConnections=true;
config.usefulConnections=returnUsefulConnections(freqBandIndices,gamma_filt,significance);

for i=1:size(freqBands,1)
    config.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBandIndices{i},gamma_filt,contactNames,significance,config);
    ylim([0 0.54])
    
    hFig=gcf; hFig.Position=[hFig.Position(1:2),900,400];
    
    if print_figures
        PrintFigure(sprintf('Pictures\\MondayMeeting\\S02_D01_R01_BarPlot_%s',freqBandLabel{i}))
    end
end


