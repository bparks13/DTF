%#ok<*UNRCH>

CCC;

print_figures=true;
PRINT_FOLDER=fullfile(get_root_path,'Pictures','InitialResults');

freqBands={4:8,8:12,12:20,20:30,30:100};
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Gamma'};

%% Load S01 (ET_CL_02)

load(fullfile(get_root_path,'Files','S01_D01_R01.mat'));

significance=calculate_significance_from_surrogate(surrogate_filt,0.01,'invariant');

%% Plot overall connectivity for Rest only, and resize

config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity (Decorrelated)',PATIENT_ID,RECORDING_DATE,RUN_ID,'Rest');
series=struct('original',x_filt.Rest,'surrogate',surrogate_filt.Rest);
config_plot.surr_params.highlightSignificance=true;
plot_connectivity(gamma_filt.Rest,series,freqRange,labels,config_plot);
hFig=gcf; hFig.Position=[hFig.Position(1:2)-300,1200,700];

if print_figures
    PrintFigure(fullfile(PRINT_FOLDER,'S01_CompleteConnectivityMatrix'));
end

%% Plot bar plots, with only significant results

config=struct;

for i=1:length(freqBands)
    config.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBands{i},gamma_filt,labels,significance,config);
    hFig=gcf; hFig.Position=[hFig.Position(1:2),950,350];
    
    if print_figures
        PrintFigure(fullfile(PRINT_FOLDER,sprintf('S01_BarPlot_%sBand',freqBandLabel{i})));
    end
end

%% Close the previous subject

CCC;

print_figures=true;
PRINT_FOLDER=fullfile(get_root_path,'Pictures','InitialResults');

freqBands={4:8,8:12,12:20,20:30,30:100};
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Gamma'};

%% Load S02 (ET_CL_04)

load(fullfile(get_root_path,'Files','S02_D01_R01.mat'));

significance=calculate_significance_from_surrogate(surrogate_filt,0.01,'invariant');

%% Plot bar plots, with only significant results

config=struct;

for i=1:length(freqBands)
    config.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBands{i},gamma_filt,labels,significance,config);
    hFig=gcf; hFig.Position=[hFig.Position(1:2),950,350];
    
    if print_figures
        PrintFigure(fullfile(PRINT_FOLDER,sprintf('S02_BarPlot_%sBand',freqBandLabel{i})));
    end
end

%% Close the previous subject

CCC;

print_figures=true;
PRINT_FOLDER=fullfile(get_root_path,'Pictures','InitialResults');

freqBands={4:8,8:12,12:20,20:30,30:100};
freqBandLabel={'Theta','Alpha','Low Beta','High Beta','Gamma'};

%% Load S03 (ET_OR_018)

load(fullfile(get_root_path,'Files','S03_D01_R01.mat'));

significance=calculate_significance_from_surrogate(surrogate_filt,0.01,'invariant');

%% Plot bar plots, with only significant results

config=struct;

for i=1:length(freqBands)
    config.title=sprintf('%s Band Connectivity',freqBandLabel{i});
    plot_bar_with_error(freqBands{i},gamma_filt,labels,significance,config);
    hFig=gcf; hFig.Position=[hFig.Position(1:2),950,350];
    
    if print_figures
        PrintFigure(fullfile(PRINT_FOLDER,sprintf('S03_BarPlot_%sBand',freqBandLabel{i})));
    end
end

