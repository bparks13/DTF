function emptyConfig=create_empty_config(type)
%% emptyConfig=create_empty_config(type)
%
%  Given a specific type (i.e. 'crit' or 'plot'), returns an empty struct that contains
%  all of the possible fields to be populated. 
%
%   Inputs:
%    - type: String defining which type of struct to create, dependent on which function
%       the config is going to.
%         'mvar': Defines an empty config that is used by mvar.m
%         'plot_connectivity': Defines an empty config that is used by plot_connectivity.m
%         'filter_serial_correlation': Defines an empty config that is used by
%           filter_serial_correlation.m 
%         'surrogate_analysis': Defines an empty config that is used by
%           surrogate_analysis.m 
%         'test_model': Defines an empty config that is used by test_model.m
%         'load_channels_labels_conditions': Defines an empty config that is used by
%           load_channels_labels_conditions.m 
%         'load_data': Defines an empty config that is used by load_data.m
%         'plot_criterion': Defines an empty config that is used by plot_criterion.m
%         'plot_psd': Defines an empty config that is used by plot_psd.m
%
%   Outputs:
%    - emptyConfig: Empty struct defining all of the possible fields that can be populated
%       by a specific config
%
%  See also: mvar, plot_connectivity, filter_serial_correlation, surrogate_analysis,
%   test_model, load_channels_labels_conditions, load_data, plot_criterion, plot_psd
%

if strcmp(type,'mvar') % mvar
    emptyConfig=struct(...
        'orderRange',[],...
        'crit',[],...
        'output',[],...
        'method',[],...
        'orderSelection',[],...
        'modelOrder',[],...
        'epsilon',[],...
        'fs',[],...
        'freqRange',[],...
        'logLikelihoodMethod',[],...
        'simulated',struct('a',[],'C',[]),...
        'spectral_range',[]);
elseif strcmp(type,'plot_connectivity') % plot_connectivity
    emptyConfig=struct(...
        'seriesType',[],...
        'fs',[],...
        'hFig',[],...
        'figTitle',[],...
        'plotType',[],...
        'h',[],...
        'freqLims',[],...
        'yLims',[],...
        'surr_params',struct(...
            'threshold',[],...
            'highlightSignifiance',[],...
            'binning',[]));
elseif strcmp(type,'filter_serial_correlation') % filter_serial_correlation
    emptyConfig=struct(...
        'maxIterations',[],...
        'modelOrders',[]);
elseif strcmp(type,'surrogate_analysis')  % surrogate_analysis
    emptyConfig=struct(...
        'cond',[],...
        'method',[],...
        'iterations',[],...
        'numSamples',[]);
elseif strcmp(type,'test_model')  % test_model
    emptyConfig=struct(...
        'overrideLags',[],...
        'lags',[],...
        'changeDOF',[],...
        'numParameters',[]);
elseif strcmp(type,'load_channels_labels_conditions')  % load_variables
    emptyConfig=struct(...
        'default',[],...
        'preset',[],...
        'custom',[]);
elseif strcmp(type,'load_data') % load_data 
    emptyConfig=struct(...
        'hpf',[],...
        'comb',[],...
        'lpf',[],...
        'notch',[],...
        'ma',[],...
        'normalize',[],...
        'downsample',[]);
elseif strcmp(type,'plot_criterion') % plot_criterion
    emptyConfig=struct(...
        'hFig',[],...
        'hAx',[],...
        'epsilon',[],...
        'modelOrder',[],...
        'figTitle',[],...
        'average',[]);
elseif strcmp(type,'plot_psd') % plot_psd
    emptyConfig=struct(...
        'inputType',[],...
        'hFig',[]);
elseif strcmp(type,'plot_bar_with_error') % plot_bar_with_error
    emptyConfig=struct(...
        'yLim',[],...
        'title',[],...
        'figTitle',[],...
        'legend',[],...
        'figHandle',[],...
        'axHandle',[],...
        'hideUselessConnections',[],...
        'usefulConnections',[]);
end

end