%% pipeline.m
%
%  Combines both the processing_pipeline.m and postprocessing_pipeline.m files, so that
%  all processing is done at the same time.
%
%  See also: dtf, mvar, load_data, load_channels_labels_conditions, test_model,
%   simplify_filename, surrogate_analysis, filter_serial_correlation, test_model
%

% Make sure that the parallel pool is initialized

CCC;
dtf_startup;

%% Definitions

meta=get_path_variables(6,1,1,'_NEW_PATIENT','RC+S recorded the Vim signals');

meta.settings.cues_only=true;
meta.settings.extrap_method='';
meta.settings.alpha=0.05;
meta.settings.freqForAnalysis=4:0.5:100;

%% Load Data

[file,file_postacq]=create_file_path(meta);

meta=load_variables(meta);

if ~meta.settings.cues_only
    datastorage_postacq=load(file_postacq,'datastorage_postacq');
    datastorage_postacq.postacq_type=meta.vars.postacq_type;
else
    datastorage_postacq=struct;
end

datastorage=load(file,'datastorage');
datastorage=datastorage.datastorage;
fs_init=extract_sampling_frequency(datastorage);

datastorage_dtf=struct;

fs_new=200;
filtering=set_filtering_parameters(fs_init,fs_new);

[datastorage_dtf.x_all,datastorage_dtf.instruct]=load_data(datastorage,meta,filtering,datastorage_postacq,[]);
fs=fs_new;

numChannels=size(datastorage_dtf.x_all,2);
numConditions=length(meta.vars.conditions);

config=struct;

config.mvar=struct(...
    'orderSelection','diff2',...
    'crit','bic',...
    'method','arfit',...
    'orderRange',1:20,...   % Used to find the optimal model order
    'modelOrder',[],...     % This is the optimal model order
    'fs',fs,...
    'epsilon',0.001,...
    'freqRange',meta.settings.freqForAnalysis,...
    'output',0,...          % Since I'm using the optimal model order, don't output the search in mvar
    'logLikelihoodMethod',2);   % 1 == Matlab, 2 == Ding, 3 == Awareness Paper (see calculate_bic.m)
config.plot=struct(...
    'hFig',[],...
    'seriesType',1,...
    'plotType','avgerr',...
    'freqLims',meta.settings.freqForAnalysis,...
    'fs',fs);
config.surrogate=struct(...
    'method','sample',...
    'numSamples',1,...
    'iterations',1000);

%% Main Loop for initial processing

fprintf('First pass beginning...\n');

for j=1:numConditions
    currCond=meta.vars.cond_labels{j};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    datastorage_dtf.x.(currCond)=load_data(datastorage,meta,filtering,datastorage_postacq,meta.vars.conditions(j));
    
    numTrials=size(datastorage_dtf.x.(currCond),3);
    numChannels=length(meta.vars.channels);

    datastorage_dtf.test.h.(currCond)=nan(numTrials,numChannels);
    datastorage_dtf.test.pVal.(currCond)=nan(numTrials,numChannels);
    datastorage_dtf.test.pass.(currCond)=nan(numTrials,1);

    fprintf('Calculating optimal order...');

    [optimal_order,datastorage_dtf.criterion.(currCond)]=calculate_optimal_model_order(datastorage_dtf.x.(currCond),config.mvar);
    config.mvar.modelOrder=optimal_order;
    
    fprintf('Optimal Model Order = %d\n',optimal_order);

    for i=1:numTrials
        fprintf('%d/%d\n',i,numTrials);
        
        %% Calculate MVAR model
        [datastorage_dtf.mvar.ar.(currCond)(i).mdl, datastorage_dtf.mvar.res.(currCond)(i).E] = mvar(squeeze(datastorage_dtf.x.(currCond)(:,:,i)), config.mvar);
        
        %% Test Whiteness
        [datastorage_dtf.test.pass.(currCond)(i), datastorage_dtf.test.h.(currCond)(i,:), datastorage_dtf.test.pVal.(currCond)(i,:)]=...
            test_model(datastorage_dtf.mvar.res.(currCond)(i).E,size(datastorage_dtf.mvar.res.(currCond)(i).E,1));
    end
end

fprintf('First pass completed.\n');

%% Decorrelation

fprintf('Decorrelating time series channels...\n');
[datastorage_dtf.x,filt_values]=filter_serial_correlation(datastorage_dtf.x,datastorage_dtf.mvar.res,datastorage_dtf.test.h,config.mvar);
fprintf('Decorrelation completed.\n');

%% Main Loop for decorrelated data

fprintf('Decorrelated data connectivity calculations beginning...\n');

for i=1:numConditions
    currCond=meta.vars.cond_labels{i};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    numTrials=size(datastorage_dtf.x.(currCond),3);
    
    datastorage_dtf.test.h.(currCond)=nan(numTrials,numChannels);
    datastorage_dtf.test.pVal.(currCond)=nan(numTrials,numChannels);
    datastorage_dtf.test.pass.(currCond)=ones(numTrials,1);

    datastorage_dtf.gamma.(currCond)=zeros(numChannels,numChannels,length(meta.settings.freqForAnalysis),numTrials);

    [optimal_order,datastorage_dtf.criterion.(currCond)]=calculate_optimal_model_order(datastorage_dtf.x.(currCond), config.mvar);
    config.mvar.modelOrder=optimal_order;

    fprintf('Optimal Model Order = %d\n',optimal_order);

    for j=1:numTrials
        if filt_values.(currCond)(j).decorrelated
            fprintf('%d/%d\n',j,numTrials);

            %% Calculate MVAR model
            [datastorage_dtf.mvar.ar.(currCond)(j).mdl,datastorage_dtf.mvar.res.(currCond)(j).E]=mvar(datastorage_dtf.x.(currCond)(:,:,j),config.mvar);

            %% Test whiteness
            [datastorage_dtf.test.pass.(currCond)(j),datastorage_dtf.test.h.(currCond)(j,:),datastorage_dtf.test.pVal.(currCond)(j,:)]=...
                test_model(datastorage_dtf.mvar.res.(currCond)(j).E,size(datastorage_dtf.mvar.res.(currCond)(j).E,1));

            %% Calculate DTF Connectivity
            datastorage_dtf.gamma.(currCond)(:,:,:,j)=directedTransferFunction(datastorage_dtf.mvar.ar.(currCond)(j).mdl,...
                meta.settings.freqForAnalysis,fs);
        else
            warning('Trial %d of %d was not properly decorrelated',j,numTrials)
        end
    end
    
    print_whiteness(datastorage_dtf.test.h.(currCond),datastorage_dtf.test.pVal.(currCond),meta.vars.labels);
end

fprintf('Decorrelated data connectivity calculations completed.\n');

%% Surrogate Analysis on Decorrelated Data

fprintf('Surrogate analysis of decorrelated data beginning...\n');
[datastorage_dtf.surrogate.data,datastorage_dtf.surrogate.distribution]=surrogate_analysis(datastorage_dtf.x,fs,meta.settings.freqForAnalysis,config.mvar,config.surrogate);
fprintf('Surrogate analyis of decorrelated data completed.\n');

%% Save all relevant variables

newFile=simplify_filename(meta);

meta.vars.contactNames=get_structure_names(meta.path.patID);

meta.date_pipeline_was_run=datetime;

save(newFile,'config','datastorage_dtf','filtering','meta');










