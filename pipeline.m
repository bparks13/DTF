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

% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
% PATIENT_ID='ET_CL_002';
% RECORDING_DATE='2018_02_01';
% RUN_ID='run12';
% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
% PATIENT_ID='ET_CL_004';
% RECORDING_DATE='2018_06_20';
% RUN_ID='run5';
% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET\\OR\\with_DBS';
% PATIENT_ID='ET_OR_STIM_018';
% RECORDING_DATE='2018_11_28';
% RUN_ID='run12';
PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_Tourette';
PATIENT_ID='TS04 Double DBS Implantation';
RECORDING_DATE='2017_03_01';
RUN_ID='run16';
% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
% PATIENT_ID='ET_CL_001';
% RECORDING_DATE='2017_05_17';
% RUN_ID='run17';
MIDPATH='preproc';
ADDON='_NO_PSA';
NOTES='Same HPF as previously used, but no Right PSA to cloud the results';

cues_only=true;
extrap_method='';
alpha=0.05;

%% Load Data

FILE=fullfile(PREPATH,PATIENT_ID,RECORDING_DATE,MIDPATH,RUN_ID);

% config_load=struct('preset',1);
% [channels,labels,conditions,cond_labels,visit_type]=load_variables(PATIENT_ID,RECORDING_DATE,RUN_ID,config_load);

[channels,labels,conditions,cond_labels,visit_type]=load_variables(PATIENT_ID,RECORDING_DATE,RUN_ID);

if isempty(channels) || isempty(labels) || isempty(conditions) || isempty(cond_labels)
    return
end

datastorage=load(FILE,'datastorage');
fs_init=extract_sampling_frequency(datastorage);

filtering=struct;

order_hp=3;
cutoff_hp=6;
[filtering.hpf.num,filtering.hpf.den]=CreateHPF_butter(fs_init,order_hp,cutoff_hp);

filtering.downsample=200;
filtering.normalize='z-score';
realizationLengthInSeconds=1;
filtering.realizations.length=floor(realizationLengthInSeconds*fs_init);

order_notch=3;
cutoff_notch=[58,62];
[filtering.notch.num,filtering.notch.den]=CreateBSF_butter(fs_init,order_notch,cutoff_notch);
[filtering.lpf.num,filtering.lpf.den]=CreateLPF_butter(fs_init,8,round(filtering.downsample/2));

[x_all,fs]=load_data(datastorage,channels,[],filtering,visit_type,[],extrap_method);

numChannels=size(x_all,2);
numConditions=length(conditions);
numRealizations=nan(numConditions,1);

x=struct;
ar=struct;
res=struct; 
crit=struct;
h=struct;
pVal=struct;
gamma=struct;
avg_psd=struct;
avg_gamma=struct;
pass=struct;

surrogate=nan;
distribution=nan;
pxx=nan;
pxx_filt=nan;

ar_filt=struct;
res_filt=struct;
crit_filt=struct;
h_filt=struct;
pVal_filt=struct;
gamma_filt=struct;

freqForAnalysis=4:0.5:100;

config_mvar=struct(...
    'orderSelection','diff1',...
    'crit','bic',...
    'method','arfit',...
    'orderRange',1:20,...   % Used to find the optimal model order
    'modelOrder',[],...     % This is the optimal model order
    'fs',fs,...
    'epsilon',0.001,...
    'freqRange',freqForAnalysis,...
    'output',0,...          % Since I'm using the optimal model order, don't output the search in mvar
    'logLikelihoodMethod',2);   % 1 == Matlab, 2 == Ding, 3 == Awareness Paper (see calculate_bic.m)
config_plot=struct(...
    'hFig',[],...
    'seriesType',1,...
    'plotType','avgerr',...
    'freqLims',freqForAnalysis,...
    'fs',fs);
config_surr=struct(...
    'method','sample',...
    'numSamples',1,...
    'iterations',1000);

%% Main Loop for initial processing

fprintf('First pass beginning...\n');

for j=1:numConditions
    currCond=cond_labels{j};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    [x.(currCond),~]=load_data(datastorage,channels,conditions(j),filtering,visit_type,cues_only);
    numRealizations(j)=size(x.(currCond),3);
    
    numTrials=size(x.(currCond),3);
    numChannels=length(channels);
    numSamples=length(x.(currCond)(:,1));

    h.(currCond)=nan(numTrials,numChannels);
    pVal.(currCond)=nan(numTrials,numChannels);
    pass.(currCond)=nan(numTrials,1);

    gamma.(currCond)=nan;
    
    fprintf('Calculating optimal order...');

    [config_mvar.modelOrder,crit.(currCond)]=calculate_optimal_model_order(x.(currCond), config_mvar);
    
    fprintf('Optimal Model Order = %d\n',config_mvar.modelOrder);

    for i=1:numTrials
        fprintf('%d/%d\n',i,numTrials);
        
        %% Calculate MVAR model
        [ar.(currCond)(i).mdl, res.(currCond)(i).E] = mvar(squeeze(x.(currCond)(:,:,i)), config_mvar);
        
        %% Test Whiteness
        [pass.(currCond)(i), h.(currCond)(i,:), pVal.(currCond)(i,:)]=...
            test_model(res.(currCond)(i).E,length(x.(currCond)(:,:,i)));
        
        %% Calculate DTF Connectivity (No DTF for correlated data)
    end
end

fprintf('First pass completed.\n');

%% Surrogate Analysis (No surrogate for correlated data)

%% Decorrelation

fprintf('Decorrelating time series channels...\n');
[x_filt,filt_values]=filter_serial_correlation(x,res,h,config_mvar);
fprintf('Decorrelation completed.\n');

%% Main Loop for decorrelated data

fprintf('Decorrelated data connectivity calculations beginning...\n');

for i=1:length(cond_labels)
    currCond=cond_labels{i};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    numTrials=size(x_filt.(currCond),3);
    
    h_filt.(currCond)=nan(numTrials,numChannels);
    pVal_filt.(currCond)=nan(numTrials,numChannels);
    pass_filt.(currCond)=ones(numTrials,1);

    gamma_filt.(currCond)=zeros(numChannels,numChannels,length(freqForAnalysis),numTrials);

    [config_mvar.modelOrder,crit_filt.(currCond)]=calculate_optimal_model_order(x_filt.(currCond), config_mvar);

    fprintf('Optimal Model Order = %d\n',config_mvar.modelOrder);

    for j=1:numTrials
        if filt_values.(currCond)(j).decorrelated
            fprintf('%d/%d\n',j,numTrials);

            %% Calculate MVAR model
            [ar_filt.(currCond)(j).mdl,res_filt.(currCond)(j).E]=mvar(squeeze(x_filt.(currCond)(:,:,j)),config_mvar);

            %% Test whiteness
            [pass_filt.(currCond)(j),h_filt.(currCond)(j,:),pVal_filt.(currCond)(j,:)]=...
                test_model(res_filt.(currCond)(j).E,length(x_filt.(currCond)(:,:,j)));

            %% Calculate DTF Connectivity
            gamma_filt.(currCond)(:,:,:,j)=dtf(ar_filt.(currCond)(j).mdl,freqForAnalysis,fs);
        else
            warning('Trial %d of %d was not properly decorrelated',j,numTrials)
        end
    end
    
    print_whiteness(h_filt.(currCond),pVal_filt.(currCond),labels);
end

fprintf('Decorrelated data connectivity calculations completed.\n');

%% Surrogate Analysis on Decorrelated Data

fprintf('Surrogate analysis of decorrelated data beginning...\n');
[surrogate_filt,distribution_filt]=surrogate_analysis(x_filt,fs,freqForAnalysis,config_mvar,config_surr);
fprintf('Surrogate analyis of decorrelated data completed.\n');

%% Save all relevant variables

[newFile,subjID,dateID,runID]=simplify_filename(PATIENT_ID,RECORDING_DATE,RUN_ID,ADDON);

contactNames=get_structure_names(subjID);

save(newFile,'ADDON','ar','channels','conditions','cond_labels','crit',...
    'FILE','freqForAnalysis','filtering','fs','fs_init','gamma','h','labels',...
    'newFile','pass','PATIENT_ID','pVal','RECORDING_DATE','res','RUN_ID','x','x_all',...
    'config_mvar','config_plot','config_surr','NOTES','surrogate','distribution','pxx',...
    'x_filt','filt_values','ar_filt','res_filt','crit_filt','h_filt','pVal_filt',...
    'gamma_filt','surrogate_filt','distribution_filt','pxx_filt','visit_type','cues_only',...
    'extrap_method','subjID','dateID','runID','contactNames','alpha',...
    'numChannels','numRealizations');









