%% pipeline.m
%
%  Combines both the processing_pipeline.m and postprocessing_pipeline.m files, so that
%  all processing is done at the same time.
%
%  See also: dtf, mvar, load_data, load_channels_labels_conditions, test_model,
%   simplify_filename, surrogate_analysis, filter_serial_correlation, test_model
%

CCC;
dtf_startup;

%% Load Data

% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET\\OR\\with_DBS';
% PATIENT_ID='ET_OR_STIM_018';
% RECORDING_DATE='2018_11_28';
% RUN_ID='run12';
PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
PATIENT_ID='ET_CL_004';
RECORDING_DATE='2018_06_20';
RUN_ID='run5';
% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_Tourette';
% PATIENT_ID='TS04 Double DBS Implantation';
MIDPATH='preproc';
ADDON='_BIC_SAMPLE_SURROGATE';
NOTES='Trying a longer realization length';
FILE=fullfile(PREPATH,PATIENT_ID,RECORDING_DATE,MIDPATH,RUN_ID);
% config=struct('default',false,'preset',2);
% [channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID,config);
[channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID);

if isempty(channels) || isempty(labels) || isempty(conditions) || isempty(cond_labels)
    return
end

fs_init=extract_sampling_frequency(FILE);

filtering=struct;
[filtering.hpf.num,filtering.hpf.den]=CreateHPF_butter(fs_init,3,2);
filtering.downsample=200;
filtering.normalize='z-score';
realizationLengthInSeconds=2.5;
filtering.realizations.length=realizationLengthInSeconds*fs_init;

order_notch=4;
cutoff_notch=[58,62];
[filtering.notch.num,filtering.notch.den]=CreateBSF_butter(fs_init,order_notch,cutoff_notch);
[filtering.lpf.num,filtering.lpf.den]=CreateLPF_butter(fs_init,8,round(filtering.downsample/2));

[x_all,fs]=load_data(FILE,channels,[],filtering);

numConditions=length(conditions);
freqRange=1:(fs/2);

x=struct;
ar=struct;
res=struct; 
crit=struct;
h=struct;
pVal=struct;
gamma=struct;
H=struct;
avg_psd=struct;
avg_gamma=struct;
pass=struct;

freqForAnalysis=4:0.5:100;

config_crit=struct(...
    'orderSelection','diff1',...
    'crit','bic',...
    'method','arfit',...
    'orderRange',1:35,...
    'fs',fs,...
    'freqRange',freqForAnalysis,...
    'logLikelihoodMethod',2);   % 1 == Matlab, 2 == Ding, 3 == Awareness Paper (see calculate_bic.m)
config_plot=struct(...
    'hFig',[],...
    'seriesType',1,...
    'plotType','avgerr',...
    'freqLims',freqForAnalysis,...
    'fs',fs);
config_surr=struct(...
    'method','sample',...
    'numSamples',1);

%% Main Loop for initial processing

for j=1:numConditions
    currCond=cond_labels{j};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    [x.(currCond),~]=load_data(FILE,channels,conditions(j),filtering);

    numTrials=size(x.(currCond),3);
    numChannels=length(channels);
    numSamples=length(x.Rest(:,1));

    h.(currCond)=nan(numTrials,numChannels);
    pVal.(currCond)=nan(numTrials,numChannels);
    pass.(currCond)=nan(numTrials,1);

    gamma.(currCond)=nan(numChannels,numChannels,length(freqRange),numTrials);
    H.(currCond)=nan(numChannels,numChannels,length(freqRange),numTrials);

    for i=1:numTrials
        fprintf('%d/%d - ',i,numTrials);
        
        %% Calculate MVAR model
        [ar.(currCond)(i).mdl, res.(currCond)(i).E, crit.(currCond)(i).(config_crit.crit)]=...
            mvar(squeeze(x.(currCond)(:,:,i)), config_crit);
        
        %% Test Whiteness
        [pass.(currCond)(i), h.(currCond)(i,:), pVal.(currCond)(i,:)]=...
            test_model(res.(currCond)(i).E,length(x.(currCond)(:,:,i)));
        
        %% Calculate DTF Connectivity
        [gamma.(currCond)(:,:,:,i), H.(currCond)(:,:,:,i)]=dtf(ar.(currCond)(i).mdl,freqRange,fs);
    end
    
    print_whiteness(h.(currCond),pVal.(currCond),labels);
end

%% Surrogate Analysis

surrogate_analysis(data,config_surr);

%% Decorrelation

[x_filt,filt_values]=filter_serial_correlation(x,res,h,config_crit);

%% Main Loop for decorrelated data

ar_filt=struct;
res_filt=struct;
crit_filt=struct;
h_filt=struct;
pVal_filt=struct;
gamma_filt=struct;

for i=1:length(cond_labels)
    currCond=cond_labels{i};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    numTrials=size(x_filt.(currCond),3);
    
    h_filt.(currCond)=nan(numTrials,numChannels);
    pVal_filt.(currCond)=nan(numTrials,numChannels);
    pass_filt.(currCond)=zeros(numTrials,1);

    gamma_filt.(currCond)=zeros(numChannels,numChannels,length(freqRange),numTrials);

    for j=1:numTrials
        if data.filt_values.(currCond)(j).decorrelated
            fprintf('%d/%d - ',j,numTrials);

            % Calculate all MVAR models
            [ar_filt.(currCond)(j).mdl,res_filt.(currCond)(j).E,crit_filt.(currCond)(j).(config_crit.crit)]=...
                mvar(squeeze(data.x_filt.(currCond)(:,:,j)),config_crit);

            % Test whiteness
            [pass_filt.(currCond),h_filt.(currCond)(j,:),pVal_filt.(currCond)(j,:)]=...
                test_model(res_filt.(currCond)(j).E,length(data.x_filt.(currCond)(:,:,j)));

            % Calculate DTF
            gamma_filt.(currCond)(:,:,:,j)=dtf(ar_filt.(currCond)(j).mdl,freqRange,fs);
        end
    end
    
    print_whiteness(h_filt.(currCond),pVal_filt.(currCond),labels);
end











