%% single_patient_all_conditions
%
%  Collect all connectivity values for all conditions in one run of one patient
%
%  See also: mvar, plot_connectivity, dtf, load_data
%

CCC;

%% Load all data

PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
PATIENT_ID='ET_CL_004';
% PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_Tourette';
% PATIENT_ID='TS04 Double DBS Implantation';
RECORDING_DATE='2018_06_20';
MIDPATH='preproc';
RUN_ID='run5';
ADDON='__PSD__Z_SCORE';
FILE=fullfile(PREPATH,PATIENT_ID,RECORDING_DATE,MIDPATH,RUN_ID);
% config=struct('default',false,'preset',2);
% [channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID,config);
[channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID);

if isempty(channels) || isempty(labels) || isempty(conditions) || isempty(cond_labels)
    return
end

fs_init=extract_sampling_frequency(FILE);

filtering=struct;
filtering.NO_FILTERING=true;
[filtering.hpf.num,filtering.hpf.den]=CreateHPF_butter(fs_init,3,2);
% filtering.downsample=400;
filtering.normalize='z-score';

% order_notch=4;
% cutoff_notch=[58,62;118,122];
% 
% for i=1:length(cutoff_notch)
%     [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(fs,order_notch,cutoff_notch(i,:));
% end

% [filtering.lpf.num,filtering.lpf.den]=CreateLPF_butter(fs,8,600);

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
avg_psd=struct;
avg_gamma=struct;
pass=struct;

freqForAnalysis=4:200;
config_crit=struct('orderSelection','min','crit','psd','orderRange',1:100,'fs',fs,...
    'freqRange',freqForAnalysis);
config_plot=struct('hFig',[],'seriesType',1,'plotType','avgerr','freqLims',freqForAnalysis,...
    'fs',fs);

%% Main Loop

for j=1:numConditions
    currCond=cond_labels{j};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    [x.(currCond),fs]=load_data(FILE,channels,conditions(j),filtering);

    numTrials=size(x.(currCond),3);
    numChannels=length(channels);
    numSamples=length(x.Rest(:,1));

    %% Calculate all MVAR models

    for i=1:numTrials
        fprintf('%d/%d - ',i,numTrials);
        [ar.(currCond)(i).mdl,res.(currCond)(i).E,crit.(currCond)(i).(config_crit.crit)]=mvar(squeeze(x.(currCond)(:,:,i)),config_crit);
    end

    %% Test whiteness

    h.(currCond)=nan(numTrials,numChannels);
    pVal.(currCond)=nan(numTrials,numChannels);
    pass.(currCond)=nan(numTrials,1);

    for i=1:numTrials
        [pass.(currCond)(i),h.(currCond)(i,:),pVal.(currCond)(i,:)]=test_model(res.(currCond)(i).E,length(x.(currCond)(:,:,i)));
    end
    
    print_whiteness(h.(currCond),pVal.(currCond),labels);

    %% Calculate DTF

    gamma.(currCond)=zeros(numChannels,numChannels,length(freqRange),numTrials);

    for i=1:numTrials
        gamma.(currCond)(:,:,:,i)=dtf(ar.(currCond)(i).mdl,freqRange,fs);
    end

    %% Plot all connectivities and PSDs

    config_plot.hFig=figure;
    config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
    
    plot_connectivity(gamma.(currCond),x.(currCond),freqRange,labels,config_plot);
end

%% Save relevant variables

count=1;

newFile=fullfile(pwd,sprintf('%s__%s__%s%s.mat',PATIENT_ID,RECORDING_DATE,RUN_ID,ADDON));

while exist(newFile,'file')
    newFile=fullfile(pwd,sprintf('%s__%s__%s%s_(%d).mat',PATIENT_ID,RECORDING_DATE,RUN_ID,ADDON,count));
    count=count+1;
end

config_crit.hFig=[];
config_plot.hFig=[];
save(newFile,'ar','avg_gamma','avg_psd','channels','conditions','cond_labels','crit',...
    'FILE','freqRange','filtering','fs','fs_init','gamma','h','labels','PATIENT_ID','pVal',...
    'RECORDING_DATE','res','RUN_ID','x','x_all','config_crit','config_plot');


















