%% single_patient_all_conditions
%
%  Collect all connectivity values for all conditions in one run of one patient
%
%  See also: mvar, plot_connectivity, dtf, load_data

CCC;

%% Load all data

PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
PATIENT_ID='ET_CL_004';
RECORDING_DATE='2018_06_20';
MIDPATH='preproc';
RUN_ID='run5';
ADDON='__ALLCOND_ALLCOMB_TESTINGYW';
FILE=fullfile(PREPATH,PATIENT_ID,RECORDING_DATE,MIDPATH,RUN_ID);
config=struct('default',false,'preset',1);
[channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID,config);
% [channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID);

if isempty(channels) || isempty(labels) || isempty(conditions) || isempty(cond_labels)
    return
end

fs=extract_sampling_frequency(FILE);

% order_notch=4;
% cutoff_notch=[54,66;114,126;176,184;236,244];

filtering=struct('NO_FILTERING',true);
% [filtering.hpf.num,filtering.hpf.den]=CreateHPF_butter(fs,4,4);

% for i=1:length(cutoff_notch)
%     [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(fs,order_notch,cutoff_notch(i,:));
% end

numConditions=length(conditions);
freqRange=2:50;

x=struct;
x_all=struct;
ar=struct;
res=struct; % residuals (E from estimate)
crit=struct;
h=struct;
pVal=struct;
gamma=struct;
avg_psd=struct;
avg_gamma=struct;
pass=struct;

for j=1:numConditions
    currCond=cond_labels{j};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    [x.(currCond),fs,x_all.(currCond)]=load_data(FILE,channels,conditions(j),filtering);

    numTrials=size(x.(currCond),3);
    numChannels=length(channels);

    %% Calculate all MVAR models

    for i=1:numTrials
        fprintf('%d/%d - ',i,numTrials);
        [ar.(currCond)(i).mdl,res.(currCond)(i).E,crit.(currCond)(i).BIC]=mvar(squeeze(x.(currCond)(:,:,i)));
    end

    %% Test whiteness

    h.(currCond)=ones(numTrials,numChannels);
    pVal.(currCond)=zeros(numTrials,numChannels);
    pass.(currCond)=0;

    for i=1:numTrials
        [pass.(currCond),h.(currCond)(i,:),pVal.(currCond)(i,:)]=test_model(res.(currCond)(i).E,length(x.(currCond)(:,:,i)));

        if ~pass.(currCond)
            fprintf('WARNING: Null hypothesis of uncorrelated errors rejected for trial %d\n',i);
            for k=1:numChannels
                fprintf('\t%s: h = %d with p = %.4f\n',labels{k},h.(currCond)(i,k),pVal.(currCond)(i,k));
            end
        end
    end

    %% Calculate DTF

    gamma.(currCond)=zeros(numChannels,numChannels,length(freqRange),numTrials);

    for i=1:numTrials
        gamma.(currCond)(:,:,:,i)=dtf(ar.(currCond)(i).mdl,freqRange,fs);
    end

    %% Plot all connectivities and PSDs stacked

    config=struct('hFig',[],'seriesType',1);
    config.hFig=figure;
    config.figTitle=sprintf('%s, %s, %s - %s: Individual',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);

    for i=1:numTrials
        plot_connectivity(gamma.(currCond)(:,:,:,i),squeeze(x.(currCond)(:,:,i)),freqRange,labels,config);
    end
    
    % FOR plot_connectivity, ADD IN A CONFIG OPTION TO PLOT ALL, PLOT AVERAGES, AND PLOT
    % SHADED ERROR BARS

    %% Calculate averages of connectivities and PSDs

    avg_psd.(currCond)=zeros(length(freqRange),numChannels);
    avg_gamma.(currCond)=zeros(numChannels,numChannels,length(freqRange));

    window=round(fs);
    overlap=round(window/2);

    for i=1:numTrials
        avg_psd.(currCond)=avg_psd.(currCond)+pwelch(x.(currCond)(:,:,i),window,overlap,freqRange,fs);
        avg_gamma.(currCond)=avg_gamma.(currCond)+gamma.(currCond)(:,:,:,i);
    end

    avg_psd.(currCond)=avg_psd.(currCond)/numTrials;
    avg_gamma.(currCond)=avg_gamma.(currCond)/numTrials;

    %% Plot averages of connectivities and PSDs
    
    config.seriesType=3;
    config.hFig=figure;
    config.figTitle=sprintf('%s, %s, %s - %s: Average',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);

    plot_connectivity(avg_gamma.(currCond),avg_psd.(currCond),freqRange,labels,config);
    
    %% Plot the criterion of all the trials
    
    config=struct;
    
    config.hFig=figure;
    
    for i=1:numTrials
        plot_criterion(crit.(currCond)(i).BIC,config);
    end
end

%% Save relevant variables

count=1;

newFile=fullfile(pwd,sprintf('%s__%s__%s%s.mat',PATIENT_ID,RECORDING_DATE,RUN_ID,ADDON));

while exist(newFile,'file')
    newFile=fullfile(pwd,sprintf('%s__%s__%s%s_(%d).mat',PATIENT_ID,RECORDING_DATE,RUN_ID,ADDON,count));
    count=count+1;
end

save(newFile,'ar','avg_gamma','avg_psd','channels','conditions','cond_labels','crit',...
    'FILE','freqRange','filtering','fs','gamma','h','labels','res','x','x_all');


















