%% single_patient_all_conditions
%
%  Collect all connectivity values for all conditions in one run of one patient
%
%   See also: mvar, plot_connectivity, dtf, load_data

CCC;

%% Load all data

PREPATH='\\gunduz-lab.bme.ufl.edu\\Study_ET_Closed_Loop';
PATIENT_ID='ET_CL_004';
RECORDING_DATE='2018_06_20';
MIDPATH='preproc';
RUN_ID='run5';
ADDON='__ALLCOND_ALLCOMB';
FILE=fullfile(PREPATH,PATIENT_ID,RECORDING_DATE,MIDPATH,RUN_ID);
[channels,labels,conditions,cond_labels]=load_channels_labels_conditions(PATIENT_ID,RECORDING_DATE,RUN_ID);
% channels=[8,7;6,5;4,3;2,1];
% labels={'Vim (3-2)','Vim (1-0)','Cort (3-2)','Cort (1-0)'};
% conditions=[1,2,3,4,5]; % 1 == Rest, 2 == cue right, 3 == cue left, 4 == move right, 5 == move left
% cond_labels={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};
fs_filtering=extract_sampling_frequency(FILE);

order_notch=4;
cutoff_notch=[54,66;114,126;176,184;236,244];

filtering=struct('notch',[]);

for i=1:length(cutoff_notch)
    [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(fs_filtering,order_notch,cutoff_notch(i,:));
end

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
    
    [x.(currCond),fs]=load_data(FILE,channels,conditions(j),filtering);
    [x_all.(currCond),~]=load_data(FILE,channels,[],filtering);

    numTrials=size(x.(currCond),3);
    numChannels=length(channels);

    %% Calculate all MVAR models

    for i=1:numTrials
        fprintf('%d - ',i);
        [ar(i).(currCond).estMdl,res(i).(currCond).E,crit(i).(currCond).BIC]=mvar(squeeze(x.(currCond)(:,:,i)));
    end

    %% Test whiteness

    h.(currCond)=ones(numTrials,numChannels);
    pVal.(currCond)=zeros(numTrials,numChannels);
    pass.(currCond)=0;

    for i=1:numTrials
        [pass.(currCond),h.(currCond)(i,:),pVal.(currCond)(i,:)]=test_model(res(i).(currCond).E,length(x.(currCond)(:,:,i)));

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
        gamma.(currCond)(:,:,:,i)=dtf(ar(i).(currCond).estMdl,freqRange,fs);
    end

    %% Plot all connectivities and PSDs stacked

    config=struct('hFig',[],'seriesType',1);
    config.hFig=figure;
    config.figTitle=sprintf('%s, %s, %s - %s: Individual',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);

    for i=1:numTrials
        plot_connectivity(gamma.(currCond)(:,:,:,i),squeeze(x.(currCond)(:,:,i)),freqRange,labels,config);
    end

    %% Calculate and plot averages of connectivities and PSDs

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

    config.seriesType=3;
    config.hFig=figure;
    config.figTitle=sprintf('%s, %s, %s - %s: Average',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);

    plot_connectivity(avg_gamma.(currCond),avg_psd.(currCond),freqRange,labels,config);
end

save(fullfile(pwd,sprintf('%s__%s__%s%s.mat',PATIENT_ID,RECORDING_DATE,RUN_ID,ADDON)),'ar',...
    'avg_gamma','avg_psd','condition','crit','FILE','freqRange','fs','gamma','h','labels','res','x');


















