%% single_patient_all_conditions
%
%  Collect all connectivity values for all conditions in one run of one patient
%
%   See also: mvar, plot_connectivity, dtf, load_data

CCC;

%% Load all data

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
channels=[8,7;4,3];
labels={'Vim (3-2)','Cortex (3-2)'};
condition=[1,2,3,4,5]; % 1 == Rest, 2 == cue right, 3 == cue left, 4 == move right, 5 == move left
cond_label={'Rest','CueRight','CueLeft','MoveRight','MoveLeft'};

order_notch=4;
cutoff_notch=[54,66;114,126;176,184;236,244];

filtering=struct('notch',[]);

for i=1:length(cutoff_notch)
    [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(2400,order_notch,cutoff_notch(i,:));
end

numConditions=length(condition);
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

for j=1:numConditions
    currCond=cond_label{j};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    [x.(currCond),fs]=load_data(FILE,channels,condition(j),filtering);
    [x_all.(currCond),~]=load_data(FILE,channels,[],filtering);

    numTrials=size(x.(currCond),3);
    numChannels=length(channels);

    %% Calculate all MVAR models

    for i=1:numTrials
        fprintf('%d - ',i);
        [ar(i).(currCond).estMdl,res(i).(currCond).E,crit(i).(currCond).BIC]=mvar(squeeze(x.(currCond)(:,:,i)));
    end

    %% Test whiteness

    h.(currCond)=ones(numTrials,1);
    pVal.(currCond)=zeros(numTrials,numChannels);

    for i=1:numTrials
        [h.(currCond)(i),pVal.(currCond)(i,:)]=test_model(res(i).(currCond).E,length(x.(currCond)(:,:,i)));

        if h.(currCond)(i)
            fprintf('WARNING: Null hypothesis of uncorrelated errors rejected for trial %d, with p = %.4f\n',i,pVal(i));
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

    plot_connectivity(avg_gamma.(currCond),avg_psd.(currCond),freqRange,labels,config);
end

save('ET_CL_004__2018_06_20_ALLCOND','ar','avg_gamma','avg_psd','condition','crit','FILE','freqRange',...
    'fs','gamma','h','labels','res','x');


















