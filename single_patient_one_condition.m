%% single_patient_one_condition
%
%  Runs DTF on all trials from one condition in a single patient
%

CCC;

%% Load all trials

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
channels=[7,6;4,3];
labels={'Vim (2-1)','Cortex (3-2)'};
condition=1; % 1 == Rest

order_lp=8;
cutoff_lp=200;
[num_lp,den_lp]=CreateLPF_butter(fs,order_lp,cutoff_lp);

filtering=struct('lpf',[]);
filtering.lpf.num=num_lp;
filtering.lpf.den=den_lp;

[x,fs]=load_data(FILE,channels,condition,filtering);

numTrials=size(x,3);
numChannels=length(channels);

%% Calculate all MVAR models

ar=struct('estMdl',[]);
res=struct('E',[]); % residuals (E from estimate)
crit=struct('BIC',[]);

for i=1:numTrials
    [ar(i).estMdl,res(i).E,crit(i).BIC]=mvar(squeeze(x(:,:,i)));
end

%% Test whiteness

h=ones(numTrials,1);

for i=1:numTrials
    h(i)=test_model(res(i).E,length(x(:,:,i)));
    
    if h(i)
        fprintf('WARNING: Null hypothesis of uncorrelated errors rejected for trial %d\n',i);
    end
end

%% Calculate DTF

freqRange=2:100;

gamma=zeros(numChannels,numChannels,length(freqRange),numTrials);

for i=1:numTrials
    gamma(:,:,:,i)=dtf(ar(i).estMdl,freqRange,fs);
end

%% Plot all connectivities and PSDs stacked

config=struct('hFig',[],'seriesType',1);
config.hFig=figure;

for i=1:numTrials
    plot_connectivity(gamma(:,:,:,i),squeeze(x(:,:,i)),freqRange,labels,config);
end

%% Calculate and plot averages of connectivities and PSDs

avg_psd=zeros(length(freqRange),numChannels);
avg_gamma=zeros(numChannels,numChannels,length(freqRange));

window=round(fs);
overlap=round(window/2);

for i=1:numTrials
    avg_psd=avg_psd+pwelch(x(:,:,i),window,overlap,freqRange,fs);
    avg_gamma=avg_gamma+gamma(:,:,:,i);
end

avg_psd=avg_psd/numTrials;
avg_gamma=avg_gamma/numTrials;

config.seriesType=3;
config.hFig=figure;

plot_connectivity(avg_gamma,avg_psd,freqRange,labels,config);