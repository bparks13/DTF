%% single_patient_one_condition
%
%  Runs DTF on all trials from one condition in a single patient
%

CCC;

%% Load all trials

FILE='\\gunduz-lab.bme.ufl.edu\Study_ET_Closed_Loop\ET_CL_004\2018_06_20\preproc\run5.mat';
channels=[8,7;4,3;2,1];
labels={'Vim (3-2)','Cortex (3-2)','Cortex (1-0)'};
condition=1; % 1 == Rest

order_notch=4;
cutoff_notch=[54,66;114,126;176,184;236,244];

filtering=struct('notch',[]);

for i=1:length(cutoff_notch)
    [filtering.notch(i).num,filtering.notch(i).den]=CreateBSF_butter(2400,order_notch,cutoff_notch(i,:));
end

order_lp=4;
cutoff_lp=800;
[num_lp,den_lp]=CreateLPF_butter(2400,order_lp,cutoff_lp);

filtering.lpf.num=num_lp;
filtering.lpf.den=den_lp;

[x,fs]=load_data(FILE,channels,condition,filtering);
[x_all,~]=load_data(FILE,channels,[],filtering);

numTrials=size(x,3);
numChannels=length(channels);

%% Calculate all MVAR models

ar=struct('estMdl',[]);
res=struct('E',[]); % residuals (E from estimate)
crit=struct('BIC',[]);

for i=1:numTrials
    fprintf('%d - ',i);
    [ar(i).estMdl,res(i).E,crit(i).BIC]=mvar(squeeze(x(:,:,i)));
end

%% Test whiteness

h=ones(numTrials,numChannels);
pVal=zeros(numTrials,numChannels);

for i=1:numTrials
    [pass,h(i,:),pVal(i,:)]=test_model(res(i).E,length(x(:,:,i)));
    
    if ~pass
        fprintf('WARNING: Null hypothesis of uncorrelated errors rejected for trial %d\n',i);
        for j=1:numChannels
            fprintf('\t%s: h = %d with p = %.4f\n',labels{j},h(i,j),pVal(i,j));
        end
    end
end

%% Calculate DTF

freqRange=2:50;

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

%% Save a file containing relevant information

save('ET_CL_004__2018_06_20_5','ar','avg_gamma','avg_psd','condition','crit','FILE','freqRange',...
    'fs','gamma','h','labels','res','x');