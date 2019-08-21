CCC;

%% Load Rat Hippocampal Data

data=load('\\gunduz-lab.bme.ufl.edu\Lab\Sarah\Classes\Multivariate\Project\Ding_Data\laminar_rat_3.7mm_0.7mA_bipolar_CA1_DG.mat');
X=data.XX;
numChannels=size(X,1);
numRealizations=size(X,2);
fs=200;

[num,den]=CreateBSF_butter(fs,3,[58 62]);

for i=1:numRealizations
    X(:,i,:)=filtfilt(num,den,squeeze(X(:,i,:))')';
    
    X(1,i,:)=X(1,i,:)./std(X(1,i,:));
    X(2,i,:)=X(2,i,:)./std(X(2,i,:));
end

freqRange=0:0.1:15;
spectral_range=0:0.05:100;
orderRange=1:50;
labels={'CA1','DG'};

config=struct(...
    'orderRange',orderRange,...
    'crit','aic',...
    'method','arfit',...
    'fs',fs,...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

mdl=struct;
E=struct;
crit=struct;

%% Calculate criterion values for each type

fprintf('Starting AIC...\n');

for i=1:numRealizations
    [~,~,crit(i).aic]=mvar(squeeze(X(:,i,:))',config);
end

config.crit='bic';

fprintf('Starting BIC...\n');

for i=1:numRealizations
    [~,~,crit(i).bic]=mvar(squeeze(X(:,i,:))',config);
end

config.crit='spectra';
config.spectral_range=spectral_range;
config.normalizeSpectra=true;

fprintf('Starting MSE...\n');

for i=1:numRealizations
    [~,~,crit(i).mse]=mvar(squeeze(X(:,i,:))',config);
end

%% Calculate average criterion for each type 

avg=struct;
tmp_sum1=zeros(length(orderRange),1);
tmp_sum2=zeros(length(orderRange),1);
tmp_sum3=zeros(length(orderRange),1);

for i=1:numRealizations
    tmp_sum1=tmp_sum1+crit(i).aic;
    tmp_sum2=tmp_sum2+crit(i).bic;
    tmp_sum3=tmp_sum3+crit(i).mse;
end

avg.aic=tmp_sum1./numRealizations;
avg.bic=tmp_sum2./numRealizations;
avg.mse=tmp_sum3./numRealizations;

%% Plot criterion values

figure; 
subplot(311); plot(orderRange,avg.aic); ylabel('AIC');
subplot(312); plot(orderRange,avg.bic); ylabel('BIC');
subplot(313); plot(orderRange,avg.mse); ylabel('MSE'); xlabel('Frequency [Hz]');

%% Calculate AR models using model order with the minimum average 

[~,min_aic]=min(avg.aic); % AIC sucks, don't use this one at all
[~,order_bic]=min(avg.bic);
[~,order_mse]=min(avg.mse);

config=struct(...
    'orderRange',order_bic,...
    'crit','bic',...
    'method','arfit',...
    'fs',fs,...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

fprintf('Starting BIC...\n');

for i=1:numRealizations
    [mdl(i).bic,E(i).bic,~]=mvar(squeeze(X(:,i,:))',config);
end

config.crit='spectra';
config.orderRange=order_mse;
% config.spectral_range=spectral_range;
config.normalizeSpectra=true;

fprintf('Starting MSE...\n');

for i=1:numRealizations
    [mdl(i).mse,E(i).mse,~]=mvar(squeeze(X(:,i,:))',config);
end

%% Calculate connectivity values for BIC and MSE

gamma=struct('bic',nan(numChannels,numChannels,length(freqRange),numRealizations),...
    'mse',nan(numChannels,numChannels,length(freqRange),numRealizations));
H=struct('bic',nan(numChannels,numChannels,length(freqRange),numRealizations),...
    'mse',nan(numChannels,numChannels,length(freqRange),numRealizations));
S=struct('bic',nan(numChannels,numChannels,length(freqRange),numRealizations),...
    'mse',nan(numChannels,numChannels,length(freqRange),numRealizations));

for i=1:numRealizations
    [gamma.bic(:,:,:,i),H.bic(:,:,:,i)]=dtf(mdl(i).bic,freqRange,fs);
    [gamma.mse(:,:,:,i),H.mse(:,:,:,i)]=dtf(mdl(i).mse,freqRange,fs);
    
    S.bic(:,:,:,i)=calculate_spectra(H.bic(:,:,:,i),mdl(i).bic.C);
    S.mse(:,:,:,i)=calculate_spectra(H.mse(:,:,:,i),mdl(i).mse.C);
end

%% Calculate average FFT values, and average Spectral values

tmp_fft_sum=zeros(size(X,3)/2+1,numChannels);
tmp_mse_sum=zeros(size(mdl(1).mse.pxx));
tmp_S_bic_sum=zeros(numChannels,numChannels,length(freqRange));
tmp_S_mse_sum=zeros(numChannels,numChannels,length(freqRange));

for i=1:numRealizations
    [tmp_fft,~]=calculate_fft(squeeze(X(:,i,:))',fs);
    tmp_fft_sum=tmp_fft_sum+tmp_fft;
    tmp_mse_sum=tmp_mse_sum+mdl(i).mse.pxx;
    tmp_S_bic_sum=tmp_S_bic_sum+S.bic(:,:,:,i);
    tmp_S_mse_sum=tmp_S_mse_sum+S.mse(:,:,:,i);
end

avg_fft=tmp_fft_sum./numRealizations;
avg_mse=tmp_mse_sum./numRealizations;

avg_S_bic=tmp_S_bic_sum./numRealizations;
avg_S_mse=tmp_S_mse_sum./numRealizations;

%% Plot the connectivity and spectral values for BIC

config_plot=struct('seriesType',7,'fs',fs,'plotType','avgerr');
plot_connectivity(gamma.bic,S.bic,freqRange,labels,config_plot);

%% Perform surrogate analysis on the realizations to calculate empirical distribution

file=struct;
file.fs=fs;
file.config_plot=config_plot;
file.config_plot.freqLims=freqRange;
file.freqRange=freqRange;
file.x.rat=permute(X,[3,1,2]);
file.config_crit=config;
file.ar.rat.mdl=mdl.bic;
file.newFile=fullfile(get_root_path,'Files','rat_hippocampus.mat');

surrogate=struct;
distribution=struct;

config_surr=struct('iterations',nchoosek(numRealizations,numChannels) * 2);

[surrogate.bic,distribution.bic]=surrogate_analysis(file,config_surr);
surrogate.bic=surrogate.bic.rat;
distribution.bic=distribution.bic.rat;

file.ar.rat.mdl=mdl.mse;

[surrogate.mse,distribution.mse]=surrogate_analysis(file,config_surr);
surrogate.mse=surrogate.mse.rat;
distribution.mse=distribution.mse.rat;

%% Plot the surrogate values for BIC

series=struct('spectral',S.bic,'surrogate',surrogate.bic);
config_plot.surr_params.highlightSignificance=true;
config_plot.surr_params.threshold=0.01;
plot_connectivity(gamma.bic,series,freqRange,labels,config_plot);

%% Plot a histogram of surrogate values obtained for BIC

% Taken from Evaluating causal relations in neural systems: Granger causality, directed transfer function and statistical assessment of significance

allGammaValues=surrogate.bic(2,1,:,:);
allGammaValues=allGammaValues(:);


