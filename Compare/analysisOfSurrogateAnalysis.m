%% Analysis of how surrogate analysis changes depending on the number of samples randomly selected

CCC;

%% Load rat hippocampal data

data=load('\\gunduz-lab.bme.ufl.edu\Lab\Sarah\Classes\Multivariate\Project\Ding_Data\laminar_rat_3.7mm_0.7mA_bipolar_CA1_DG.mat');
X=permute(data.XX,[3 1 2]);
numSamples=size(X,1);
numChannels=size(X,2);
numRealizations=size(X,3);
fs=200;

[num,den]=CreateBSF_butter(fs,3,[58 62]);

for i=1:numRealizations
    X(:,:,i)=filtfilt(num,den,X(:,:,i));
    
    X(:,1,i)=X(:,1,i)./std(X(:,1,i));
    X(:,2,i)=X(:,2,i)./std(X(:,2,i));
end

% Calculate values based on pre-determined model orders

% order=18; % Model order of 18 for BIC
freqRange=0:0.1:15;

config=struct(...
    'orderRange',18,...
    'crit','bic',...
    'method','arfit',...
    'fs',fs,...
    'orderSelection','min',...
    'freqRange',freqRange,...
    'logLikelihoodMethod',2 ...
    );

mdl=struct;
E=struct;
gamma=struct('bic',nan(numChannels,numChannels,length(freqRange),numRealizations));
H=struct('bic',nan(numChannels,numChannels,length(freqRange),numRealizations));
S=struct('bic',nan(numChannels,numChannels,length(freqRange),numRealizations));

for i=1:numRealizations
    [mdl(i).bic,E(i).bic,~]=mvar(X(:,:,i),config);
    [gamma.bic(:,:,:,i),H.bic(:,:,:,i)]=dtf(mdl(i).bic,freqRange,fs);
    S.bic(:,:,:,i)=calculate_spectra(H.bic(:,:,:,i),mdl(i).bic.C);
end

%% Save this file setup

% save(fullfile(get_root_path,'Files','rat_data_for_surrogate_analysis.mat'),'-v7.3');

% load(fullfile(get_root_path,'Files','rat_data_for_surrogate_analysis.mat'));

%% Set values for surrogate analysis with randomized sample replacement

config_plot=struct('seriesType',5,'fs',fs,'plotType','avgerr'); 
config_plot.surr_params.highlightSignificance=true;
labels={'CA1','DG'};
config.orderRange=1:15;

file=struct;
file.fs=fs;
file.config_plot=config_plot;
file.config_plot.freqLims=freqRange;
file.freqRange=freqRange;
file.x.rat=X;
file.config_crit=config;
file.ar.rat.mdl=mdl.bic;
file.newFile=fullfile(get_root_path,'Files','rat_data_for_surrogate_analysis.mat');

surrogate=struct;
distribution=struct;
pxx=struct;

config_surr=struct('iterations',1000,'numSamples',1,'method','sample');

%% Run surrogate analysis with randomized sample replacement

for i=1:numSamples
    config_surr.numSamples=i;
    [tmp_surr,tmp_dist,tmp_pxx]=surrogate_analysis(file,config_surr);
    surrogate(i).bic=tmp_surr.rat;
    distribution(i).bic=tmp_dist.rat;
    pxx(i).bic=tmp_pxx.rat;
end

%% Run surrogate analysis with time shifting

config_surr=struct('method','shift','numIterations',1000);

[shift_surr,shift_dist,shift_pxx]=surrogate_analysis(file,config_surr);

%% Check for uniformity in the distribution

% figure;
% histogram2(tmp_dist.rat(:,1),tmp_dist.rat(:,2),'FaceColor','flat')
% colorbar;

%% Plot the time-shifted connectivity measures

% plot_connectivity(gamma.bic,tmp_surr.rat,freqRange,labels,config_plot);

%%

% config_plot.seriesType=8;

% 
% for i=1:numSamples
%     series=struct('surrogate',surrogate(i).bic,'original_psd',pxx(i).bic);
%     plot_connectivity(gamma.bic,series,freqRange,labels,config_plot);
% %     plot_connectivity(gamma.bic,surrogate(i).bic,freqRange,labels,config_plot);
% end