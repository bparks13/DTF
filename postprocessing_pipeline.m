%% postprocessing_pipeline
%
%  Using the file already created, runs additional processing steps on the data, such as
%  surrogate analysis and decorrelating residuals. Then, runs the connectivity analysis on
%  the filtered data, and runs surrogate on that filtered data 
%
%  See also: processing_pipeline, surrogate_analysis, filter_serial_correlation, mvar,
%   test_model, dtf
%

CCC;

%% Set file

file=uigetfile(fullfile(get_root_path,'Files','*.mat'));

if file~=0
    data=load(fullfile(get_root_path,'Files',file));
end

%% Surrogate analysis

% surrogate_analysis(fullfile(get_root_path,'Files',file));   % run and save the surrogate analysis for all conditions
if isfield(data,'config_surr')
    surrogate_analysis(data,data.config_surr); 
elseif exist('config_surr','var') == 1
    surrogate_analysis(data,config_surr);   % run and save the surrogate analysis for all conditions
else
    surrogate_analysis(data);
end

%% Decorrelation

% [x_filt,filt_values]=filter_serial_correlation(file); % filter all conditions
[x_filt,filt_values]=filter_serial_correlation(data); % filter all conditions
save(fullfile(get_root_path,'Files',file),'-append','x_filt','filt_values');
data.x_filt=x_filt;
data.filt_values=filt_values;

%% Run analyses on filtered data

ar_filt=struct;
res_filt=struct;
crit_filt=struct;
h_filt=struct;
pVal_filt=struct;
gamma_filt=struct;

fields=fieldnames(data.x_filt);
numChannels=size(data.x_filt.(fields{1}),2);
labels=data.labels;
freqRange=data.freqRange;
fs=data.fs;

config_crit=data.config_crit;

% Main Loop

for i=1:length(fields)
    currCond=fields{i};
    
    fprintf('Beginning condition ''%s''\n',currCond);
    
    numTrials=size(data.x_filt.(currCond),3);
    
    h_filt.(currCond)=nan(numTrials,numChannels);
    pVal_filt.(currCond)=nan(numTrials,numChannels);
    pass_filt.(currCond)=nan(numTrials,1);

    gamma_filt.(currCond)=zeros(numChannels,numChannels,length(freqRange),numTrials);

    for j=1:numTrials
        if data.filt_values.(currCond)(j).decorrelated
            fprintf('%d/%d - ',j,numTrials);

            % Calculate all MVAR models
            [ar_filt.(currCond)(j).mdl,res_filt.(currCond)(j).E,crit_filt.(currCond)(j).(config_crit.crit)]=...
                mvar(squeeze(data.x_filt.(currCond)(:,:,j)),config_crit);

            % Test whiteness
            [~,h_filt.(currCond)(j,:),pVal_filt.(currCond)(j,:)]=...
                test_model(res_filt.(currCond)(j).E,length(data.x_filt.(currCond)(:,:,j)));

            % Calculate DTF
            gamma_filt.(currCond)(:,:,:,j)=dtf(ar_filt.(currCond)(j).mdl,freqRange,fs);
        end
    end
    
    print_whiteness(h_filt.(currCond),pVal_filt.(currCond),labels);
end

save(fullfile(get_root_path,'Files',file),'-append','ar_filt','res_filt','crit_filt',...
    'h_filt','pVal_filt','gamma_filt')

%% Surrogate analysis on filtered data

if isfield(data,'config_surr')
    config_surr=data.config_surr;
    config_surr.signal='decorr';
elseif exist('config_surr','var') == 1
    config_surr.signal='decorr';
else
    config_surr=struct('signal','decorr');
end

surrogate_analysis(file,config_surr);

