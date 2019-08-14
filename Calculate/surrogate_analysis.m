function [surrogate,distribution]=surrogate_analysis(file,config)
%% [results,distribution]=surrogate_analysis(file,config)
%
%  Performs surrogate analysis on a pre-existing file, creating an empirical distribution
%  of samples and then calculating if the values found for connectivity are significant
%  based on this distribution
%
%   Inputs:
%    - file: String containing the filename to be opened. If this is empty, uigetfile is
%       called so the user can pick the file to open. Note that this is not the absolute
%       path, only the filename itself
%    - config: Optional struct containing additional parameters
%       cond: A cell array containing one or more conditions to focus on, running
%           surrogate on each individual trial with no combining of trials
%       method: String indicating which method to use. All surrogate analysis uses
%           trial-shuffling, but this can be constrained to a single condition ['single']
%           (default) or can shuffle trials from all conditions together ['combine']
%       signal: String defining if the original signal is used ['original', default], or the
%           decorrelated signal ['decorr']
%
%   Outputs:
%    - surrogate: Struct containing all of the values created by the surrogate analysis
%    - distribution: Struct containing the specific trials used in each condition, to
%       check that the surrogate analysis has a uniform distribution
%  
%  See also: mvar, dtf
%

%  Based on the trial-shuffling procedure for surrogate found in
%  10.1016/j.neuroimage.2004.09.036 

method='single';
bool_original=true;
bool_structGiven=false;

%% Load file

if nargin==0 || isempty(file)
    file=uigetfile(fullfile(get_root_path,'Files','*.mat'));
    
    if file == 0
        surrogate=[];
        return
    end
elseif isstruct(file)
    bool_structGiven=true;
end

if ~bool_structGiven
    data=load(fullfile(get_root_path,'Files',file));
else
    data=file;
    
    [~,name,ext] = fileparts(data.newFile);
    file=[name,ext];
end

fs=data.fs;
freqForAnalysis=data.config_plot.freqLims;
freqRange=data.freqRange;
fields=fieldnames(data.x);

if nargin == 2 && isstruct(config)
    if isfield(config,'cond') && ~isempty(config.cond)
        fields=config.cond;
    end
    
    if isfield(config,'method') && ~isempty(config.method)
        method=config.method;
    end
    
    if isfield(config,'signal') && ~isempty(config.signal)
        if strcmp(config.signal,'original')
            bool_original=true;
        elseif strcmp(config.signal,'decorr') && isfield(data,'x_filt')
            bool_original=false;
        end
    end
end

numIterations=1000;

config_crit=struct(...
    'orderSelection',data.config_crit.orderSelection,...
    'crit',data.config_crit.crit,...
    'orderRange',[],...
    'method',data.config_crit.method,...
    'fs',fs,...
    'freqRange',freqForAnalysis,...
    'output',0);

surrogate=struct;
distribution=struct;

%% Surrogate analysis

if strcmp(method,'single')
    numFields=length(fields);
    totalOperations=numFields*numIterations;

    for k=1:numFields
        currCond=fields{k};

        if bool_original
            x=data.x.(currCond);
            ar=data.ar.(currCond);
        else
            x=data.x_filt.(currCond);
            ar=data.ar_filt.(currCond);
        end

        config_crit.orderRange=round(summarize_model_orders(ar));
        numSamples=size(x,1);
        numChannels=size(x,2);
        numTrials=size(x,3);

        tmp_x=zeros(numSamples,numChannels);
        gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

        distribution.(currCond)=nan(numIterations,numChannels);

        for i=1:numIterations
            randomTrials=randperm(numTrials,numChannels);

            for j=1:numChannels
                tmp_x(:,j)=x(:,j,randomTrials(j));
            end

            [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
            gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

            distribution.(currCond)(i,:)=randomTrials;

            if mod((k-1)*numIterations+i,floor(totalOperations/20)) == 0
                fprintf('%d%%\n',floor(((k-1)*numIterations+i)/totalOperations*100));
            end
        end

        surrogate.(currCond)=gamma_dist;
    end
elseif strcmp(method,'combine')
    minLength=inf;
    numTrials=0;
    
    if bool_original
        xVar='x';
    else
        xVar='x_filt';
    end
    
    numChannels=size(data.(xVar).(fields{1}),2);
    
    for i=1:length(fields)
        if minLength > size(data.(xVar).(fields{i}),1)
            minLength=size(data.(xVar).(fields{i}),1);
        end
        
        numTrials=numTrials+size(data.(xVar).(fields{i}),3);
    end
    
    x=nan(minLength,numChannels,numTrials);
    names=cell(numTrials,1);
    count=1;
    
    for i=1:length(fields)
        for j=1:size(data.(xVar).(fields{i}),3)
            x(:,:,count)=data.(xVar).(fields{i})(1:minLength,:,j);
            names{count}=sprintf('%s%d',fields{i},j);
            count=count+1;
        end
    end
    
    tmp_x=zeros(minLength,numChannels);
    gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

    tmp_dist=cell(numIterations,numChannels);
    config_crit.orderRange=round(summarize_model_orders(data.ar));

    for i=1:numIterations
        randomTrials=randperm(numTrials,numChannels);
        
        for j=1:numChannels
            tmp_x(:,j)=x(:,j,randomTrials(j));
            tmp_dist{i,j}=names{randomTrials(j)};
        end
        
        [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
        gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);

        if mod(i,floor(numIterations/20)) == 0
            fprintf('%d%%\n',floor(i/numIterations*100));
        end
    end
    
    for i=1:length(fields)
        surrogate.(fields{i})=gamma_dist;
        distribution.(fields{i})=tmp_dist;
    end
end
%% Either save the file, or return the struct

if nargout==0
    if bool_original
        save(fullfile(get_root_path,'Files',file),'-append','surrogate','distribution')
        clear surrogate distribution
    else
        surrogate_filt=surrogate; %#ok<*NASGU>
        distribution_filt=distribution;
        save(fullfile(get_root_path,'Files',file),'-append','surrogate_filt','distribution_filt')
        clear surrogate_filt surrogate distribution distribution_filt
    end
end

end