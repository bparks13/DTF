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
    if isfield(config,'cond')
        fields=config.cond;
    end
    
    if isfield(config,'method')
        
    end
    
    if isfield(config,'signal')
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

numFields=length(fields);
totalOperations=numFields*numIterations;

for k=1:numFields
    currCond=fields{k};
    
    if bool_original
        x=data.x.(currCond);
    else
        x=data.x_filt.(currCond);
    end
    
    config_crit.orderRange=round(summarize_model_orders(data.ar.(currCond)));
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

%% Either save the file, or return the struct

if nargout==0
    if bool_original
        save(fullfile(get_root_path,'Files',file),'-append','surrogate')
        clear surrogate
    else
        surrogate_filt=surrogate; %#ok<NASGU>
        save(fullfile(get_root_path,'Files',file),'-append','surrogate_filt')
        clear surrogate_filt surrogate
    end
end

end