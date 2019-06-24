function [surrogate]=surrogate_analysis(file,config)
%% [results]=surrogate_analysis(file,config)
%
%  Performs surrogate analysis on a pre-existing file, creating an empirical distribution
%  of samples and then calculating if the values found for connectivity are significant
%  based on this distribution
%
%   Inputs:
%    - file: String containing the filename to be opened. If this is empty, uigetfile is
%       called so the user can pick the file to open
%    - config: Optional struct containing additional parameters
%       cond: A cell array containing one or more conditions to focus on, running
%           surrogate on each individual trial with no combining of trials
%       method: String indicating which method to use. All surrogate analysis uses
%           trial-shuffling, but this can be constrained to a single condition ['single']
%           (default) or can shuffle trials from all conditions together ['combine']
%
%   Outputs:
%    - surrogate: Struct containing all of the values created by the surrogate analysis
%       
%
%  Based on the trial-shuffling procedure for surrogate found in
%  10.1016/j.neuroimage.2004.09.036 
%

%% Load file

if nargin==0
    file=uigetfile(fullfile(get_root_path,'Files','*.mat'));
end

if file~=0
    data=load(fullfile(get_root_path,'Files',file));
end

fs=data.fs;
freqForAnalysis=data.config_plot.freqLims;
freqRange=data.freqRange;

config_crit=struct(...
    'orderSelection',data.config_crit.orderSelection,...
    'crit',data.config_crit.crit,...
    'orderRange',[],...
    'fs',fs,...
    'freqRange',freqForAnalysis,...
    'output',0);

surrogate=struct;

%% Surrogate on a single condition

x=data.x.Rest;
config_crit.orderRange=round(summarize_model_orders(data.ar.Rest));
numSamples=size(x,1);
numChannels=size(x,2);
numTrials=size(x,3);

numIterations=1000;

tmp_x=zeros(numSamples,numChannels);
gamma_dist=zeros(numChannels,numChannels,length(freqRange),numIterations);

for i=1:numIterations
    randomTrials=randperm(numTrials,numChannels);
    
    for j=1:numChannels
        tmp_x(:,j)=x(:,j,randomTrials(j));
    end
    
    [tmp_mdl,~,~]=mvar(tmp_x,config_crit);
    gamma_dist(:,:,:,i)=dtf(tmp_mdl,freqRange,fs);
    
    if mod(i,floor(numIterations/10)) == 0
        fprintf('%d%%\n',floor(i/numIterations*100));
    end
end

surrogate.Rest=gamma_dist;

if nargout==0
    save(fullfile(get_root_path,'Files',regexprep(file,'.mat','__SURROGATE.mat')),'surrogate')
    clear surrogate
end

end