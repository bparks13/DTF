function [surrogate]=surrogate_analysis(file)
%% [results]=surrogate_analysis(file)
%
%  Performs surrogate analysis on a pre-existing file, creating an empirical distribution
%  of samples and then calculating if the values found for connectivity are significant
%  based on this distribution
%
%   Inputs:
%    - file: String containing the filename to be opened. If this is empty, uigetfile is
%       called so the user can pick the file to open
%
%   Outputs:
%    - surrogate: Struct containing all of the values created by the surrogate analysis
%       
%
%  Based on the trial-shuffling procedure for surrogate found in
%  10.1016/j.neuroimage.2004.09.036 
%

%% Load file

% FILE='ET_CL_004__2018_06_20__run5__PSD__Z_SCORE.mat';

if nargin==0
    file=uigetfile(fullfile(get_root_path(),'Files','*.mat'));
end

if file~=0
    data=load(fullfile(get_root_path(),'Files',file));
end

fs=data.fs;
freqForAnalysis=data.config_plot.freqLims;
freqRange=data.freqRange;

config_crit=struct('orderSelection',data.config_crit.orderSelection,...
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
    disp(i)
end

surrogate.Rest=gamma_dist;

if nargout==0
    save(fullfile(get_root_path(),'Files',file),'surrogate')
else
    clear results
end

end