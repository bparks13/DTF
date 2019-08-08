function plot_connectivity_decorr(file)
%% plot_connectivity_decorr(file)
%
%  Helper function that takes in a file name, and plots all the connectivity figures for
%  the decorrelated data, including any surrogate analysis that was run
%
%   Inputs:
%    - file: String containing the filename, but not the path to the file. If not given,
%       the function calls uigetfile to grab the correct filename
%
%   Outputs:
%    One figure for each condition that is run for the file
%
%  See also: plot_connectivity
%

if nargin==0
    file=uigetfile(fullfile(get_root_path,'Files','*.mat'));
end

data=load(fullfile(get_root_path,'Files',file));

cond_labels=data.cond_labels;
config_plot=data.config_plot;
PATIENT_ID=data.PATIENT_ID;
RECORDING_DATE=data.RECORDING_DATE;
RUN_ID=data.RUN_ID;
x=data.x_filt;
freqRange=data.freqRange;
labels=data.labels;
gamma=data.gamma_filt;

if isfield(data,'surrogate_filt')
    surrogate=data.surrogate_filt;
end

for i=1:length(cond_labels)
    currCond=cond_labels{i};
    config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
    if exist('surrogate','var') == 1
        series=struct('original',x.(currCond),'surrogate',surrogate.(currCond));
        plot_connectivity(gamma.(currCond),series,freqRange,labels,config_plot);
    else
        plot_connectivity(gamma.(currCond),x.(currCond),freqRange,labels,config_plot);
    end
end

end