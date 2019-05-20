%% plot_file
%
%  Opens up a previously saved file, and plots the relevant figures
%
%   See also: plot_connectivity, plot_criterion, mvar
%

CCC;

%% Load data

FILE='ET_CL_004__2018_06_20__run5__ALLCOND_ALLCOMB.mat';

data=load(FILE);

%% Plot connectivity figure

config=struct('seriesType',1,'plotType','avgerr');

for i=1:length(data.cond_labels)
    currCond=data.cond_labels{i};
    config.figTitle=sprintf('%s, %s, %s - %s: Connectivity',data.PATIENT_ID,data.RECORDING_DATE,data.RUN_ID,currCond);
    plot_connectivity(data.gamma.(currCond),data.x.(currCond),data.freqRange,data.labels,config);
end

return

%% Plot criterion

config=struct; %#ok<UNRCH>

for i=1:length(data.cond_labels)
    currCond=data.cond_labels{i};
    config.hFig=figure;
    
    for j=1:length(data.crit.(currCond))
        config.modelOrder=data.ar.(currCond)(j).mdl.order;
        config.figTitle=sprintf('%s, %s, %s - %s: Criterion',data.PATIENT_ID,data.RECORDING_DATE,data.RUN_ID,currCond);
        plot_criterion(data.crit.(currCond)(j).BIC,config);
    end
end 


