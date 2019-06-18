%% plot_file
%
%  Opens up a previously saved file, and plots the relevant figures
%
%   See also: plot_connectivity, plot_criterion, mvar
%

CCC;

%% Load data

FILE='ET_CL_004__2018_06_20__run5__DOWNSAMPLE_600Hz.mat';

load(FILE);

%% Plot connectivity figure

for i=1:length(cond_labels)
    currCond=cond_labels{i};
%     config.h=data.h.(currCond);
    config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
    plot_connectivity(gamma.(currCond),x.(currCond),freqRange,labels,config_plot);
end

return

%#ok<*UNRCH>

%% Plot criterion

config=struct; 

for i=1:length(cond_labels)
    currCond=cond_labels{i};
    config.hFig=figure;
    
    for j=1:length(crit.(currCond))
        config.modelOrder=ar.(currCond)(j).mdl.order;
        config.figTitle=sprintf('%s, %s, %s - %s: Criterion',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
        plot_criterion(crit.(currCond)(j).(config_crit.crit),config);
    end
end 

%% Plot original, estimated, and residuals

figure;
numChannels=size(x_all,2);
currCond='Rest';
trialNum=1;
t=(0:(length(x.(currCond))-1))/fs;
colors=linspecer(3);
modelOrder=ar.(currCond)(trialNum).mdl.order;

for i=1:numChannels
    subplot(numChannels,1,i);
    plot(t,x.(currCond)(:,i,trialNum),'Color',colors(1,:)); hold on;
    plot(t(modelOrder+1:end),ar.(currCond)(trialNum).mdl.x_hat(:,i),'Color',colors(2,:));
    plot(t(modelOrder+1:end),res.(currCond)(trialNum).E(:,i),'Color',colors(3,:));
    xlim([t(1) t(end)])
    legend('x','x\_hat','error');
end

%% Comparing the PSD of the original signal and the estimated signal

figure;
currCond='Rest';
trialNum=5;
numChannels=size(x_all,2);

for i=1:numChannels
    subplot(numChannels,1,i);
    pxx_sig=pwelch(x.(currCond)(:,i,trialNum),fs,fs/2,1:(fs/2),fs); 
    plot(freqRange,10*log10(pxx_sig(freqRange))); hold on;
    plot(freqRange,10*log10(ar.(currCond)(trialNum).mdl.pxx(freqRange,i)));
end


