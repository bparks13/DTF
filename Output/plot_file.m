%% plot_file
%
%  Opens up a previously saved file, and plots the relevant figures
%
%   See also: plot_connectivity, plot_criterion, mvar
%

CCC;
dtf_startup;

%% Load data

% FILE='ET_CL_004__2018_06_20__run5__200Hz__Z_SCORE__BIC_(1).mat';

file=uigetfile(fullfile(get_root_path,'Files','*.mat'));

load(fullfile(get_root_path,'Files',file));

%% Plot connectivity figure - Original

% for i=1:length(cond_labels)
%     currCond=cond_labels{i};
%     config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
%     if exist('surrogate','var') == 1
%         series=struct('original',x.(currCond),'surrogate',surrogate.(currCond));
%         config_plot.surr_params.highlightSignificance=true;
%         plot_connectivity(gamma.(currCond),series,freqRange,labels,config_plot);
%     else
%         plot_connectivity(gamma.(currCond),x.(currCond),freqRange,labels,config_plot);
%     end
% end

%% Plot connectivity figure - Decorrelated

for i=1:length(cond_labels)
    currCond=cond_labels{i};
    config_plot.figTitle=sprintf('%s, %s, %s - %s: Connectivity (Decorrelated)',PATIENT_ID,RECORDING_DATE,RUN_ID,currCond);
    if exist('surrogate_filt','var') == 1
        series=struct('original',x_filt.(currCond),'surrogate',surrogate_filt.(currCond));
        config_plot.surr_params.highlightSignificance=true;
        plot_connectivity(gamma_filt.(currCond),series,freqRange,labels,config_plot);
    else
        plot_connectivity(gamma_filt.(currCond),x_filt.(currCond),freqRange,labels,config_plot);
    end
end

%% Plot criterion

config=struct; 
config.crit=config_crit.crit;

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
    xlim([t(1) t(end)]); xlabel('Time [s]'); ylabel('Z-Score'); 
    legend('x','x\_hat','error');
end

%% Plot the ACF of the residuals

margins=0.04;

for i=1:numConditions
    currCond=cond_labels{i};
    figure('Name',sprintf('ACF - %s',currCond));
    numTrials=size(x.(currCond),3);
    
    for j=1:numTrials
        for k=1:numChannels
            currSubPlot=(j-1)*numChannels+k;
            ax=subplot_tight(numTrials,numChannels,currSubPlot,margins);
            plot_acf(res.(currCond)(j).E(:,k),[],[],ax); ylabel(''); xlabel(''); title(''); 
            ylim([-0.2 0.2]); xlim([1 10]); 
            
            if h.(currCond)(j,k)
                ax.Color=[.6 .6 .6];
            end
        end
    end
    
    drawnow;
end

%% Comparing the PSD of the original signal and the estimated signal

figure;
currCond='Rest';
trialNum=5;
numChannels=size(x_all,2);

for i=1:numChannels
    subplot(numChannels,1,i);
    pxx_sig=pwelch(x.(currCond)(:,i,trialNum),fs,fs/2,1:(fs/2),fs); 
    if isempty(ar.(currCond)(trialNum).mdl.pxx)
       pxx_ar=pwelch(ar.(currCond)(trialNum).mdl.x_hat,fs,fs/2,1:(fs/2),fs); 
    else
        pxx_ar=ar.(currCond)(trialNum).mdl.pxx;
    end
    plot(freqRange,10*log10(pxx_sig(freqRange)),'b'); hold on;
    plot(freqRange,10*log10(pxx_ar(freqRange,i)),'r');
end




