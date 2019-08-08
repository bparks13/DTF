function plot_trial(currCond,trial,file,plots)
%% plot_trial(currCond,trial,file,plots)
%
%  Helper function designed to show the various plots for a single trial, to better
%  isolate differences across trials
%
%   Inputs:
%    - cond: String containing the pertinent condition 
%    - trial: Int defining which trial to look at
%    - file: Optional string containing the absolute path to the file. If not given,
%       uigetfile is called to pick the file
%    - plots: Optional cell array containing the specific plots to show. If empty, all
%       will be shown except for 'crit'
%       'conn': Connectivity
%       'acf': ACF
%       'crit': Criterion
%       'time': Plots the original signal and the estimated signal, as well as the
%           residuals
%       'comp': Compare the original signal with the decorrelated signal
%
%   Outputs:
%    Figure(s) containing all the different plots, including connectivity, ACF, criterion,
%       etc. 
%
%  See also: plot_connectivity, plot_acf, plot_criterion, plot_time_series

bool_plotCrit=false;
bool_plotConn=false;
bool_plotAcf=false;
bool_plotTime=false;
bool_plotCompare=false;

if nargin==2 || isempty(plots)
    bool_plotAll=true;
else
    bool_plotAll=false;
    
    if iscell(plots)
        for i=1:length(plots)
            switch plots{i}
                case 'conn'
                    bool_plotConn=true;
                    break
                    
                case 'acf'
                    bool_plotAcf=true;
                    break
                    
                case 'crit'
                    bool_plotCrit=true;
                    break
                    
                case 'time'
                    bool_plotTime=true;
                    
                case 'comp'
                    bool_plotCompare=true;
            end
        end
    end
end

if nargin == 2
    file=uigetfile(fullfile(get_root_path,'Files','*.mat'));
end

if file ~= 0
    data=load(fullfile(get_root_path,'Files',file));
else
    return
end

conds=fieldnames(data.x);

if ~any(strcmp(currCond,conds))
    error('ERROR: Condition ''%s'' is not in the chosen file',cond);
end

%% Plot connectivity

if bool_plotAll || bool_plotConn
    c_plot=data.config_plot;
    g=data.gamma.(currCond);
    x=data.x.(currCond);
    fr=data.freqRange;
    l=data.labels;

    c_plot.figTitle=sprintf('%s, %s, %s - %s, Trial #%d: Connectivity',...
        data.PATIENT_ID,data.RECORDING_DATE,data.RUN_ID,currCond,trial);
    c_plot.plotType='ind';

    plot_connectivity(g(:,:,:,trial),x(:,:,trial),fr,l,c_plot);
end

%% Plot ACF

if bool_plotAll || bool_plotAcf
    r=data.res.(currCond)(trial).E;

    numChannels=size(data.x.(currCond),2);

    figure('Name',sprintf('ACF - %s, Trial #%d',currCond,trial));

    for i=1:numChannels
        ax=subplot(numChannels,1,i);
        maxLag=round(log(size(data.res.(currCond)(trial).E,1)));
        plot_acf(r(:,i),maxLag,[],ax); ylabel(''); xlabel(''); title(''); 

        if data.h.(currCond)(trial,i)
            ax.Color=[.6 .6 .6];
        end
    end
end

%% Plot Criterion

if bool_plotAll && bool_plotCrit
    c_plot.hFig=figure;
    c_plot.figTitle=sprintf('%s, %s, %s - %s, Trial #%d: Criterion',data.PATIENT_ID,data.RECORDING_DATE,...
        data.RUN_ID,currCond,trial);

    c=data.crit.(currCond).(data.config_crit.crit);

    plot_criterion(c,c_plot);
end

%% Plot Time

if bool_plotAll || bool_plotTime
    numChannels=size(data.x.(currCond),2);
    fs=data.fs;
    t=(0:(length(x)-1))/fs;
    colors=linspecer(3);
    m=data.ar.(currCond)(trial).mdl.order;
    x_hat=data.ar.(currCond)(trial).mdl.x_hat;
    e=data.res.(currCond)(trial).E;

    figure;

    for i=1:numChannels
        subplot(numChannels,1,i);
        plot(t,x(:,i,trial),'Color',colors(1,:)); hold on;
        plot(t(m+1:end),x_hat(:,i),'Color',colors(2,:));
        plot(t(m+1:end),e(:,i),'Color',colors(3,:));
        xlim([t(1) t(end)])
        legend('x','x\_hat','error');
    end
end

%% Plot Compare

if bool_plotAll || bool_plotCompare
    % Add in a function call here to the comparator i am about to make
end

end